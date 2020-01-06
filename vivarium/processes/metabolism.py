from __future__ import absolute_import, division, print_function

import os
from scipy import constants
import numpy as np
import copy

from vivarium.actor.process import Process, deep_merge
from vivarium.utils.units import units
from vivarium.utils.cobra_fba import CobraFBA
from vivarium.actor.process import convert_to_timeseries, plot_simulation_output
import vivarium.utils.regulation_logic as rl
from vivarium.utils.dict_utils import flatten_role_dicts



class Metabolism(Process):
    '''
    A general metabolism process, which sets the FBA problem based on input configuration data.
    initial_parameters (dict) configures the process with the following keys/values:
        - initial_state (dict) -- the default state, with a dict for internal and external:
            {'external': external_state, 'internal': internal_state}
        - stoichiometry (dict) -- {reaction_id: stoichiometry dict}
        - objective (dict) -- stoichiometry dict to be optimized
        - external_molecules (list) -- the external molecules
        - reversible_reactions (list)
    '''
    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1/units.mol
        self.density = 1100 * units.g/units.L

        # initialize fba
        self.fba = CobraFBA(initial_parameters)
        self.reaction_ids = self.fba.reaction_ids()  #list(self.fba.stoichiometry.keys())

        # additional options
        self.constrained_reaction_ids = initial_parameters.get('constrained_reaction_ids', [])
        self.initial_state = initial_parameters.get('initial_state', {})
        self.default_upper_bound = initial_parameters.get('default_upper_bound', 1000.0)
        self.regulation = initial_parameters.get('regulation', {})

        # get molecules in objective
        self.objective_molecules = []
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                self.objective_molecules.append(mol_id)

        # assign internal and external roles
        self.internal_state_ids = self.objective_molecules + ['volume', 'mass']
        roles = {
            'external': self.fba.external_molecules,
            'internal': self.internal_state_ids,
            'reactions': self.reaction_ids,
            'exchange': self.fba.external_molecules,
            'flux_bounds': self.constrained_reaction_ids}

        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0.0 for state_id in self.internal_state_ids}
        external = {state_id: 0.0 for state_id in self.fba.external_molecules}
        external.update(self.fba.minimal_external)
        external.update(self.initial_state.get('external', {}))

        default_state = {
            'external': external,
            'internal': deep_merge(dict(internal), self.initial_state.get('internal', {})),
            'reactions': {state_id: 0 for state_id in self.reaction_ids},
            'exchange': {state_id: 0 for state_id in self.fba.external_molecules},
            'flux_bounds': {state_id: self.default_upper_bound
                for state_id in self.constrained_reaction_ids}
            }

        # default emitter keys
        default_emitter_keys = {
            'internal': self.objective_molecules,
            'external': self.fba.external_molecules,
            'reactions': self.reaction_ids,
            'exchange': self.fba.external_molecules,
        }

        # default updaters
        set_internal_states = {state_id: 'set'
            for state_id in self.reaction_ids + ['volume']}
        accumulate_internal_states = {state_id: 'accumulate'
            for state_id in self.objective_molecules + ['mass']}
        default_updaters = {
            'internal': deep_merge(dict(set_internal_states), accumulate_internal_states),
            'external': {mol_id: 'accumulate' for mol_id in self.fba.external_molecules},
            'reactions': {rxn_id: 'set' for rxn_id in self.reaction_ids},
            'exchange': {rxn_id: 'set' for rxn_id in self.fba.external_molecules},
            'flux_bounds': {rxn_id: 'set' for rxn_id in self.constrained_reaction_ids},
            }

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

    def next_update(self, timestep, states):

        ## get the state
        internal_state = states['internal']
        external_state = states['external']  # TODO -- constrain metabolism by external state
        constrained_reaction_bounds = states['flux_bounds']  # (units.mmol / units.L / units.s)
        mass = internal_state['mass'] * units.fg
        volume = mass.to('g') / self.density

        ## conversion factors
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume

        ## get flux constraints
        # exchange_constraints based on external availability
        exchange_constraints = {mol_id: 0.0
            for mol_id, conc in external_state.items() if conc <= 0.0}

        # availibility is boolean state of all molecules present
        # regulation_state determines state of regulated reactions (True/False)
        availability = flatten_role_dicts({role: states[role]
            for role in ('internal', 'external')})
        availability = {state: value > 0
            for state, value in availability.items()}
        regulation_state = {rxn_id: logic(availability)
            for rxn_id, logic in self.regulation.items()}

        ## apply flux constraints
        # first, add exchange constraints
        self.fba.set_exchange_bounds(exchange_constraints)

        # next, add constraints coming from flux_bounds (typically transport)
        # to constrain exchange fluxes, add the suffix 'EX_' to the external molecule ID
        self.fba.constrain_flux(constrained_reaction_bounds)

        # finally, turn reactions on/off based on regulation
        self.fba.regulate_flux(regulation_state)

        ## solve the fba problem
        objective_exchange = self.fba.optimize() * timestep  # (units.mmol / units.L / units.s)
        exchange_reactions = self.fba.read_exchange_reactions()
        exchange_fluxes = self.fba.read_exchange_fluxes()  # (units.mmol / units.L / units.s)
        internal_fluxes = self.fba.read_internal_fluxes()  # (units.mmol / units.L / units.s)

        # timestep dependence
        exchange_fluxes.update((mol_id, flux * timestep) for mol_id, flux in exchange_fluxes.items())
        internal_fluxes.update((mol_id, flux * timestep) for mol_id, flux in internal_fluxes.items())

        # update internal counts from objective flux
        # calculate the new mass from the objective molecules' molecular weights
        objective_count = (objective_exchange * mmol_to_count).magnitude
        added_mass = 0.0
        internal_state_update = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                internal_state_update[mol_id] = int(-coeff1 * coeff2 * objective_count)

                # added biomass
                mol_mw = self.fba.molecular_weights.get(mol_id, 0.0) * (units.g / units.mol)
                mol_mass = volume * mol_mw.to('g/mmol') * objective_exchange * (units.mmol / units.L)
                added_mass += mol_mass.to('fg').magnitude # to fg

        internal_state_update.update({'mass': added_mass})

        # convert exchange fluxes to counts with mmol_to_count
        exchange_deltas = {
            reaction: int((flux * mmol_to_count).magnitude)
            for reaction, flux in exchange_fluxes.items()}

        return {
            'exchange': exchange_deltas,
            'internal': internal_state_update,
            'reactions': {**internal_fluxes, **exchange_reactions},
        }



# tests and analyses
def simulate_metabolism(config):

    # set the simulation
    metabolism = config['process']
    total_time = config.get('total_time', 3600)
    transport_kinetics = config.get('transport_kinetics', {})
    env_volume = config.get('environment_volume', 1e-12) * units.L
    timeline = config.get('timeline', [(total_time, {})])
    end_time = timeline[-1][0]

    # get initial state and parameters
    settings = metabolism.default_settings()
    state = settings['state']
    internal_updaters = settings['updaters']['internal']
    density = metabolism.density
    nAvogadro = metabolism.nAvogadro

    # initialize saved data
    saved_state = {}

    ## run simulation
    time = 0
    timestep = 1  # sec
    saved_state[time] = state
    while time < end_time:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.items():
                    state[key].update(change)
                timeline.pop(0)

        # set flux bounds from transport kinetics
        flux_bounds = {}
        for rxn_id, rate_law in transport_kinetics.items():
            flux = rate_law(state['external'])
            flux_bounds[rxn_id] = flux
        state['flux_bounds'].update(flux_bounds)

        # get update
        update = metabolism.next_update(timestep, state)

        # reactions are set as is
        state['reactions'].update(update['reactions'])

        # apply internal update
        for state_id, state_update in update['internal'].items():
            if internal_updaters.get(state_id, 'set') is 'accumulate':
                state['internal'][state_id] += state_update
            else:
                state['internal'][state_id] = state_update

        # update volume
        new_mass = state['internal']['mass'] * units.fg
        new_volume = new_mass.to('g') / density
        state['internal']['volume'] = new_volume.to('L').magnitude

        # apply external update -- use exchange without growth rate
        mmol_to_count = (nAvogadro.to('1/mmol') * env_volume).to('L/mmol').magnitude
        for mol_id, exchange_counts in update['exchange'].items():
            exchange_conc = exchange_counts / mmol_to_count  # TODO -- per second?
            state['external'][mol_id] += exchange_conc
            state['exchange'][mol_id] = 0.0 # reset exchange
            if state['external'][mol_id] < 0.0:  # this shouldn't be needed
                state['external'][mol_id] = 0.0

        # save state
        saved_state[time] = copy.deepcopy(state)

    return saved_state

def save_network(config, out_dir='out'):
    # TODO -- make this function into an analysis
    import math
    from vivarium.utils.make_network import make_network, save_network

    # initialize the process
    metabolism = config['process']
    stoichiometry = metabolism.fba.stoichiometry
    reaction_ids = list(stoichiometry.keys())
    external_mol_ids = metabolism.fba.external_molecules
    objective = metabolism.fba.objective

    data = simulate_metabolism(config)
    timeseries = convert_to_timeseries(data)
    reactions =  timeseries['reactions']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        flux = abs(np.mean(reactions[rxn_id][1:]))
        reaction_fluxes[rxn_id] = math.log(1000 * flux + 1.1)

    # define node type
    node_types = {rxn_id: 'reaction' for rxn_id in reaction_ids}
    node_types.update({mol_id: 'external_mol' for mol_id in external_mol_ids})
    node_types.update({rxn_id: 'objective' for rxn_id in objective.keys()})
    info = {
        'node_types': node_types,
        'reaction_fluxes': reaction_fluxes}

    nodes, edges = make_network(stoichiometry, info)
    save_network(nodes, edges, out_dir)

def get_toy_configuration():
    stoichiometry = {
        'R1': {'A': -1, 'ATP': -1, 'B': 1},
        'R2a': {'B': -1, 'ATP': 2, 'NADH': 2, 'C': 1},
        'R2b': {'C': -1, 'ATP': -2, 'NADH': -2, 'B': 1},
        'R3': {'B': -1, 'F': 1},
        'R4': {'C': -1, 'G': 1},
        'R5': {'G': -1, 'C': 0.8, 'NADH': 2},
        'R6': {'C': -1, 'ATP': 2, 'D': 3},
        'R7': {'C': -1, 'NADH': -4, 'E': 3},
        'R8a': {'G': -1, 'ATP': -1, 'NADH': -2, 'H': 1},
        'R8b': {'G': 1, 'ATP': 1, 'NADH': 2, 'H': -1},
        'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1},
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2']

    objective = {'v_biomass': 1.0}

    reversible = ['R6', 'R7', 'Rres']  # stoichiometry.keys()

    default_reaction_bounds = 1000.0

    exchange_bounds = {
        'A': 0.02,
        'D': -0.01,
        'E': -0.01,
        'F': 0.005,
        'H': 0.005,
        'O2': 0.1}

    mass = 1339 * units.fg
    density = 1100 * units.g/units.L
    volume = mass.to('g') / density
    initial_state = {
        'internal': {
            'mass': mass.magnitude,  # fg
            'volume': volume.magnitude},  # L
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0,
        }}

    # molecular weight units are (units.g / units.mol)
    molecular_weights = {
        'A': 500.0,
        'B': 500.0,
        'C': 500.0,
        'D': 500.0,
        'E': 500.0,
        'F': 50000.0,
        'H': 1.00794,
        'O2': 31.9988,
        'ATP': 507.181,
        'NADH': 664.425}

    config = {
        'stoichiometry': stoichiometry,
        'reversible': reversible,
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state,
        'exchange_bounds': exchange_bounds,
        'default_upper_bound': default_reaction_bounds,
        'molecular_weights': molecular_weights,
        }

    return config

def test_config(get_config=get_toy_configuration):
    # configure metabolism process
    config = get_config()
    metabolism = Metabolism(config)

    print('MODEL: {}'.format(metabolism.fba.model))
    print('REACTIONS: {}'.format(metabolism.fba.model.reactions))
    print('METABOLITES: {}'.format(metabolism.fba.model.metabolites))
    print('GENES: {}'.format(metabolism.fba.model.genes))
    print('COMPARTMENTS: {}'.format(metabolism.fba.model.compartments))
    print('SOLVER: {}'.format(metabolism.fba.model.solver))
    print('EXPRESSION: {}'.format(metabolism.fba.model.objective.expression))

    print(metabolism.fba.optimize())
    print(metabolism.fba.model.summary())
    print('internal: {}'.format(metabolism.fba.internal_reactions()))
    print('external: {}'.format(metabolism.fba.external_reactions()))
    print(metabolism.fba.reaction_ids())
    print(metabolism.fba.get_reactions())
    print(metabolism.fba.get_reaction_bounds())
    print(metabolism.fba.read_exchange_fluxes())

def kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def random_rate(mu=0, sigma=1):
    def rate(state):
        return np.random.normal(loc=mu, scale=sigma)
    return rate



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # toy functions
    def toy_transport():
        # process-like function for transport kinetics
        transport_kinetics = {
            "EX_A": kinetic_rate("A", -1e-1, 5),  # A import
            # "R1": kinetic_rate("A", 2.5e-2, 5),  # A import
            # "EX_E": random_rate(1e-2, 1e-3)  # (0, 1e-1),   #  E exchange, requires 'EX_' prefix
        }
        return transport_kinetics

    def toy_regulation():
        regulation = {
            'R4': rl.build_rule('IF (O2_external) and not (F_external)'),
        }
        return regulation

    # configure toy model
    toy_config = get_toy_configuration()
    transport = toy_transport()
    regulation = toy_regulation()
    toy_config['constrained_reaction_ids'] = list(transport.keys())
    toy_config['regulation'] = regulation

    toy_metabolism = Metabolism(toy_config)

    # simulate toy model
    timeline = [
        (300, {'external': {
            'A': 1}
        }),
        (600, {'external': {
            'F': 0}
        }),
        (2500, {})]

    simulation_config = {
        'process': toy_metabolism,
        'timeline': timeline,
        'transport_kinetics': transport,
        'environment_volume': 5e-13}

    plot_settings = {
        'skip_roles': ['exchange'],
        'overlay': {
            'reactions': 'flux_bounds'}}

    saved_data = simulate_metabolism(simulation_config)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # make flux network from toy model
    network_config = {
        'process': toy_metabolism,
        'total_time': 10,
        'environment_volume': 5e-13}
    save_network(network_config, out_dir)
