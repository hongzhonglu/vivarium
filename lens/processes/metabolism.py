from __future__ import absolute_import, division, print_function

import os
from scipy import constants
import numpy as np

from lens.actor.process import Process, deep_merge
from lens.utils.units import units
from lens.utils.cobra_fba import CobraFBA


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

        # get FBA configuration
        self.stoichiometry = initial_parameters['stoichiometry']
        self.objective = initial_parameters['objective']
        self.external_molecules = initial_parameters['external_molecules']
        self.reaction_ids = self.stoichiometry.keys()

        # additional options
        self.initial_state = initial_parameters.get('initial_state', {})
        self.molecular_weights = initial_parameters.get('molecular_weights', {})
        reversible = initial_parameters.get('reversible', [])
        flux_bounds = initial_parameters.get('flux_bounds', {})
        exchange_bounds = initial_parameters.get('exchange_bounds', {})
        self.default_upper_bound = initial_parameters.get('default_upper_bound', 1000.0)

        # initialize fba
        self.fba = CobraFBA(dict(
            stoichiometry=self.stoichiometry,
            reversible=reversible,
            external_molecules=self.external_molecules,
            objective=self.objective,
            initial_state=self.initial_state,
            flux_bounds = flux_bounds,
            default_upper_bound=self.default_upper_bound))

        # print(self.fba.get_reaction_bounds())

        # set bounds on exchange fluxes
        self.fba.constrain_exchange_flux(exchange_bounds)

        # get molecules in objective
        self.objective_molecules = []
        for reaction_id, coeff1 in self.objective.items():
            for mol_id, coeff2 in self.stoichiometry[reaction_id].items():
                self.objective_molecules.append(mol_id)

        # assign internal and external roles
        self.internal_state_ids = self.objective_molecules + ['volume', 'mass']
        roles = {
            'external': self.external_molecules,
            'internal': self.internal_state_ids,
            'reactions': self.reaction_ids,
            'exchange': self.external_molecules,
            'flux_bounds': self.reaction_ids}

        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0.0 for state_id in self.internal_state_ids}
        external = {state_id: 10.0 for state_id in self.external_molecules}
        default_state = {
            'external':  deep_merge(dict(external), self.initial_state.get('external', {})),
            'internal': deep_merge(dict(internal), self.initial_state.get('internal', {})),
            'reactions': {state_id: 0 for state_id in self.reaction_ids},
            'exchange': {state_id: 0 for state_id in self.external_molecules},
            'flux_bounds': {state_id: self.default_upper_bound for state_id in self.reaction_ids}
            }

        # default emitter keys
        default_emitter_keys = {
            'internal': self.objective_molecules,
            'external': self.fba.external_molecules,
            'reactions': self.reaction_ids,
            'exchange': self.external_molecules,
        }

        # default updaters
        set_internal_states = {state_id: 'set' for state_id in self.reaction_ids + ['volume']}
        accumulate_internal_states = {state_id: 'accumulate' for state_id in self.objective_molecules + ['mass']}
        default_updaters = {
            'internal': deep_merge(dict(set_internal_states), accumulate_internal_states),
            'external': {mol_id: 'accumulate' for mol_id in self.external_molecules},
            'reactions': {rxn_id: 'set' for rxn_id in self.reaction_ids},
            'exchange': {rxn_id: 'set' for rxn_id in self.external_molecules},
            'flux_bounds': {rxn_id: 'set' for rxn_id in self.reaction_ids},
            }

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

    def next_update(self, timestep, states):

        internal_state = states['internal']
        external_state = states['external']
        reaction_flux_bounds = states['flux_bounds']
        mass = internal_state['mass'] * units.fg
        volume = mass.to('g') / self.density

        # conversion factors
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume

        # set flux constraints.
        self.fba.constrain_flux(reaction_flux_bounds)

        # solve the fba problem
        objective_exchange = self.fba.optimize()  # (units.mmol / units.L)
        exchange_fluxes = self.fba.read_exchange_fluxes()  # (units.mmol / units.L)
        internal_fluxes = self.fba.read_internal_fluxes()  # (units.mmol / units.L)

        # update internal counts from objective flux
        # calculate the new mass from the objective molecules' molecular weights
        objective_count = (objective_exchange * mmol_to_count).magnitude
        added_mass = 0.0
        internal_state_update = {}
        for reaction_id, coeff1 in self.objective.items():
            for mol_id, coeff2 in self.stoichiometry[reaction_id].items():
                internal_state_update[mol_id] = int(-coeff1 * coeff2 * objective_count)

                # added mass
                mol_mw = self.molecular_weights.get(mol_id, 0.0) * (units.g / units.mol)
                mol_mass = volume * mol_mw.to('g/mmol') * objective_exchange * (units.mmol / units.L)
                added_mass += mol_mass.to('fg').magnitude # to fg

        internal_state_update.update({'mass': added_mass})

        # convert exchange fluxes to counts with mmol_to_count
        environment_deltas = {
            reaction: int((flux * mmol_to_count).magnitude)
            for reaction, flux in exchange_fluxes.items()}

        # return update
        return {
            'internal': internal_state_update,
            'external': environment_deltas,
            'reactions': internal_fluxes,
            'exchange': exchange_fluxes,
        }


# tests and analyses
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

    reversible = stoichiometry.keys()  # ['R6', 'R7', 'Rres']

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
        'F': 500.0,
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

def toy_transport_kinetics():
    # transport kinetics
    transport_kinetics = {
        "R1": kinetic_rate("A", 1e-1, 5),   # A import
        "R3": kinetic_rate("F", 1e-1, 5),  # F export
        # "R6": kinetic_rate("D", 1e-3, 12),  # D export
        # "R7": kinetic_rate("E", 1e-3, km),  # E export
        "R8a": kinetic_rate("H", 1e-1, 5), # H export
        # "R8b": kinetic_rate("H", 1e-5, 0.1),  # H import
        # "Rres": kinetic_rate("O2", 4e-1, 200), # O2 import
    }

    return transport_kinetics

def simulate_metabolism(metabolism, config):

    # set the simulation
    total_time = config.get('total_time', 3600)
    transport_kinetics = config.get('transport_kinetics', {})
    env_volume = config.get('environment_volume', 1e-12)  * units.L

    # get initial state and parameters
    settings = metabolism.default_settings()
    state = settings['state']
    internal_updaters = settings['updaters']['internal']
    density = metabolism.density
    nAvogadro = metabolism.nAvogadro
    reaction_ids = metabolism.reaction_ids
    exchange_ids = metabolism.external_molecules

    saved_data = {
        'internal': {state_id: [value] for state_id, value in state['internal'].items()},
        'external': {state_id: [value] for state_id, value in state['external'].items()},
        'reactions': {rxn_id: [0] for rxn_id in reaction_ids},
        'exchange': {rxn_id: [0] for rxn_id in exchange_ids},
        'flux_bounds': {rxn_id: [0] for rxn_id, rxn_fun in transport_kinetics.items()},
        'time': [0]}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < total_time:

        time += timestep

        # set flux bounds from transport kinetics
        flux_bounds = {}
        for rxn_id, rate_law in transport_kinetics.items():
            flux = rate_law(state['external'])
            flux_bounds[rxn_id] = flux
        state['flux_bounds'].update(flux_bounds)

        # get update
        update = metabolism.next_update(timestep, state)

        # exchanges and reactions are set as is
        state['reactions'] = update['reactions']
        state['exchange'] = update['exchange']

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
        for mol_id, exchange_counts in update['external'].items():
            exchange_conc = exchange_counts / mmol_to_count  # TODO -- per second?
            state['external'][mol_id] += exchange_conc
            if state['external'][mol_id] < 0.0:  # this shouldn't be needed
                state['external'][mol_id] = 0.0

        # save state
        saved_data['time'].append(time)
        for role in ['internal', 'external', 'reactions', 'exchange', 'flux_bounds']:
            for state_id, value in state[role].items():
                saved_data[role][state_id].append(value)

    return saved_data

def plot_output(saved_state, out_dir='out', filename='metabolism'):
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import math

    skip_keys = ['time', 'flux_bounds']

    # remove series with all zeros
    zero_state = []
    for key1 in saved_state.iterkeys():
        if key1 not in skip_keys:
            for key2, series in saved_state[key1].items():
                if all(v == 0 for v in series):
                    zero_state.append((key1, key2))
    for (key1, key2) in zero_state:
        del saved_state[key1][key2]

    max_rows = 25
    data_keys = [key for key in saved_state.keys() if key not in skip_keys]
    flux_bounds_data = saved_state.get('flux_bounds')
    time_vec = [t/60./60. for t in saved_state['time']]  # convert to hours
    n_zeros = len(zero_state)
    n_data = [len(saved_state[key].keys()) for key in data_keys]

    # limit number of rows to max_rows
    n_cols = len(n_data)
    n_rows_base = max(n_data)
    for role in n_data:
        new_rows = role / max_rows
        if new_rows > 1:
            n_cols += int(new_rows)
            n_rows_base = max_rows
    n_rows = n_rows_base + int(math.ceil(n_zeros/20.0))

    # make figure
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 2))
    plt.rcParams.update({'font.size': 12})
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.4, hspace=1.5)

    row_idx = 0
    col_idx = 0
    for key in data_keys:
        for state_id, series in sorted(saved_state[key].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])
            ax.plot(time_vec, series)

            # if state has target, plot series in red
            if state_id in flux_bounds_data.keys():
                target_series = flux_bounds_data[state_id]
                ax.plot(time_vec, target_series, 'r', label='bound')
                ax.legend()

            ax.title.set_text(str(key) + ': ' + state_id)
            ax.ticklabel_format(style='sci',axis='y')
            ax.set_xlabel('time (hr)')

            row_idx += 1

            # limit number of rows
            if row_idx > max_rows:
                col_idx += 1
                row_idx = 0

        col_idx += 1
        row_idx = 0

    # # additional data as text
    # if zero_state:
    #     zeros = ['{}[{}]'.format(state, role) for (role, state) in zero_state]
    #     zeros_text = 'states with all zeros: {}'.format(zeros)
    #     ax = fig.add_subplot(grid[n_rows_base:, :])
    #     ax.text(0.01, 0.99, zeros_text, ha='left', va='top', wrap=True)
    #     ax.axis('off')

    # make figure output directory and save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.5, hspace=0.9)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')

def save_network(metabolism, total_time=60, out_dir='out'):
    # TODO -- make this function into an analysis
    import math
    from lens.utils.make_network import make_network, save_network

    # initialize the process
    stoichiometry = metabolism.stoichiometry
    reaction_ids = stoichiometry.keys()
    external_mol_ids = metabolism.external_molecules
    objective = metabolism.objective

    # run test to get simulation output
    simulation_config = {'total_time': total_time}
    saved_state = simulate_metabolism(metabolism, simulation_config)
    reactions = saved_state['reactions']

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

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## set up metabolism with a toy configuration
    get_config = get_toy_configuration
    metabolism = Metabolism(get_config())

    ## test toy model
    test_config(get_config)

    ## simulate toy model
    simulation_config = {
        'total_time': 7200,
        'transport_kinetics': toy_transport_kinetics(),
        'environment_volume': 5e-15}
    saved_data = simulate_metabolism(metabolism, simulation_config)
    plot_output(saved_data, out_dir)

    ## make flux network from toy model
    save_network(metabolism, 10, out_dir)
