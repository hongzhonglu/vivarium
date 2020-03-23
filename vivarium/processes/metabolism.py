from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process_with_environment,
    plot_simulation_output
)
from vivarium.utils.units import units
from vivarium.utils.cobra_fba import CobraFBA
from vivarium.utils.dict_utils import tuplify_port_dicts, deep_merge
from vivarium.utils.regulation_logic import build_rule
from vivarium.processes.derive_globals import AVOGADRO

# external concentrations lower than exchange threshold are considered depleted
EXCHANGE_THRESHOLD = 1e-6
GLOBALS = ['volume', 'mass', 'mmol_to_counts']


class Metabolism(Process):
    """
    A general metabolism process, which sets the FBA problem based on input configuration data.
    initial_parameters (dict) configures the process with the following keys/values:
        - initial_state (dict) -- the default state, with a dict for internal and external:
            {'external': external_state, 'internal': internal_state}
        - stoichiometry (dict) -- {reaction_id: stoichiometry dict}
        - objective (dict) -- stoichiometry dict to be optimized
        - external_molecules (list) -- the external molecules
        - reversible_reactions (list)
    """
    def __init__(self, initial_parameters={}):
        self.nAvogadro = AVOGADRO

        # initialize FBA
        self.fba = CobraFBA(initial_parameters)
        self.reaction_ids = self.fba.reaction_ids()

        # additional FBA options
        self.constrained_reaction_ids = initial_parameters.get('constrained_reaction_ids', [])
        self.default_upper_bound = initial_parameters.get('default_upper_bound', 1000.0)

        # get regulation functions
        regulation_logic = initial_parameters.get('regulation', {})
        self.regulation = {
            reaction: build_rule(logic) for reaction, logic in regulation_logic.items()}

        # get molecules from fba objective
        self.objective_composition = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if mol_id in self.objective_composition:
                    self.objective_composition[mol_id] += coeff1 * coeff2
                else:
                    self.objective_composition[mol_id] = coeff1 * coeff2

        ## Get initial state from state specified in parameters
        default_state = initial_parameters.get('initial_state', {})
        global_state = default_state.get('global', {})
        internal_state = default_state.get('internal', {})
        external_state = default_state.get('external', {})
        initial_mass = global_state.get('mass', 0)
        mmol_to_counts = global_state.get('mmol_to_counts', 1)
        density = global_state.get('density') * units.g / units.L
        mw = self.fba.molecular_weights

        ## update initial states to match global mass
        ## internal state
        # get initial internal pools based on objective composition, molecular weights, and initial mass
        composition = {mol_id: (-coeff if coeff < 0 else 0)
            for mol_id, coeff in self.objective_composition.items()}
        composition_mass = sum([coeff * mw.get(mol_id, 0.0)
            for mol_id, coeff in composition.items()])
        scale_mass = initial_mass / composition_mass
        updated_internal_state = {mol_id: int(scale_mass * coeff * mmol_to_counts)
            for mol_id, coeff in composition.items()}
        updated_internal_state = deep_merge(dict(updated_internal_state), internal_state)

        ## update global state
        updated_mass = sum([count / mmol_to_counts * mw.get(mol_id, 0.0)
            for mol_id, count in updated_internal_state.items()]) * units.fg
        updated_volume = (updated_mass.to('g') / density)
        updated_mmol_to_counts = (AVOGADRO * updated_volume)

        updated_global_state = {
            'mass': updated_mass.magnitude,
            'volume': updated_volume.to('fL').magnitude,
            'mmol_to_counts': updated_mmol_to_counts.to('L/mmol').magnitude,
            'density': density.magnitude}

        ## external state
        updated_external_state = {state_id: 0.0 for state_id in self.fba.external_molecules}
        updated_external_state.update(self.fba.minimal_external)  # optimal minimal media from fba
        updated_external_state.update(external_state)

        # save initial state
        self.initial_state = {
            'external': updated_external_state,
            'internal': updated_internal_state,
            'reactions': {state_id: 0 for state_id in self.reaction_ids},
            'exchange': {state_id: 0 for state_id in self.fba.external_molecules},
            'flux_bounds': {state_id: self.default_upper_bound
                            for state_id in self.constrained_reaction_ids},
            'global': updated_global_state}

        ## assign ports
        self.internal_state_ids = list(self.objective_composition.keys())
        ports = {
            'external': self.fba.external_molecules,
            'internal': self.internal_state_ids,
            'reactions': self.reaction_ids,
            'exchange': self.fba.external_molecules,
            'flux_bounds': self.constrained_reaction_ids,
            'global': GLOBALS}

        ## parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(ports, parameters)

    def default_settings(self):

        # default emitter keys
        default_emitter_keys = {
            'internal': self.internal_state_ids,
            'external': self.fba.external_molecules,
            'reactions': self.reaction_ids,
            'flux_bounds': self.constrained_reaction_ids,
            'global': ['mass']}

        # schema
        schema = {
            'reactions': {rxn_id: {
                    'updater': 'set',
                    'divide': 'set'}
                for rxn_id in self.reaction_ids},
            'global': {'mass': {
                    'updater': 'accumulate'}}}

        return {
            'state': self.initial_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

    def next_update(self, timestep, states):

        ## get the state
        external_state = states['external']
        constrained_reaction_bounds = states['flux_bounds']  # (units.mmol / units.L / units.s)
        volume = states['global']['volume'] * units.fL
        mmol_to_counts = states['global']['mmol_to_counts'] * units.L / units.mmol

        ## get flux constraints
        # exchange_constraints based on external availability
        exchange_constraints = {mol_id: 0.0
            for mol_id, conc in external_state.items() if conc <= EXCHANGE_THRESHOLD}

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
        regulation_state = {}
        for reaction_id, reg_logic in self.regulation.items():
            regulation_state[reaction_id] = reg_logic(flattened_states)

        ## apply flux constraints
        # first, add exchange constraints
        self.fba.set_exchange_bounds(exchange_constraints)

        # next, add constraints coming from flux_bounds
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
        # calculate added mass from the objective molecules' molecular weights
        objective_count = (objective_exchange * mmol_to_counts).magnitude
        added_mass = 0.0
        internal_state_update = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                internal_state_update[mol_id] = int(-coeff1 * coeff2 * objective_count)

                # added biomass
                mol_mw = self.fba.molecular_weights.get(mol_id, 0.0) * (units.g / units.mol)
                mol_mass = volume * mol_mw.to('g/mmol') * objective_exchange * (units.mmol / units.L)
                added_mass += mol_mass.to('fg').magnitude

        global_update = {'mass': added_mass}

        # convert exchange fluxes to counts
        # TODO -- use derive_counts for exchange
        exchange_deltas = {
            reaction: int((flux * mmol_to_counts).magnitude)
            for reaction, flux in exchange_fluxes.items()}

        all_fluxes = {}
        all_fluxes.update(internal_fluxes)
        all_fluxes.update(exchange_reactions)

        return {
            'exchange': exchange_deltas,
            'internal': internal_state_update,
            'reactions': all_fluxes,
            'global': global_update}



# tests and analyses
def plot_exchanges(timeseries, sim_config, out_dir):
    # plot focused on exchanges

    nAvogadro = AVOGADRO
    env_volume = sim_config['environment_volume']
    timeline = sim_config['timeline']  # TODO -- add tickmarks for timeline events
    external_ts = timeseries['external']
    internal_ts = timeseries['internal']
    global_ts = timeseries['global']

    # pull volume and mass out from internal
    volume = global_ts.pop('volume') * units.fL
    mass = global_ts.pop('mass')

    # conversion factor
    mmol_to_counts = [nAvogadro.to('1/mmol') * vol.to('L') for vol in volume]

    # plot results
    cols = 1
    rows = 3
    plt.figure(figsize=(cols * 6, rows * 1.5))

    # define subplots
    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)

    # plot external state
    for mol_id, series in external_ts.items():
        ax1.plot(series, label=mol_id)
    ax1.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=2)
    ax1.title.set_text('environment: {} (L)'.format(env_volume))
    ax1.set_ylabel('concentrations')
    ax1.set_yscale('log')

    # plot internal counts
    for mol_id, counts_series in internal_ts.items():
        # conc_series = [(count / conversion).to('mmol/L').magnitude
        #    for count, conversion in zip(counts_series, mmol_to_counts)]
        ax2.plot(counts_series, label=mol_id)

    # # plot internal concentrations
    # for mol_id, counts_series in internal_ts.items():
    #     conc_series = [(count / conversion).to('mmol/L').magnitude
    #        for count, conversion in zip(counts_series, mmol_to_counts)]
    #     ax2.plot(conc_series, label=mol_id)

    # ax2.legend(loc='center left', bbox_to_anchor=(1.6, 0.5), ncol=3)
    # ax2.title.set_text('internal metabolites')
    # ax2.set_ylabel('conc (mM)')

    ax2.legend(loc='center left', bbox_to_anchor=(1.6, 0.5), ncol=3)
    ax2.title.set_text('internal metabolites')
    ax2.set_ylabel('delta counts')

    # plot mass
    ax3.plot(mass, label='mass')
    ax3.set_ylabel('mass (fg)')

    # adjust axes
    for axis in [ax1, ax2, ax3]:
        axis.spines['right'].set_visible(False)
        axis.spines['top'].set_visible(False)

    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xlabel('time (s)', fontsize=12)

    # save figure
    fig_path = os.path.join(out_dir, 'exchanges')
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')

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

    settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 10}

    timeseries = simulate_process_with_environment(metabolism, settings)
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

# configs
def get_initial_global_state():
    # initial state
    mass = 1339 * units.fg
    density = 1100 * units.g / units.L
    volume = mass.to('g') / density
    mmol_to_counts = (AVOGADRO * volume).to('L/mmol')

    globals = {
        'density': density.magnitude,
        'mass': mass.magnitude,  # fg
        'volume': volume.to('fL').magnitude,
        'mmol_to_counts': mmol_to_counts.magnitude}

    return {
        'global': globals}

def get_e_coli_core_config():
    metabolism_file = os.path.join('models', 'e_coli_core.json')
    initial_state = get_initial_global_state()
    return {
        'model_path': metabolism_file,
        'initial_state': initial_state}

def get_iAF1260b_config():
    metabolism_file = os.path.join('models', 'iAF1260b.json')
    initial_state = get_initial_global_state()
    return {
        'model_path': metabolism_file,
        'initial_state': initial_state}

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
        'A': -0.02,
        'D': 0.01,
        'E': 0.01,
        'F': -0.005,
        'H': -0.005,
        'O2': -0.1}

    # mass = 1339 * units.fg
    # density = 1100 * units.g/units.L
    # volume = mass.to('g') / density
    global_state = get_initial_global_state()
    initial_state = {
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0}}
    initial_state.update(global_state)

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
        'molecular_weights': molecular_weights}

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

# toy functions
def make_kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def make_random_rate(mu=0, sigma=1):
    def rate(state):
        return np.random.normal(loc=mu, scale=sigma)
    return rate

def toy_transport():
    # process-like function for transport kinetics, used by simulate_metabolism
    transport_kinetics = {
        "EX_A": make_kinetic_rate("A", -1e-1, 5),  # A import
    }
    return transport_kinetics

# tests
def test_BiGG_metabolism(settings={}):
    metabolism_config = get_iAF1260b_config()
    # metabolism_config = get_e_coli_core_config()
    metabolism = Metabolism(metabolism_config)

    # simulate metabolism
    sim_settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'timeline': [(10, {})]}
    sim_settings.update(settings)

    saved_data = simulate_process_with_environment(metabolism, sim_settings)

    return saved_data

def test_toy_metabolism(out_dir='out'):

    regulation_logic = {
        'R4': 'if (external, O2) > 0.1 and not (external, F) < 0.1'}

    toy_config = get_toy_configuration()
    transport = toy_transport()

    toy_config['constrained_reaction_ids'] = list(transport.keys())
    toy_config['regulation'] = regulation_logic
    toy_metabolism = Metabolism(toy_config)

    # simulate toy model
    timeline = [
        (10, {'external': {
            'A': 1}}),
        (20, {'external': {
            'F': 0}}),
        (30, {})]

    settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'environment_volume': 1e-14,  # L
        'timestep': 1,
        'timeline': timeline}

    saved_data = simulate_process_with_environment(toy_metabolism, settings)
    return saved_data

def make_network():
    # configure toy model
    toy_config = get_toy_configuration()
    toy_metabolism = Metabolism(toy_config)

    # make flux network from toy model
    network_config = {
        'process': toy_metabolism,
        'total_time': 10,
        'environment_volume': 5e-13}
    save_network(network_config, out_dir)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_dir_BiGG = os.path.join('out', 'tests', 'metabolism_BiGG')
    if not os.path.exists(out_dir_BiGG):
        os.makedirs(out_dir_BiGG)

    # run BiGG metabolism
    timeline = [(2520, {})]
    sim_settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'environment_volume': 1e-13,  # L
        'timestep': 1,
        'timeline': timeline}

    plot_settings = {
        'max_rows': 30,
        'remove_zeros': True,
        'skip_ports': ['exchange', 'flux_bounds', 'reactions'],
        'overlay': {
            'reactions': 'flux_bounds'}}

    timeseries = test_BiGG_metabolism(sim_settings) # 2520 sec (42 min) is the expected doubling time in minimal media
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    plot_simulation_output(timeseries, plot_settings, out_dir_BiGG)
    plot_exchanges(timeseries, sim_settings, out_dir_BiGG)

    # run toy model
    plot_settings = {
        'skip_ports': ['exchange'],
        'overlay': {
            'reactions': 'flux_bounds'}}

    timeseries = test_toy_metabolism()
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # make_network()
