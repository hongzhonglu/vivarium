from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process_with_environment,
    save_timeseries,
    flatten_timeseries,
    load_timeseries,
    assert_timeseries_close,
    REFERENCE_DATA_DIR,
    TEST_OUT_DIR,
    plot_simulation_output
)
from vivarium.utils.make_network import (
    get_reactions,
    make_network,
    save_network
)
from vivarium.data.synonyms import get_synonym
from vivarium.utils.units import units
from vivarium.utils.cobra_fba import CobraFBA
from vivarium.utils.dict_utils import tuplify_port_dicts, deep_merge
from vivarium.utils.regulation_logic import build_rule
from vivarium.processes.derive_globals import AVOGADRO


NAME = 'metabolism_process'
GLOBALS = ['volume', 'mass', 'mmol_to_counts']



def get_fg_from_counts(counts_dict, mw):
    composition_mass = sum([
        coeff / AVOGADRO * mw.get(mol_id, 0.0) * (units.g / units.mol)
        for mol_id, coeff in counts_dict.items()])  # g
    return composition_mass.to('fg')



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

    defaults = {
        'constrained_reaction_ids': [],
        'default_upper_bound': 1000.0,
        'regulation': {},
        'initial_state': {},
        'exchange_threshold': 1e-6, # external concs lower than exchange_threshold are considered depleted
        'initial_mass': 1339,  # fg
    }

    def __init__(self, initial_parameters={}):
        self.nAvogadro = AVOGADRO

        # initialize FBA
        self.fba = CobraFBA(initial_parameters)
        self.reaction_ids = self.fba.reaction_ids()
        self.exchange_threshold = self.defaults['exchange_threshold']

        # additional FBA options
        self.constrained_reaction_ids = initial_parameters.get(
            'constrained_reaction_ids', self.defaults['constrained_reaction_ids'])
        self.default_upper_bound = initial_parameters.get(
            'default_upper_bound', self.defaults['default_upper_bound'])

        # get regulation functions
        regulation_logic = initial_parameters.get(
            'regulation', self.defaults['regulation'])
        self.regulation = {
            reaction: build_rule(logic)
            for reaction, logic in regulation_logic.items()}

        # get molecules from fba objective
        self.objective_composition = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if mol_id in self.objective_composition:
                    self.objective_composition[mol_id] += coeff1 * coeff2
                else:
                    self.objective_composition[mol_id] = coeff1 * coeff2

        ## Get initial internal state from initial_mass
        initial_metabolite_mass = initial_parameters.get('initial_mass', self.defaults['initial_mass'])
        mw = self.fba.molecular_weights
        composition = {
            mol_id: (-coeff if coeff < 0 else 0)
            for mol_id, coeff in self.objective_composition.items()}
        composition_mass = get_fg_from_counts(composition, mw)
        scaling_factor = (initial_metabolite_mass / composition_mass).magnitude
        internal_state = {mol_id: int(coeff * scaling_factor)
            for mol_id, coeff in composition.items()}
        initial_mass = get_fg_from_counts(internal_state, mw)
        print('metabolism initial mass: {}'.format(initial_mass))

        ## Get external state from minimal_external fba solution
        external_state = {state_id: 0.0 for state_id in self.fba.external_molecules}
        external_state.update(self.fba.minimal_external)  # optimal minimal media from fba

        # save initial state
        self.initial_state = {
            'external': external_state,
            'internal': internal_state,
            'reactions': {state_id: 0 for state_id in self.reaction_ids},
            'exchange': {state_id: 0 for state_id in self.fba.external_molecules},
            'flux_bounds': {state_id: self.default_upper_bound
                            for state_id in self.constrained_reaction_ids},
        }

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
            'internal': {mol_id: {
                    'mass': self.fba.molecular_weights[mol_id]}
                for mol_id in self.internal_state_ids},
            'reactions': {rxn_id: {
                    'updater': 'set',
                    'divide': 'set'}
                for rxn_id in self.reaction_ids}}

        # derivers
        deriver_setting = [{
            'type': 'mass',
            'source_port': 'internal',
            'derived_port': 'global',
            'keys': self.internal_state_ids},
            {
            'type': 'globals',
            'source_port': 'global',
            'derived_port': 'global',
            'keys': []}
        ]

        return {
            'state': self.initial_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'deriver_setting': deriver_setting,
            'time_step': 2}

    def next_update(self, timestep, states):

        ## get the state
        external_state = states['external']
        constrained_reaction_bounds = states['flux_bounds']  # (units.mmol / units.L / units.s)
        mmol_to_counts = states['global']['mmol_to_counts'] * units.L / units.mmol

        ## get flux constraints
        # exchange_constraints based on external availability
        exchange_constraints = {mol_id: 0.0
            for mol_id, conc in external_state.items() if conc <= self.exchange_threshold}

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

        # timestep dependence on fluxes
        exchange_fluxes.update((mol_id, flux * timestep) for mol_id, flux in exchange_fluxes.items())
        internal_fluxes.update((mol_id, flux * timestep) for mol_id, flux in internal_fluxes.items())

        # update internal counts from objective flux
        # calculate added mass from the objective molecules' molecular weights
        objective_count = (objective_exchange * mmol_to_counts).magnitude
        internal_state_update = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if coeff2 < 0:  # pull out molecule if it is USED to make biomass (negative coefficient)
                    added_count = int(-coeff1 * coeff2 * objective_count)
                    internal_state_update[mol_id] = added_count

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
        }



# configs
def get_e_coli_core_config():
    metabolism_file = os.path.join('models', 'e_coli_core.json')
    return {'model_path': metabolism_file}

def get_iAF1260b_config():
    metabolism_file = os.path.join('models', 'iAF1260b.json')
    return {'model_path': metabolism_file}

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

    reversible = ['R6', 'R7', 'Rres']

    default_reaction_bounds = 1000.0

    exchange_bounds = {
        'A': -0.02,
        'D': 0.01,
        'E': 0.01,
        'F': -0.005,
        'H': -0.005,
        'O2': -0.1}

    initial_state = {
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0}}

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

# sim functions
def run_sim_save_network(config=get_toy_configuration(), out_dir='out/network'):
    metabolism = Metabolism(config)

    # initialize the process
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
    reactions = timeseries['reactions']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        flux = abs(np.mean(reactions[rxn_id][1:]))
        reaction_fluxes[rxn_id] = np.log(1000 * flux + 1.1)

    # define node type
    node_types = {rxn_id: 'reaction' for rxn_id in reaction_ids}
    node_types.update({mol_id: 'external_mol' for mol_id in external_mol_ids})
    node_types.update({rxn_id: 'objective' for rxn_id in objective.keys()})
    info = {
        'node_types': node_types,
        'reaction_fluxes': reaction_fluxes}

    nodes, edges = make_network(stoichiometry, info)
    save_network(nodes, edges, out_dir)

def run_metabolism(metabolism, settings):
    sim_settings = default_sim_settings
    sim_settings.update(settings)
    return simulate_process_with_environment(metabolism, sim_settings)

# plots
def plot_exchanges(timeseries, sim_config, out_dir='out', filename='exchanges'):
    # plot focused on exchanges with the environment

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
    ax1.set_ylabel('concentrations (logs)')
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
    ax2.set_ylabel('delta counts (log)')
    ax2.set_yscale('log')

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
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')

# energy carriers in BiGG models
BiGG_energy_carriers = [
    'atp_c',
    'gtp_c',
    'nad_c',
    'nadp_c',
    'fad_c',
]

def energy_synthesis_plot(timeseries, settings, out_dir, figname='energy_use'):
    # plot the synthesis of energy carriers in BiGG model output

    energy_reactions = settings.get('reactions', {})
    saved_reactions = timeseries['reactions']
    time_vec = timeseries['time']

    # get each energy carrier's total flux
    carrier_use = {}
    for reaction_id, coeffs in energy_reactions.items():
        reaction_ts = saved_reactions[reaction_id]

        for mol_id, coeff in coeffs.items():

            # save if energy carrier is used
            if coeff < 0:
                added_flux = [-coeff*ts for ts in reaction_ts]
                if mol_id not in carrier_use:
                    carrier_use[mol_id] = added_flux
                else:
                    carrier_use[mol_id] = [
                        sum(x) for x in zip(carrier_use[mol_id], added_flux)]

    # make the figure
    n_cols = 1
    n_rows = 1
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 2))
    grid = plt.GridSpec(n_rows, n_cols)

    # first subplot
    ax = fig.add_subplot(grid[0, 0])
    for mol_id, series in carrier_use.items():
        ax.plot(time_vec, series, label=mol_id)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_title('energy use')
    ax.set_xlabel('time ($s$)')
    ax.set_ylabel('$(mmol*L^{{{}}}*s^{{{}}})$'.format(-1, -1))  # TODO -- use unit schema in figures

    # save figure
    fig_path = os.path.join(out_dir, figname)
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')


# tests
default_sim_settings = {
    'environment_port': 'external',
    'exchange_port': 'exchange',
    'environment_volume': 1e-6,  # L
    'timestep': 1,
    'timeline': [(10, {})]}

def test_toy_metabolism():
    regulation_logic = {
        'R4': 'if (external, O2) > 0.1 and not (external, F) < 0.1'}

    toy_config = get_toy_configuration()
    transport = toy_transport()  # TODO -- this is no longer running in the test.

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

    settings = default_sim_settings
    settings.update({'timeline': timeline})
    return simulate_process_with_environment(toy_metabolism, settings)

def test_BiGG_metabolism(config=get_iAF1260b_config(), settings={}):
    metabolism = Metabolism(config)
    run_metabolism(metabolism, settings)

def test_config(config=get_toy_configuration()):
    # configure metabolism process
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


reference_sim_settings = {
    'environment_port': 'external',
    'exchange_port': 'exchange',
    'environment_volume': 1e-5,  # L
    'timestep': 1,
    'timeline': [(20, {})]}

def test_metabolism_similar_to_reference():
    config = get_iAF1260b_config()
    metabolism = Metabolism(config)
    timeseries = run_metabolism(metabolism, reference_sim_settings)
    reference = load_timeseries(
        os.path.join(REFERENCE_DATA_DIR, NAME + '.csv'))
    assert_timeseries_close(timeseries, reference)



if __name__ == '__main__':
    out_dir = os.path.join(TEST_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # configure BiGG metabolism
    config = get_iAF1260b_config()
    metabolism = Metabolism(config)

    # simulation settings
    timeline = [(2520, {})] # 2520 sec (42 min) is the expected doubling time in minimal media
    sim_settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'environment_volume': 1e-5,  # L
        'timestep': 1,
        'timeline': timeline}

    # run simulation
    timeseries = run_metabolism(metabolism, sim_settings)
    save_timeseries(timeseries, out_dir)

    volume_ts = timeseries['global']['volume']
    mass_ts = timeseries['global']['mass']
    print('volume growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    print('mass growth: {}'.format(mass_ts[-1] / mass_ts[0]))

    # plot settings
    plot_settings = {
        'max_rows': 30,
        'remove_zeros': True,
        'skip_ports': ['exchange', 'flux_bounds', 'reactions'],
        'overlay': {
            'reactions': 'flux_bounds'}}

    # make plots from simulation output
    plot_simulation_output(timeseries, plot_settings, out_dir, 'BiGG_simulation')
    plot_exchanges(timeseries, sim_settings, out_dir)

    # make plot of energy reactions
    stoichiometry = metabolism.fba.stoichiometry
    energy_carriers = [get_synonym(mol_id) for mol_id in BiGG_energy_carriers]
    energy_reactions = get_reactions(stoichiometry, energy_carriers)
    energy_plot_settings = {'reactions': energy_reactions}
    energy_synthesis_plot(timeseries, energy_plot_settings, out_dir)

    # make a gephi network
    run_sim_save_network(get_iAF1260b_config(), out_dir)
