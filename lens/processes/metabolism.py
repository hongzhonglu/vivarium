from __future__ import absolute_import, division, print_function

import os
from scipy import constants
import numpy as np

from lens.actor.process import Process, deep_merge
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.environment.lattice_compartment import add_str_in_list, remove_str_in_list, add_str_to_keys


# helper functions
def get_reverse(stoichiometry, reversible_reactions, reverse_key):
    '''
    stoichiometry (dict) -- {mol_id (str): stoich (dict)}
    reversible_reactions (list) -- reactions that need a reverse stoichiometry
    '''
    reverse_stoichiometry = {}
    for rxn_id in reversible_reactions:
        forward_stoich = stoichiometry[rxn_id]
        reverse_stoichiometry[rxn_id + reverse_key] = {
            mol_id: -1 * coeff for mol_id, coeff in forward_stoich.iteritems()}
    return reverse_stoichiometry

def get_molecules_from_stoich(stoichiometry):
    molecules = set()
    for reaction, stoich in stoichiometry.iteritems():
        molecules.update(stoich.keys())
    return list(molecules)


class Metabolism(Process):
    '''
    A general metabolism process, which sets the FBA problem based on input configuration data.
    initial_parameters (dict) configures the process with the following keys/values:
        - initial_state (dict) -- the default state, with a dict for internal and external:
            {'external': external_state, 'internal': internal_state}
        - stoichiometry (dict) -- {reaction_id: stoichiometry dict}
        - objective (dict) -- stoichiometry dict to be optimized
        - external_key (str) -- the key (typically [e]) designating external states in the FBA problem
        - external molecules (list) -- the external molecules, without the added external_key
        - target_key (TODO)
        - flux_targets (TODO)
    '''
    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1/units.mol
        self.density = 1100 * units.g/units.L
        self.initial_state = initial_parameters['initial_state']
        self.stoichiometry = initial_parameters['stoichiometry']
        self.objective = initial_parameters['objective']

        # make reverse reactions
        self.reverse_key = '_reverse'
        self.reversible_reactions = initial_parameters.get('reversible_reactions', [])
        reverse_stoichiometry = get_reverse(self.stoichiometry, self.reversible_reactions, self.reverse_key)
        self.stoichiometry.update(reverse_stoichiometry)
        self.reaction_ids = self.stoichiometry.keys()

        # transport
        # self.transport_limits = initial_parameters.get('transport_limits')
        self.external_key = initial_parameters.get('external_key', '')
        self.target_key = initial_parameters.get('target_key', '')
        self.external_molecule_ids = initial_parameters['external_molecules']
        external_molecule_ids_e = add_str_in_list(self.external_molecule_ids, self.external_key)
        self.flux_targets = add_str_in_list(initial_parameters.get('flux_targets', []), self.target_key)

        # regulation
        self.regulation = initial_parameters.get('regulation')

        self.fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=external_molecule_ids_e,
            objective=self.objective,
            objectiveType="standard",
            solver="glpk-linear")

        # assign internal and external state ids
        roles = {
            'external': self.external_molecule_ids,
            'internal': self.stoichiometry.keys() + ['volume', 'mass']}
        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0 for state_id in self.stoichiometry.iterkeys()}
        flux_targets = {state_id: 0 for state_id in self.flux_targets}
        internal.update(flux_targets)
        default_state = {
            'external':  self.initial_state.get('external'),
            'internal': deep_merge(dict(internal), self.initial_state.get('internal'))}

        # default emitter keys
        default_emitter_keys = {
            'internal': self.fba.getReactionIDs(),
            'external': self.fba.getExternalMoleculeIDs()}

        # default updaters
        set_internal_states = self.stoichiometry.keys() + ['mass']
        default_updaters = {
            'internal': {state_id: 'set' for state_id in set_internal_states},
            'external': {mol_id: 'accumulate' for mol_id in self.external_molecule_ids}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        ''' timestep is assumed to be 1 second '''

        internal_state = states['internal']
        external_state = add_str_to_keys(states['external'], self.external_key)
        mass = internal_state['mass'] * units.fg
        volume = mass.to('g') / self.density

        # # wcEcoli exchange state assumes mmol/g.h. coefficient in g.s/L. FBA takes mmol/L
        # # Coefficient converts between flux (mol/g DCW/hr) basis and concentration (M) basis.
        # # assumes dryMass / cellMass = 0.3
        # coefficient = 0.3 * self.density * (timestep * units.s)

        # TODO -- define constrained, unconstrained (np.inf), and 0 exchanges
        # set external concentrations.
        external_molecule_ids = self.fba.getExternalMoleculeIDs()
        # self.fba.setExternalMoleculeLevels([external_state[molID] for molID in external_molecule_ids])
        self.fba.setExternalMoleculeLevels([external_state[molID] / 3600 for molID in external_molecule_ids])  # TODO -- handle unit conversions (x/3600) separately

        ## Regulation
        # get the regulatory state of the reactions TODO -- are reversible reactions regulated?
        total_state = deep_merge(dict(internal_state), external_state)
        boolean_state = {mol_id: (value>1e-5) for mol_id, value in total_state.iteritems()}  # TODO -- generalize the threshold
        regulatory_state = {rxn_id: regulatory_logic(boolean_state)
                            for rxn_id, regulatory_logic in self.regulation.iteritems()}

        # set reaction flux bounds, based on default bounds and regulatory_state
        # TODO -- set default_flux_bounds in init
        # TODO -- if negative, lower bounds should set bound on reverse reactions
        default_flux_bounds = np.inf
        flux_bounds = np.array([[0.0, default_flux_bounds] for rxn in self.reaction_ids])
        for rxn_index, rxn_id in enumerate(self.reaction_ids):
            if regulatory_state.get(rxn_id) is False:
                flux_bounds[rxn_index] = [0.0, 0.0]

        # pass flux bounds to FBA
        self.fba.setReactionFluxBounds(
            self.reaction_ids,
            lowerBounds=flux_bounds[:,0],
            upperBounds=flux_bounds[:,1])

        # # set flux bounds on targets.
        # # TODO -- only update fluxes that are different from the previous timestep
        # for flux_id in self.flux_targets:
        #     reaction_id = flux_id.replace(self.target_key, '')
        #     flux_target = internal_state[flux_id]
        #     if flux_target >= 0:
        #         self.fba.setReactionFluxBounds(reaction_id, 0, flux_target)
        #     elif flux_target < 0:
        #         # if negative, set -(flux_target) on reverse reaction
        #         if reaction_id in self.reversible_reactions:
        #             rev_reaction_id = reaction_id + self.reverse_key
        #             self.fba.setReactionFluxBounds(rev_reaction_id, 0, -flux_target)

        # get exchanges and growth
        exchange_fluxes = self.fba.getExternalExchangeFluxes() * timestep #* mass.to('g').magnitude  # TODO -- scale by mass
        growth_rate = self.fba.getBiomassReactionFlux()
        new_mass = {'mass': (mass.magnitude * np.exp(growth_rate * timestep))}

        # get delta counts for external molecules
        mmolToCount = self.nAvogadro.to('1/mmol') * volume  # convert volume fL to L
        delta_exchange_counts = (mmolToCount * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(remove_str_in_list(external_molecule_ids, self.external_key), delta_exchange_counts))

        # get reaction flux
        rxn_ids = self.fba.getReactionIDs()
        rxn_fluxes = self.fba.getReactionFluxes()
        rxn_dict = dict(zip(rxn_ids, rxn_fluxes))

        update = {
            'internal': deep_merge(dict(new_mass), rxn_dict),
            'external': environment_deltas}
        return update


# tests and analyses of process
def kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def get_configuration():
    stoichiometry = {
        "R1": {"A": -1, "ATP": -1, "B": 1},
        "R2a": {"B": -1, "ATP": 2, "NADH": 2, "C": 1},
        "R2b": {"C": -1, "ATP": -2, "NADH": -2, "B": 1},
        "R3": {"B": -1, "F": 1},
        "R4": {"C": -1, "G": 1},
        "R5": {"G": -1, "C": 0.8, "NADH": 2},
        "R6": {"C": -1, "ATP": 2, "D": 3},
        "R7": {"C": -1, "NADH": -4, "E": 3},
        "R8a": {"G": -1, "ATP": -1, "NADH": -2, "H": 1},
        "R8b": {"G": 1, "ATP": 1, "NADH": 2, "H": -1},
        "Rres": {"NADH": -1, "O2": -1, "ATP": 1}}

    external_molecules = ["A", "F", "D", "E", "H", "O2"]

    objective = {"v_biomass": {"C": 1, "F": 1, "H": 1, "ATP": 10}}

    # transportLimits = {
    #     "A": 21.,
    #     "F": 5.0,
    #     "D": -12.0,
    #     "E": -12.0,
    #     "H": 5.0,
    #     "O2": 15.0,
    # }

    initial_state = {
        'internal': {
            'mass': 1339,
            'volume': 1E-15},
        'external': {
            "A": 21.0,
            "F": 5.0,
            "D": 12.0,
            "E": 12.0,
            "H": 5.0,
            "O2": 100.0}}

    config = {
        'stoichiometry': stoichiometry,
        # 'reversible_reactions': stoichiometry.keys(),
        'external_molecules': external_molecules,
        'objective': objective['v_biomass'],
        'initial_state': initial_state}
    return config

def test_metabolism(total_time=3600):
    import lens.utils.regulation_logic as rl

    target_key = '__target'

    # transport kinetics
    transport_rates = {
        "R1": kinetic_rate("A", 1.2e-2, 20),   # A import
        "R3": kinetic_rate("F", 2e-3, 5),  # F export
        # "R6": kinetic_rate("D", 1e-3, 12),  # D export
        # "R7": kinetic_rate("E", 1e-3, km),  # E export
        "R8a": kinetic_rate("H", 2e-3, 4.5), # H export
        # "R8b": kinetic_rate("H", 1e-5, 0.1),  # H import
        # "Rres": kinetic_rate("O2", 4e-1, 200), # O2 import
    }

    # regulation
    regulation = {
        "R1": rl.build_rule('IF not (F)'),
    }

    # configure process
    config = get_configuration()
    config['flux_targets'] = transport_rates.keys()
    config['target_key'] = target_key
    config['regulation'] = regulation

    metabolism = Metabolism(config)
    target_rxn_ids = metabolism.flux_targets

    # get initial state and parameters
    settings = metabolism.default_settings()
    state = settings['state']
    density = metabolism.density
    nAvogadro = metabolism.nAvogadro

    saved_state = {
        'internal': {state_id: [value] for state_id, value in state['internal'].iteritems()},
        'external': {state_id: [value] for state_id, value in state['external'].iteritems()},
        'time': [0]}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < total_time:
        mass_t0 = state['internal']['mass']
        volume_t0 = (mass_t0 * units.fg).to('g') / density
        time += timestep

        # get update
        update = metabolism.next_update(timestep, state)
        state['internal'] = update['internal']
        growth_rate = metabolism.fba.getObjectiveValue()

        # apply external update
        mmolToCount = (nAvogadro.to('1/mmol') * volume_t0).to('L/mmol').magnitude
        for mol_id, exchange in update['external'].iteritems():
            exchange_rate = exchange / mmolToCount  # TODO -- per second?
            delta_conc = exchange_rate / growth_rate * mass_t0 * (np.exp(growth_rate * timestep) - 1)
            state['external'][mol_id] += delta_conc * 0.001  # TODO -- scaling?
            if state['external'][mol_id] < 0:  # this shouldn't be needed
                state['external'][mol_id] = 0

        # get new flux targets
        target_fluxes = {}
        for rxn_id, rate_law in transport_rates.iteritems():
            target_flux = rate_law(state['external'])
            target_fluxes[rxn_id + target_key] = target_flux
        state['internal'].update(target_fluxes)

        # save state
        saved_state['time'].append(time)
        saved_state['internal']['volume'].append(volume_t0.magnitude)  # TODO -- get new volume
        for state_id, value in state['internal'].iteritems():
            saved_state['internal'][state_id].append(value)
        for state_id, value in state['external'].iteritems():
            saved_state['external'][state_id].append(value)

    # get select timeseries
    mass_ts = saved_state['internal']['mass'][2:]
    # O2_ts = saved_state['external']['O2'][2:]
    F_ts = saved_state['external']['F'][2:]
    H_ts = saved_state['external']['H'][2:]
    # A_ts = saved_state['external']['A'][2:]
    # D_ts = saved_state['external']['D']  # can vary
    # E_ts = saved_state['external']['E']  # can vary

    # check values are strictly increasing
    assert all(i < j for i, j in zip(mass_ts, mass_ts[1:]))

    # check values are strictly decreasing
    assert all(i > j for i, j in zip(F_ts, F_ts[1:]))
    assert all(i > j for i, j in zip(H_ts, H_ts[1:]))
    # assert all(i > j for i, j in zip(O2_ts, O2_ts[1:]))
    # assert all(i > j for i, j in zip(A_ts, A_ts[1:]))

    data = {
        'saved_state': saved_state,
        'target_rxn_ids': target_rxn_ids,
        'target_key': target_key}
    return data

def plot_metabolism_output(data, out_dir='out'):
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import math

    saved_state = data['saved_state']
    target_rxn_ids = data.get('target_rxn_ids')
    target_key = data.get('target_key')
    targeted_states = remove_str_in_list(target_rxn_ids, target_key)

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [t/60./60. for t in saved_state['time']]  # convert to hours
    n_data = [len(saved_state[key].keys()) for key in data_keys]

    n_col = 3
    n_rows = math.ceil(sum(n_data) / n_col) + 1

    # make figure
    fig = plt.figure(figsize=(n_col * 8, n_rows * 2.5))
    plot_idx = 1
    for key in data_keys:
        for state_id, series in sorted(saved_state[key].iteritems()):
            if state_id not in target_rxn_ids:
                ax = fig.add_subplot(n_rows, n_col, plot_idx)
                ax.plot(time_vec, series)
                ax.title.set_text(str(key) + ': ' + state_id)
                ax.ticklabel_format(style='sci',axis='y')
                ax.set_xlabel('time (hr)')
                plot_idx += 1

                # if state has target, plot series in red
                if state_id in targeted_states:
                    target_id = state_id + target_key
                    target_series = saved_state[key][target_id]
                    ax.plot(time_vec, target_series, 'r', label='target')
                    ax.legend()

    # make figure output directory and save figure
    fig_path = os.path.join(out_dir, 'metabolism_dFBA')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')

def save_metabolic_network(out_dir='out'):
    # TODO -- make this function into an analysis
    import math
    from lens.utils.make_network import make_network, save_network

    flux_simstep = 10 # simulation step for flux edge weights

    # # initialize process
    config = get_configuration()
    metabolism = Metabolism(config)
    stoichiometry = metabolism.stoichiometry
    reaction_ids = stoichiometry.keys()
    external_mol_ids = config['external_molecules']
    objective = config['objective']

    # run test to get simulation output
    data = test_metabolism(60)
    saved_state = data['saved_state']
    external = saved_state['external']
    internal = saved_state['internal']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        flux = internal[rxn_id][flux_simstep]
        reaction_fluxes[rxn_id] = math.log(1000 * flux + 1.1)

    # transport node type
    node_types = {rxn_id: 'reaction' for rxn_id in reaction_ids}
    node_types.update({mol_id: 'external_mol' for mol_id in external_mol_ids})
    node_types.update({mol_id: 'objective' for mol_id in objective.keys()})
    info = {
        'node_types': node_types,
        'reaction_fluxes': reaction_fluxes}

    nodes, edges = make_network(stoichiometry, info)
    save_network(nodes, edges, out_dir)

if __name__ == '__main__':
    saved_data = test_metabolism()
    out_dir = os.path.join('out', 'tests', 'metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_metabolism_output(saved_data, out_dir)
    save_metabolic_network(out_dir)