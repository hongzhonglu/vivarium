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
        self.initial_state = initial_parameters['initial_state']
        self.stoichiometry = initial_parameters['stoichiometry']
        self.reversible = initial_parameters['reversible']
        self.objective = initial_parameters['objective']
        self.external_molecules = initial_parameters['external_molecules']
        self.reaction_ids = self.stoichiometry.keys()

        # initialize fba
        self.fba = CobraFBA(dict(
            stoichiometry=self.stoichiometry,
            reversible=self.reversible,
            external_molecules=self.external_molecules,
            objective=self.objective,
            initial_state=self.initial_state))

        # assign internal and external roles
        roles = {
            'external': self.external_molecules,
            'internal': list(self.stoichiometry.keys()) + ['volume', 'mass'],
            'flux_bounds': self.reaction_ids}

        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0 for state_id in self.stoichiometry.keys()}
        default_state = {
            'external':  self.initial_state.get('external'),
            'internal': deep_merge(dict(internal), self.initial_state.get('internal')),
            'flux_bounds': [],
        }

        # default emitter keys
        default_emitter_keys = {
            'internal': self.reaction_ids,
            'external': self.fba.external_molecules,
        }

        # default updaters
        set_internal_states = list(self.stoichiometry.keys()) + ['mass']
        default_updaters = {
            'internal': {state_id: 'set' for state_id in set_internal_states},
            'external': {mol_id: 'accumulate' for mol_id in self.external_molecules},
            'flux_bounds': {rxn_id: 'set' for rxn_id in self.reaction_ids},
        }

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

    def next_update(self, timestep, states):

        internal_state = states['internal']
        external_state = states['external'] # add_str_to_keys(states['external'], self.external_key)
        mass = internal_state['mass'] * units.fg
        volume = mass.to('g') / self.density

        # Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
        # coefficient = dry_mass / cell_mass * self.density * (timestep * units.s)
        coefficient = mass * self.density * (timestep * units.s)
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume  # convert volume fL to L

        # set external constraints.
        # v_ex </= [external] (mmol/L) / biomass (gDCW/L) / timestep (s)  --- flux in (mmol/g/s)
        # mass is in g rather than g/L. So use density (g/L)
        external_molecule_ids = self.fba.external_molecules
        external_flux_constraint = {
            molecule: external_state[molecule] / (self.density.magnitude * timestep)
            for molecule in external_molecule_ids}

        # print('external_flux_constraint: {}'.format(external_flux_constraint))
        self.fba.constrain_external_flux(external_flux_constraint)

        ## solve the problem!
        growth_rate = self.fba.optimize()
        exchange_fluxes = self.fba.read_external_fluxes() # (units.mmol / units.L)
        internal_fluxes = self.fba.read_internal_fluxes() # (units.mmol / units.L)

        # calculate the new mass
        # new_mass = {'mass': mass.magnitude}
        new_mass = {'mass': (mass.magnitude * np.exp(growth_rate * timestep))}

        # calculate delta counts for external molecules based on exchange flux and volume
        environment_deltas = {
            reaction: (flux * mmol_to_count).magnitude
            for reaction, flux in exchange_fluxes.items()}

        # return update
        return {
            'internal': deep_merge(dict(new_mass), internal_fluxes),
            'external': environment_deltas}



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
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'BIOMASS': 1}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2', 'BIOMASS']

    objective = {'v_biomass': 1.0}

    initial_state = {
        'internal': {
            'mass': 1.0, #1339,
            'volume': 1E-15},
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0,
            'BIOMASS': 30.0
        }}

    config = {
        'stoichiometry': stoichiometry,
        'reversible': stoichiometry.keys(),
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state}
    return config


def test_config(get_config=get_toy_configuration):
    # configure metabolism process
    config = get_config()
    metabolism = Metabolism(config)

    # TODO -- constrain exchanges before testing

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
    print(metabolism.fba.read_external_fluxes())


def simulate_metabolism(metabolism, total_time=3600):

    # get initial state and parameters
    settings = metabolism.default_settings()
    state = settings['state']
    density = metabolism.density
    nAvogadro = metabolism.nAvogadro

    saved_data = {
        'internal': {state_id: [value] for state_id, value in state['internal'].items()},
        'external': {state_id: [value] for state_id, value in state['external'].items()},
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
        growth_rate = metabolism.fba.objective_value()

        # apply external update
        mmol_to_count = (nAvogadro.to('1/mmol') * volume_t0).to('L/mmol').magnitude
        for mol_id, exchange in update['external'].items():
            exchange_rate = exchange / mmol_to_count  # TODO -- per second?
            delta_conc = exchange_rate / growth_rate * mass_t0 * (np.exp(growth_rate * timestep) - 1)

            state['external'][mol_id] += delta_conc
            if state['external'][mol_id] < 0.0:  # this shouldn't be needed
                state['external'][mol_id] = 0.0

        # save state
        saved_data['time'].append(time)
        saved_data['internal']['volume'].append(volume_t0.magnitude)  # TODO -- get new volume
        for state_id, value in state['internal'].items():
            saved_data['internal'][state_id].append(value)
        for state_id, value in state['external'].items():
            saved_data['external'][state_id].append(value)

    return saved_data

def plot_output(saved_state, out_dir='out'):
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import math

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [t/60./60. for t in saved_state['time']]  # convert to hours
    n_data = [len(saved_state[key].keys()) for key in data_keys]

    n_col = 3
    n_rows = math.ceil(sum(n_data) / n_col) + 1

    # make figure
    fig = plt.figure(figsize=(n_col * 6, n_rows * 2))
    plot_idx = 1
    for key in data_keys:
        for state_id, series in sorted(saved_state[key].items()):
            ax = fig.add_subplot(n_rows, n_col, plot_idx)
            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + state_id)
            ax.ticklabel_format(style='sci',axis='y')
            ax.set_xlabel('time (hr)')
            plot_idx += 1

    # make figure output directory and save figure
    fig_path = os.path.join(out_dir, 'metabolism_dFBA')
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
    saved_state = simulate_metabolism(metabolism, total_time)

    external = saved_state['external']
    internal = saved_state['internal']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        print('rxn_id: {}'.format(rxn_id))
        flux = np.mean(internal[rxn_id][1:])  #[flux_simstep]
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
    saved_data = simulate_metabolism(metabolism, 1200)
    plot_output(saved_data, out_dir)

    ## make flux network from toy model
    save_network(metabolism, 20, out_dir)
