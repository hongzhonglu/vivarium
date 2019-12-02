from __future__ import absolute_import, division, print_function

import os
from scipy import constants
import numpy as np

from lens.actor.process import Process, deep_merge
from lens.utils.units import units
from lens.utils.cobra_fba import CobraFBA


# helper functions
# TODO -- make reverse in cobra_fba
def get_reverse(stoichiometry, reversible_reactions, reverse_key):
    '''
    stoichiometry (dict) -- {mol_id (str): stoich (dict)}
    reversible_reactions (list) -- reactions that need a reverse stoichiometry
    '''
    reverse_stoichiometry = {}
    for rxn_id in reversible_reactions:
        forward_stoich = stoichiometry[rxn_id]
        reverse_stoichiometry[rxn_id + reverse_key] = {
            mol_id: -1 * coeff for mol_id, coeff in forward_stoich.items()}
    return reverse_stoichiometry


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
        self.initial_state = initial_parameters['initial_state']
        self.stoichiometry = initial_parameters['stoichiometry']
        self.objective = initial_parameters['objective']

        # make reverse reactions
        # TODO -- make reverse in cobra_fba
        reverse_key = '_reverse'
        self.reversible_reactions = initial_parameters.get('reversible_reactions', [])
        reverse_stoichiometry = get_reverse(self.stoichiometry, self.reversible_reactions, reverse_key)
        self.stoichiometry.update(reverse_stoichiometry)
        self.reaction_ids = self.stoichiometry.keys()

        # transport
        self.external_molecule_ids = initial_parameters['external_molecules']

        # initialize fba
        self.fba = CobraFBA(dict(
            stoichiometry=self.stoichiometry,
            external_molecules=self.external_molecule_ids,
            objective=self.objective,
            initial_state=self.initial_state))

        # assign internal and external roles
        roles = {
            'external': self.external_molecule_ids,
            'internal': list(self.stoichiometry.keys()) + ['volume', 'mass']}

        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0 for state_id in self.stoichiometry.keys()}
        default_state = {
            'external':  self.initial_state.get('external'),
            'internal': deep_merge(dict(internal), self.initial_state.get('internal'))}

        # default emitter keys
        default_emitter_keys = {
            'internal': self.fba.reaction_ids(),
            'external': self.fba.external_molecules}

        # default updaters
        set_internal_states = list(self.stoichiometry.keys()) + ['mass']
        default_updaters = {
            'internal': {state_id: 'set' for state_id in set_internal_states},
            'external': {mol_id: 'accumulate' for mol_id in self.external_molecule_ids}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

    def next_update(self, timestep, states):

        internal_state = states['internal']
        external_state = states['external'] # add_str_to_keys(states['external'], self.external_key)
        mass = internal_state['mass'] * units.fg
        volume = mass.to('g') / self.density

        # set external constraints.
        # v_external </= [external] / ([biomass] * timestep)
        external_molecule_ids = self.fba.external_molecules
        external_flux_constraint = {
            molecule: external_state[molecule] / (mass.magnitude * timestep)
            for molecule in external_molecule_ids}

        self.fba.constrain_external_flux(external_flux_constraint)

        ## solve the problem!
        growth_rate = self.fba.optimize()
        exchange_fluxes = self.fba.read_external_fluxes(timestep)
        internal_fluxes = self.fba.read_internal_fluxes()

        # calculate the new mass
        # new_mass = {'mass': mass.magnitude}
        new_mass = {'mass': (mass.magnitude * np.exp(growth_rate * timestep))}

        # calculate delta counts for external molecules based on exchange flux and volume
        mmolToCount = self.nAvogadro.to('1/mmol') * volume  # convert volume fL to L
        environment_deltas = {
            reaction: (flux * mmolToCount).magnitude
            for reaction, flux in exchange_fluxes.items()}

        # return update
        return {
            'internal': deep_merge(dict(new_mass), internal_fluxes),
            'external': environment_deltas}


# tests and analyses of process
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

    # transportLimits = {
    #     'A': 21.,
    #     'F': 5.0,
    #     'D': -12.0,
    #     'E': -12.0,
    #     'H': 5.0,
    #     'O2': 15.0,
    # }

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
            'BIOMASS': 30.0}}

    config = {
        'stoichiometry': stoichiometry,
        # 'reversible_reactions': stoichiometry.keys(),
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state}
    return config


def test_toy(total_time=100):

    # configure metabolism process
    config = get_toy_configuration()
    metabolism = Metabolism(config)

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


        print('t = {} ------------------------'.format(time))


        import ipdb; ipdb.set_trace()


        # apply external update
        mmolToCount = (nAvogadro.to('1/mmol') * volume_t0).to('L/mmol').magnitude
        for mol_id, exchange in update['external'].items():
            exchange_rate = exchange / mmolToCount  # TODO -- per second?


            # TODO -- is growth rate needed here?
            delta_conc = exchange_rate / growth_rate * mass_t0 * (np.exp(growth_rate * timestep) - 1)


            print('{} external: {}'.format(mol_id, update['external']))
            print('{} delta_conc: {}'.format(mol_id, delta_conc))


            state['external'][mol_id] += delta_conc #* 0.00000001  # TODO -- scaling?
            if state['external'][mol_id] < 1e-9:  # this shouldn't be needed
                state['external'][mol_id] = 1e-9


        # save state
        saved_data['time'].append(time)
        saved_data['internal']['volume'].append(volume_t0.magnitude)  # TODO -- get new volume
        for state_id, value in state['internal'].items():
            saved_data['internal'][state_id].append(value)
        for state_id, value in state['external'].items():
            saved_data['external'][state_id].append(value)

    return saved_data

def plot_metabolism_output(saved_state, out_dir='out'):
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

def save_metabolic_network(time_total=60, out_dir='out'):
    # TODO -- make this function into an analysis
    import math
    from lens.utils.make_network import make_network, save_network

    # TODO -- make weigh based on average over time range?
    # flux_simstep = 1 #10 # simulation step for flux edge weights

    # # initialize process
    config = get_toy_configuration()  # TODO -- make get_config() configurable
    metabolism = Metabolism(config)
    stoichiometry = metabolism.stoichiometry
    reaction_ids = stoichiometry.keys()
    external_mol_ids = config['external_molecules']
    objective = config['objective']

    # run test to get simulation output
    saved_state = test_toy(time_total)  # TODO -- make test(sec) configurable
    # saved_state = data['saved_state']
    external = saved_state['external']
    internal = saved_state['internal']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        print('rxn_id: {}'.format(rxn_id))
        flux = np.mean(internal[rxn_id][1:])  #[flux_simstep]
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
    out_dir = os.path.join('out', 'tests', 'metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## test toy model
    saved_data = test_toy(5)
    plot_metabolism_output(saved_data, out_dir)

    ## make network of toy model
    # save_metabolic_network(2, out_dir)
