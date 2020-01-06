from __future__ import absolute_import, division, print_function

import os
import csv

from vivarium.actor.process import Process, deep_merge
from vivarium.data.spreadsheets import load_tsv
from vivarium.data.helper import get_mols_from_reg_logic
import vivarium.utils.regulation_logic as rl
from vivarium.environment.lattice_compartment import remove_str_in_list, add_str_to_keys

TSV_DIALECT = csv.excel_tab

DATA_DIR = os.path.join('vivarium', 'data', 'flat')
LIST_OF_FILENAMES = (
    "covert2002_regulatory_proteins.tsv",
    )

def get_reverse(reactions):
    reverse_stoichiometry = {}
    for reaction in reactions:
        if reaction['Reversible']:
            reaction_id = reaction['Reaction']
            stoich = {mol_id: -1 * coeff
                      for mol_id, coeff in reaction['Stoichiometry'].items()}
            reverse_stoichiometry[reaction_id + '_reverse'] = stoich
    return reverse_stoichiometry

def get_molecules_from_reactions(stoichiometry):
    molecules = set()
    for reaction, stoich in stoichiometry.items():
        molecules.update(stoich.keys())
    return list(molecules)


class Regulation(Process):
    def __init__(self, initial_parameters={}):
        self.external_key = '[e]'
        self.internal, self.external = self.load_data()

        roles = {
            'internal': self.internal,
            'external': self.external}
        parameters = {}
        parameters.update(initial_parameters)

        super(Regulation, self).__init__(roles, parameters)


    def default_settings(self):

        # default state
        internal_molecules = {key: 0 for key in self.internal}
        external_molecules = {key: 0 for key in self.external}
        internal = deep_merge(internal_molecules, {'volume': 1})
        default_state = {
            'external': external_molecules,
            'internal': internal}

        # default emitter keys
        default_emitter_keys = {
            'internal': self.internal,
            'external': self.external}

        # default updaters
        default_updaters = {
            'internal': {state_id: 'set' for state_id in self.regulation_logic.keys()},  # set updater for boolean values
            'external': {state_id: 'accumulate' for state_id in self.external}}  # all external values use default 'delta' udpater

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = add_str_to_keys(states['external'], self.external_key)
        total_state = deep_merge(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.items()}

        regulatory_state = {mol_id: regulatory_logic(boolean_state)
                            for mol_id, regulatory_logic in self.regulation_logic.items()}

        return {
            'internal': regulatory_state,
            'external': {}}

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            path = os.path.join(DATA_DIR, filename)
            data[attrName] = load_tsv(path)

        # make regulatory logic functions
        self.regulation_logic = {}
        for protein in data['covert2002_regulatory_proteins']:
            protein_id = protein['Protein']
            rule = rl.build_rule(protein['Regulatory Logic'])
            if rule({}):
                self.regulation_logic[protein_id] = rule

        # get all molecules listed in "Regulatory Logic"
        all_molecules = get_mols_from_reg_logic(data['covert2002_regulatory_proteins'])

        # remove external molecules from internal_molecules
        external_molecules = [mol_id for mol_id in all_molecules if self.external_key in mol_id]
        internal_molecules = [mol_id for mol_id in all_molecules if mol_id not in external_molecules]
        external_molecules = remove_str_in_list(external_molecules, self.external_key)

        return internal_molecules, external_molecules



def test_covert2002_regulation():
    timeline = [
        (0, {'external': {
            'GLC': 1}
        }),
        (100, {'external': {
            'GLC': 0,
            'OXYGEN-MOLECULE': 1}
        }),
        (200, {'external': {
            'SUC': 0}
        }),
        (300, {'external': {
            'ACET': 0,
            'OXYGEN-MOLECULE': 0}
        }),
        (500, {}),
    ]

    # configure process
    regulation = Regulation({})

    # get initial state and parameters
    settings = regulation.default_settings()
    state = settings['state']

    saved_state = {'internal': {}, 'external': {}, 'time': []}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.items():
                    state[key].update(change)

        update = regulation.next_update(timestep, state)
        saved_state['time'].append(time)

        # update external state
        for state_id, value in state['external'].items():
            if state_id in saved_state['external'].keys():
                saved_state['external'][state_id].append(value)
            else:
                saved_state['external'][state_id] = [value]

        # update internal state from update
        for state_id, value in update['internal'].items():
            if state_id in saved_state['internal'].keys():
                saved_state['internal'][state_id].append(value)
            else:
                saved_state['internal'][state_id] = [value]

    return saved_state

def bool_to_int(series):
    int_series = []
    for v in series:
        if v is True:
            int_series.append(1)
        elif v is False:
            int_series.append(0)
        else:
            int_series.append(v)
    return int_series

def plot_regulation_output(saved_state, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [float(t) / 3600 for t in saved_state['time']]  # convert to hours

    # make figure, with grid for subplots
    n_data = [len(saved_state[key].keys()) for key in data_keys]
    n_rows = sum(n_data)
    fig = plt.figure(figsize=(8, n_rows * 2.5))
    grid = plt.GridSpec(n_rows + 1, 1, wspace=0.4, hspace=1.5)

    # plot data
    plot_idx = 0
    for key in data_keys:
        for mol_id, series in sorted(saved_state[key].items()):
            ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)
            numeric_series = bool_to_int(series)

            ax.plot(time_vec, numeric_series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            ax.set_xlabel('time (hrs)')

            if key is 'internal':
                ax.set_yticks([0.0, 1.0])
                ax.set_yticklabels(["False", "True"])

            plot_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'covert2002_regulation')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_covert2002_regulation()
    out_dir = os.path.join('out', 'tests', 'CovertPalsson2002_regulation')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_regulation_output(saved_state, out_dir)
