from __future__ import absolute_import, division, print_function

import os
import csv

from lens.actor.process import Process, dict_merge
from lens.data.spreadsheets import load_tsv
from lens.data.helper import get_mols_from_reg_logic
import lens.utils.regulation_logic as rl
from lens.environment.lattice_compartment import remove_str_in_list, add_str_to_keys

TSV_DIALECT = csv.excel_tab

DATA_DIR = os.path.join('lens', 'data', 'flat')
LIST_OF_FILENAMES = (
    "covert2002_regulatory_proteins.tsv",
    )

def get_reverse(reactions):
    reverse_stoichiometry = {}
    for reaction in reactions:
        if reaction['Reversible']:
            reaction_id = reaction['Reaction']
            stoich = {mol_id: -1 * coeff
                      for mol_id, coeff in reaction['Stoichiometry'].iteritems()}
            reverse_stoichiometry[reaction_id + '_reverse'] = stoich
    return reverse_stoichiometry

def get_molecules_from_reactions(stoichiometry):
    molecules = set()
    for reaction, stoich in stoichiometry.iteritems():
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

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''

        # TODO -- which states should be boolean?
        internal_molecules = {key: 0 for key in self.internal}
        external_molecules = {key: 0 for key in self.external}
        internal = dict_merge(internal_molecules, {'volume': 1})

        return {
            'external': external_molecules,
            'internal': internal}

    def default_emitter_keys(self):
        keys = {
            'internal': self.internal,
            'external': self.external
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''

        updater_types = {
            'internal': {state_id: 'set' for state_id in self.regulation_logic.keys()},  # set updater for boolean values
            'external': {state_id: 'accumulate' for state_id in self.external}}  # all external values use default 'delta' udpater

        return updater_types

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = add_str_to_keys(states['external'], self.external_key)
        total_state = dict_merge(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.iteritems()}

        regulatory_state = {mol_id: regulatory_logic(boolean_state)
                            for mol_id, regulatory_logic in self.regulation_logic.iteritems()}

        return {
            'internal': regulatory_state,
            'external': {}}

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

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

def test_covert2002_regulation(total_time=3600):
    # configure process
    regulation = Regulation({})

    # get initial state and parameters
    state = regulation.default_state()

    saved_state = {
        'internal': {state_id: [value] for state_id, value in state['internal'].iteritems()},
        'external': {state_id: [value] for state_id, value in state['external'].iteritems()},
        'time': [0]}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < total_time:
        time += timestep

        # get update
        update = regulation.next_update(timestep, state)



        # save state
        saved_state['time'].append(time)
        # saved_state['internal']['volume'].append(volume_t0.magnitude)  # TODO -- get new volume
        for state_id, value in state['internal'].iteritems():
            saved_state['internal'][state_id].append(value)
        for state_id, value in state['external'].iteritems():
            saved_state['external'][state_id].append(value)


    data = {'saved_state': saved_state}
    return data

def plot_regulation_output(data, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    saved_state = data['saved_state']


    import ipdb; ipdb.set_trace()


    # save figure
    fig_path = os.path.join(out_dir, 'covert2002_regulation')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_covert2002_regulation()
    out_dir = os.path.join('out', 'CovertPalsson2002_regulation')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_regulation_output(saved_state, out_dir)
