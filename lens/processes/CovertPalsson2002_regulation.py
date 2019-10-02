from __future__ import absolute_import, division, print_function

import os
import re
import csv

from lens.actor.process import Process
from lens.data.spreadsheets import load_tsv
import lens.utils.regulation_logic as rl

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

def merge_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

class Regulation(Process):
    def __init__(self, initial_parameters={}):
        self.e_key = '[e]'
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
        internal = merge_dicts(internal_molecules, {'volume': 1})

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
        external_state = self.add_e_to_dict(states['external'])
        total_state = merge_dicts(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.iteritems()}

        regulatory_state = {mol_id: regulatory_logic(boolean_state)
                            for mol_id, regulatory_logic in self.regulation_logic.iteritems()}

        import ipdb; ipdb.set_trace()

        return {'internal': regulatory_state, 'external': {}}

    def add_e_to_dict(self, molecules_dict):
        ''' convert external state to compatible format by adding e_key'''
        e_dict = {}
        for key, value in molecules_dict.iteritems():
            if self.e_key in key:
                e_dict[key] = value
            else:
                e_dict[key + self.e_key] = value
        return e_dict

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        self.regulation_logic = {reaction['Protein']: rl.build_rule(reaction['Regulatory Logic'])
            for reaction in data['covert2002_regulatory_proteins']}

        # get all molecules listed in "Regulatory Logic" TODO -- make this a re-usable function (for metabolism too)
        internal_molecules = set()
        regex_split = 'action is complex|surplus |active |IF |not |or |and |\(|\)| '
        for reaction in data['covert2002_regulatory_proteins']:
            protein = set([reaction['Protein']])
            reg_logic = reaction['Regulatory Logic']
            reg_logic_parsed = re.split(regex_split, reg_logic)  # split string to remove patterns in regex_split
            reg_logic_parsed = set(filter(None, reg_logic_parsed))  # remove empty strings

            # add protein and regulatory molecules to set
            internal_molecules.update(protein)
            internal_molecules.update(reg_logic_parsed)

        # remove external molecules from internal_molecules
        external_molecules = set(mol_id for mol_id in internal_molecules if self.e_key in mol_id)
        internal_molecules.difference_update(external_molecules)
        external_molecules = [mol.replace(self.e_key,'') for mol in external_molecules]


        return list(internal_molecules), external_molecules
