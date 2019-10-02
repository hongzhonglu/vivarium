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
    # TODO -- merge with rate_law_utilities.get_molecules_from_stoich
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

        internal_molecules = {key: True for key in self.internal}
        external_molecules = {key: True for key in self.external}

        return {
            'external': external_molecules,
            'internal': internal_molecules}

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        total_state = merge_dicts(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.iteritems()}

        # TODO -- passing boolean_state sometimes returns TypeError: 'bool' object is not callable
        regulatory_state = {mol_id: regulatory_logic(boolean_state) for mol_id, regulatory_logic in self.regulation_logic.iteritems()}

        update = regulatory_state
        return update

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        self.regulation_logic = {reaction['Protein']: rl.build_rule(reaction['Regulatory Logic'])
            for reaction in data['covert2002_regulatory_proteins']}

        import ipdb;
        ipdb.set_trace()

        # get all molecules listed in "Regulatory Logic" TODO -- make this a re-usable function (for metabolism too)q
        internal_molecules = set()
        regex_split = 'action is complex|surplus |active |IF |not |or |and |\(|\)| ' # TODO -- remove 'action is complex', 'surplus' from file
        for reaction in data['covert2002_regulatory_proteins']:
            reg_logic = reaction['Regulatory Logic']
            reg_logic_parsed = re.split(regex_split, reg_logic)  # split string to remove patterns in regex_split
            reg_logic_parsed = list(filter(None, reg_logic_parsed))  # remove empty strings
            internal_molecules.update(reg_logic_parsed)

        # remove external molecules from internal_molecules
        external_molecules = set([mol_id
            for mol_id in internal_molecules if '[e]' in mol_id])
        internal_molecules.difference_update(external_molecules)

        external_molecules = [mol.replace('[e]','') for mol in external_molecules]

        return list(internal_molecules), list(external_molecules)
