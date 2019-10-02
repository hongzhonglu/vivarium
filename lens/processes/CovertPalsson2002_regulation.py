from __future__ import absolute_import, division, print_function

import os
import csv

from lens.actor.process import Process, dict_merge
from lens.data.spreadsheets import load_tsv
from lens.data.helper import mols_from_reg_logic
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
        external_state = self.add_e_to_dict(states['external'])
        total_state = dict_merge(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.iteritems()}

        regulatory_state = {mol_id: regulatory_logic(boolean_state)
                            for mol_id, regulatory_logic in self.regulation_logic.iteritems()}

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

    def add_e_key(self, molecule_ids):
        return [mol_id + self.e_key for mol_id in molecule_ids]

    def remove_e_key(self, molecule_ids):
        return [mol_id.replace(self.e_key, '') for mol_id in molecule_ids]

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
        all_molecules = mols_from_reg_logic(data['covert2002_regulatory_proteins'])

        # remove external molecules from internal_molecules
        external_molecules = [mol_id for mol_id in all_molecules if self.e_key in mol_id]
        internal_molecules = [mol_id for mol_id in all_molecules if mol_id not in external_molecules]
        external_molecules = self.remove_e_key(external_molecules)

        return internal_molecules, external_molecules
