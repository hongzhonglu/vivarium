from __future__ import absolute_import, division, print_function

import os
import csv

from lens.actor.process import Process
from lens.data.spreadsheets import load_tsv
from lens.utils.regulation_logic import RegulatoryLogic

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


class Regulation(Process):
    def __init__(self, initial_parameters={}):
        self.load_data()
        all_molecule_ids = get_molecules_from_reactions(self.stoichiometry)
        self.internal_molecule_ids = [mol_id
                                      for mol_id in all_molecule_ids if
                                      mol_id not in self.external_molecule_ids + ['Biomass']]

        roles = {'internal': self.internal_molecule_ids + ['volume']}
        parameters = {}
        parameters.update(initial_parameters)

        super(Regulation, self).__init__(roles, parameters)

    def default_state(self):
        # TODO -- get default activity of all proteins. -- state requires all molecules listed in regulatory_proteins.tsv
        pass

    def next_update(self, timestep, states):
        pass

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        rc = RegulatoryLogic()
        self.regulation_logic = {reaction['Protein']: rc.get_logic_function(reaction['Regulatory Logic'])
            for reaction in data['covert2002_regulatory_proteins']}
