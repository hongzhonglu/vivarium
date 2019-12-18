from __future__ import absolute_import, division, print_function

import re

def get_mols_from_stoich(stoichiometry):
    molecules = set()
    for reaction, stoich in stoichiometry.items():
        molecules.update(stoich.keys())
    return list(molecules)

def get_mols_from_reg_logic(data):
    # get all molecules listed in "Regulatory Logic"
    regulation_molecules = set()
    regex_split = 'action is complex|surplus |active |IF |not |or |and |\(|\)| '
    for reaction in data:
        protein = set([reaction['Protein']])
        reg_logic = reaction['Regulatory Logic']
        reg_logic_parsed = re.split(regex_split, reg_logic)  # split string to remove patterns in regex_split
        reg_logic_parsed = set(filter(None, reg_logic_parsed))  # remove empty strings

        # add protein and regulatory molecules to set
        regulation_molecules.update(protein)
        regulation_molecules.update(reg_logic_parsed)

    # remove empty strings
    regulation_molecules.discard('')

    return list(regulation_molecules)