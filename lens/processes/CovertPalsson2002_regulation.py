from __future__ import absolute_import, division, print_function

import os
import csv
from itertools import ifilter

from lens.actor.process import Process
from lens.data.spreadsheets import JsonReader

TSV_DIALECT = csv.excel_tab


DATA_DIR = os.path.join('lens', 'data', 'CovertPalsson2002')
LIST_OF_FILENAMES = (
    "reactions.tsv",
    "regulatory_proteins.tsv",
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

def load_tsv(dir_name, file_name):
    # TODO (Eran) -- use load_tsv from utils/data/knowledge_base
    file_path = os.path.join(dir_name, file_name)
    with open(file_path, 'rU') as tsvfile:
        reader = JsonReader(
            ifilter(lambda x: x.lstrip()[0] != "#", tsvfile),  # Strip comments
            dialect=TSV_DIALECT)
        attr_list = []
        fieldnames = reader.fieldnames
        for row in reader:
            attr_list.append({field: row[field] for field in fieldnames})
    return attr_list

def get_molecules_from_reactions(stoichiometry):
    # TODO -- merge with rate_law_utilities.get_molecules_from_reactions
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

        roles = {
            'internal': self.internal_molecule_ids + ['volume']}
        parameters = {}
        parameters.update(initial_parameters)

        super(Regulation, self).__init__(roles, parameters)

    def next_update(self, timestep, states):
        pass


    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}

        # TODO -- get regulatory logic from regulatory_proteins.tsv, reactions.tsv




        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        self.stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['reactions']}
        transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['transport']}
        maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['maintenance_biomass_fluxes']}
        self.stoichiometry.update(transport_stoichiometry)
        self.stoichiometry.update(maintenance_stoichiometry)

        reverse_stoichiometry = get_reverse(data['reactions'])
        reverse_transport_stoichiometry = get_reverse(data['transport'])
        self.stoichiometry.update(reverse_stoichiometry)
        self.stoichiometry.update(reverse_transport_stoichiometry)

        # TODO -- remove growth exchange flux from file?
        self.external_molecule_ids = [reaction['Stoichiometry'].keys()[0]
            for reaction in data['exchange_fluxes'] if reaction['Reaction'] != 'Growth']

        self.objective = {"Biomass": 1.0}

        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.external_molecule_ids}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']
