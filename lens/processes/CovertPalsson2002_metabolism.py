from __future__ import absolute_import, division, print_function

import os
import csv
from scipy import constants
from itertools import ifilter

from lens.actor.process import Process
from lens.reconstruction.spreadsheets import JsonReader
from lens.utils import units
from lens.utils.modular_fba import FluxBalanceAnalysis

TSV_DIALECT = csv.excel_tab

DATA_DIR = os.path.join('lens', 'reconstruction', 'CovertPalsson2002')

LIST_OF_FILENAMES = (
    "reactions.tsv",
    "regulatory_proteins.tsv",
    "maintenance_biomass_fluxes.tsv",
    "transport.tsv",
    "exchange_fluxes.tsv",
    "GLC_G6P_initial.tsv",
    "GLC_G6P_flux_bounds.tsv",
    )


COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

# helper functions
def get_reverse(reactions):
    reverse_stoichiometry = {}
    for reaction in reactions:
        if reaction['Reversible']:
            reaction_id = reaction['Reaction']
            stoich = {mol_id: -1 * coeff for mol_id, coeff in reaction['Stoichiometry'].iteritems()}
            reverse_stoichiometry[reaction_id + '_reverse'] = stoich
    return reverse_stoichiometry

def get_molecules_from_reactions(stoichiometry):
    # TODO -- merge with rate_law_utilities.get_molecules_from_reactions
    molecules = set()
    for reaction, stoich in stoichiometry.iteritems():
        molecules.update(stoich.keys())
    return list(molecules)

def load_tsv(dir_name, file_name):
    # TODO -- use load_tsv from utils/reconstruction/knowledge_base
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


class Metabolism(Process):
    def __init__(self, initial_parameters={}):

        # parameters  # TODO -- pass paramters in?
        self.nAvogadro = constants.N_A
        self.cell_density = 1100.0 * (units.g / units.L)

        # initialize mass
        self.dry_mass = 403.0 * units.fg
        self.cell_mass = 1339.0 * units.fg

        self.load_data()
        all_molecule_ids = get_molecules_from_reactions(self.stoichiometry)
        self.internal_molecule_ids = [mol_id for mol_id in all_molecule_ids if mol_id not in self.external_molecule_ids + ['Biomass']]
        # TODO -- what about Biomass?

        # initialize FBA
        self.fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=self.external_molecule_ids,
            objective=self.objective,
            objectiveType="standard",
            solver="glpk-linear",
        )

        roles = {
            'external': self.external_molecule_ids,
            'internal': self.internal_molecule_ids + ['volume']}
        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_state(self):
        glc_g6p = True
        glc_lct = False

        if glc_g6p:
            # GLC-G6P media
            media = {
                # Biomass = 0.032,    # [g / l]
                'ACxt': 0.0,
                'CO2xt': 100.0,  # "units.mmol / units.L"
                'ETHxt': 0.0,
                'FORxt': 0.0,
                'GLCxt': 1.2209,  # "units.mmol / units.L"
                'GLxt': 0.0,
                'LACxt': 0.0,
                'LCTSxt': 0.0,
                'O2xt': 100.0,  # "units.mmol / units.L"
                'PIxt': 100.0,  # "units.mmol / units.L"
                'PYRxt': 0.0,
                'RIBxt': 0.0,
                'SUCCxt': 0.0,
            }
        elif glc_lct:
            # GLC-LCT media
            media = {
                # Biomass: 0.032,    # [g / l]
                'G6Pxt': 0,  # [m mol / L]
                'GLCxt': 1.2209,  # [m mol / L]
                'RIBxt': 0,  # [m mol / L]
                'GLxt': 0,  # [m mol / L]
                'SUCCxt': 0,  # [m mol / L]
                'PYRxt': 0,  # [m mol / L]
                'LACxt': 0,  # [m mol / L]
                'LCTSxt': 3.4034,  # [m mol / L]
                'FORxt': 0,  # [m mol / L]
                'ETHxt': 0,  # [m mol / L]
                'ACxt': 0,  # [m mol / L]
                'PIxt': 100,  # [m mol / L]
                'CO2xt': 100,  # [m mol / L]
                'O2xt': 100,  # [m mol / L]
            }

        external_molecules_changes = [key + '_change' for key in media.keys()]

        # declare the states
        environment_state = media
        environment_state.update({key: 0 for key in external_molecules_changes})
        environment_state['volume'] = 10
        cell_state = {
            'volume': 1,
            'Biomass': 0.032,  # [g / l]
        }

        return {
            'external_molecules': external_molecules_changes,
            'external': environment_state,
            'internal': cell_state}

    def next_update(self, timestep, states):

        volume = states['internal']['volume']
        external_state = states['external']

        # TODO -- set transport bounds
        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [external_state.get(molID, 0.0) for molID in exchange_molecules]
        self.fba.setExternalMoleculeLevels(external_concentrations)

        exchange_fluxes = self.fba.getExternalExchangeFluxes() # TODO * timestep

        # update state based on internal and external concentrations
        # TODO -- should this be millimolar? Get units in environment
        countsToMolar = 1 / (self.nAvogadro * volume * 1e-15)  # convert volume, fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.external_molecule_ids, delta_exchange_counts))

        # TODO -- update internal state -- Biomass/volume?
        update = {
            'external': {mol_id + '_change': delta for mol_id, delta in environment_deltas.iteritems()}, # TODO (Eran) -- pass in this '_change' substring
        }

        return update

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
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

if __name__ == '__main__':
    Metabolism()
