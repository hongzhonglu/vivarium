from __future__ import absolute_import, division, print_function

import os
import csv
from scipy import constants
from itertools import ifilter

from lens.actor.process import Process
from lens.data.spreadsheets import JsonReader
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis

TSV_DIALECT = csv.excel_tab

DATA_DIR = os.path.join('lens', 'data', 'flat')
FILENAME_STR = 'covert2002_'
LIST_OF_FILENAMES = (
    "reactions.tsv",
    "maintenance_biomass_fluxes.tsv",
    "transport.tsv",
    "exchange_fluxes.tsv",
    # "GLC_G6P_initial.tsv",
    "GLC_G6P_flux_bounds.tsv",
    )

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

## Default States
# GLC-G6P media
GLC_G6P_EXTERNAL = {
    'ACET': 0.0,
    'CO+2': 100.0,  # "units.mmol / units.L"
    'ETOH': 0.0,
    'FORMATE': 0.0,
    'GLC': 1.2209,  # "units.mmol / units.L"
    'GLYCEROL': 0.0,
    'LAC': 0.0,
    'LCTS': 0.0,
    'OXYGEN-MOLECULE': 100.0,  # "units.mmol / units.L"
    'PI': 100.0,  # "units.mmol / units.L"
    'PYR': 0.0,
    'RIB': 0.0,
    'SUC': 0.0,
}

# GLC-LCT media
GLC_LCT_EXTERNAL = {
    'G6P': 0,  # [m mol / L]
    'GLC': 1.2209,  # [m mol / L]
    'RIB': 0,  # [m mol / L]
    'GLYCEROL': 0,  # [m mol / L]
    'SUC': 0,  # [m mol / L]
    'PYR': 0,  # [m mol / L]
    'LAC': 0,  # [m mol / L]
    'LCTS': 3.4034,  # [m mol / L]
    'FORMATE': 0,  # [m mol / L]
    'ETOH': 0,  # [m mol / L]
    'ACET': 0,  # [m mol / L]
    'PI': 100,  # [m mol / L]
    'CO+2': 100,  # [m mol / L]
    'OXYGEN-MOLECULE': 100,  # [m mol / L]
}

INTERNAL = {'mass': 0.032}

# helper functions
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
    # TODO -- merge with rate_law_utilities.get_molecules_from_reactions
    molecules = set()
    for reaction, stoich in stoichiometry.iteritems():
        molecules.update(stoich.keys())
    return list(molecules)

def load_tsv(dir_name, file_name):
    # TODO -- use load_tsv from utils/data/knowledge_base
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

def merge_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


class Metabolism(Process):
    def __init__(self, initial_parameters={}):

        self.exchange_key = initial_parameters['exchange_key']

        # parameters  # TODO -- pass paramters in?
        self.nAvogadro = constants.N_A
        self.cell_density = 1100.0 * (units.g / units.L)

        # # initialize mass
        # self.dry_mass = 403.0 * units.fg
        # self.cell_mass = 1339.0 * units.fg

        self.load_data()
        all_molecule_ids = get_molecules_from_reactions(self.stoichiometry)
        self.internal_molecule_ids = [mol_id
            for mol_id in all_molecule_ids if mol_id not in self.external_molecule_ids + ['mass']]

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
        # TODO -- reconcile these with environment
        glc_g6p = True
        glc_lct = False

        if glc_g6p:
            external = GLC_G6P_EXTERNAL
        elif glc_lct:
            external = GLC_LCT_EXTERNAL

        internal = INTERNAL
        environment_deltas = [key + self.exchange_key for key in external.keys()]

        # declare the states
        external_molecules = merge_dicts(external,{key: 0 for key in environment_deltas})
        internal_molecules = merge_dicts(internal, {'volume': 1})  # fL TODO -- get volume with deriver?

        return {
            'environment_deltas': environment_deltas,
            'external': external_molecules,
            'internal': internal_molecules}

    def next_update(self, timestep, states):

        volume = states['internal']['volume']
        external_state = states['external']

        # TODO -- set transport bounds
        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [external_state.get(molID, 0.0) for molID in exchange_molecules]
        self.fba.setExternalMoleculeLevels(external_concentrations)

        # TODO -- check units on exchange flux
        exchange_fluxes = self.fba.getExternalExchangeFluxes() # TODO * timestep

        # update state based on internal and external concentrations
        countsToMolar = 1 / (self.nAvogadro * volume * 1e-15)  # convert volume fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.external_molecule_ids, delta_exchange_counts))

        # TODO -- update internal state mass
        update = {
            'external': {mol_id + self.exchange_key: delta
                for mol_id, delta in environment_deltas.iteritems()},
        }

        return update

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            full_filename = FILENAME_STR + filename
            data[attrName] = load_tsv(DATA_DIR, full_filename)

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

        self.objective = {"mass": 1.0}

        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.external_molecule_ids}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']

if __name__ == '__main__':
    Metabolism()
