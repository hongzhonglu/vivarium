from __future__ import absolute_import, division, print_function

import time
import os
import csv
from scipy import constants

from itertools import ifilter
from lens.data.spreadsheets import JsonReader
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.actor.inner import Simulation

TSV_DIALECT = csv.excel_tab

DATA_DIR = os.path.join('lens', 'data', 'CovertPalsson2002')

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

DEFAULT_COLOR = [color/255 for color in [255, 128, 0]]

class Metabolism(Simulation):
    '''
    Metabolism process based on "Covert & Palsson (2002) Transcription Regulation in
    Constraints-based Metabolic Models of Escherichia coli."
    TODO -- add regulation
    TODO -- add transport (Kremling, Bettenbrock, & Gilles 2007)
    '''

    def __init__(self, state):
        self.initial_time = state.get('time', 0.0)
        self.local_time = 0.0
        self.timestep = 1.0
        self.nAvogadro = constants.N_A
        self.volume = 1.0  # in fL
        self.division_time = 100
        self.division = []

        self.cell_density = 1100.0 * (units.g / units.L)

        # initialize mass
        self.dry_mass = 403.0 * units.fg
        self.cell_mass = 1339.0 * units.fg

        self.load_data()

        # initialize FBA
        self.fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=self.exchange_molecules,
            objective=self.objective,
            objectiveType="standard",
            solver="glpk-linear",
        )

        self.environment_change = {mol_id: 0.0 for mol_id in self.exchange_molecules}

    def update_state(self, timestep):
        # TODO -- set transport bounds
        # for reaction_id, stoich in self.transport_stoichiometry.iteritems():
        #     self.fba.setReactionFluxBounds(reaction_id, lowerBounds=0, upperBounds=1.0)

        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [self.external_concentrations.get(molID, 0.0) for molID in exchange_molecules]
        self.fba.setExternalMoleculeLevels(external_concentrations)

        exchange_fluxes = self.fba.getExternalExchangeFluxes() # TODO * timestep

        # update state based on internal and external concentrations
        # TODO -- should this be millimolar? Get units in environment
        countsToMolar = 1 / (self.nAvogadro * self.volume * 1e-15) # convert volume, fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.exchange_molecules, delta_exchange_counts))

        # Accumulate in environment_change
        self.accumulate_deltas(environment_deltas)

    def accumulate_deltas(self, environment_deltas):
        for molecule_id, count in environment_deltas.iteritems():
            self.environment_change[molecule_id] += count

    def time(self):
        return self.local_time

    def apply_outer_update(self, update):
        self.external_concentrations = update['concentrations']
        self.environment_change = {}
        for molecule in self.external_concentrations.iterkeys():
            self.environment_change[molecule] = 0

    def run_incremental(self, run_until):
        # update state once per message exchange
        while self.time() < run_until:
            self.local_time += self.timestep
            self.update_state(self.timestep)
        self.local_time = run_until
        time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': DEFAULT_COLOR,
            }

    def get_reverse(self, reactions):
        reverse_stoichiometry = {}
        for reaction in reactions:
            if reaction['Reversible']:
                reaction_id = reaction['Reaction']
                stoich = {mol_id: -1 * coeff for mol_id, coeff in reaction['Stoichiometry'].iteritems()}
                reverse_stoichiometry[reaction_id + '_reverse'] = stoich
        return reverse_stoichiometry

    # TODO -- use load_tsv from utils/data/knowledge_base
    def load_tsv(self, dir_name, file_name):
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

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = self.load_tsv(DATA_DIR, filename)

        self.stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['reactions']}
        transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['transport']}
        maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['maintenance_biomass_fluxes']}
        self.stoichiometry.update(transport_stoichiometry)
        self.stoichiometry.update(maintenance_stoichiometry)

        reverse_stoichiometry = self.get_reverse(data['reactions'])
        reverse_transport_stoichiometry = self.get_reverse(data['transport'])
        self.stoichiometry.update(reverse_stoichiometry)
        self.stoichiometry.update(reverse_transport_stoichiometry)

        # TODO -- remove growth exchange flux from file?
        self.exchange_molecules = [reaction['Stoichiometry'].keys()[0]
            for reaction in data['exchange_fluxes'] if reaction['Reaction'] != 'Growth']

        self.objective = {"Biomass": 1.0}

        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.exchange_molecules}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']
