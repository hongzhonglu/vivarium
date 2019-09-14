from __future__ import absolute_import, division, print_function

import time
import os
import csv
import numpy as np
from scipy import constants

from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.actor.inner import Simulation

TSV_DIALECT = csv.excel_tab
EXTERNAL_MOLECULES_FILE = os.path.join('lens', 'environment', 'condition', 'environment_molecules.tsv')

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 128, 0]]


KCAT_MAX = 1.4e6

stoichiometry = {
    'R1': {'A':-1, 'ATP':-1, 'B':1},
    'R2a': {'B':-1, 'ATP':2, 'NADH':2, 'C':1},
    'R2b': {'C':-1, 'ATP':-2, 'NADH':-2, 'B':1},
    'R3': {'B':-1, 'F':1},
    'R4': {'C':-1, 'G':1},
    'R5': {'G':-1, 'C':0.8, 'NADH':2},
    'R6': {'C':-1, 'ATP':2, 'D':3},
    'R7': {'C':-1, 'NADH':-4, 'E':3},
    'R8a': {'G':-1, 'ATP':-1, 'NADH':-2, 'H':1},
    'R8b': {'G':1, 'ATP':1, 'NADH':2, 'H':-1},
    'Rres': {'NADH':-1, 'O2':-1, 'ATP':1},
}

biomass_stoichiometry = {
    'v_biomass': {'C':1, 'F':1, 'H':1, 'ATP':10}
}

# in mmol/gDCW/hr
transport_limits = {
    'A': 21.0 * (units.mmol / units.g / units.h),
    'F': 5.0 * (units.mmol / units.g / units.h),
    'D': -12.0 * (units.mmol / units.g / units.h),
    'E': -12.0 * (units.mmol / units.g / units.h),
    'H': 5.0 * (units.mmol / units.g / units.h),
    'O2': 15.0 * (units.mmol / units.g / units.h),
}

enzymes = {
    'R1':'E1',
    'R2a':'E2a',
    'R2b':'E2b',
    'R3':'E3',
    'R4':'E4',
    'R5':'E5',
    'R6':'E6',
    'R7':'E7',
    'R8a':'E8a',
    'R8b':'E8b',
    'Rres':'Eres',
}

initial_concentrations = {
    'E1':10.0,
    'E2a':10.0,
    'E2b':10.0,
    'E3':10.0,
    'E4':10.0,
    'E5':10.0,
    'E6':10.0,
    'E7':10.0,
    'E8a':10.0,
    'E8b':10.0,
    'Eres':10.0,
}


class Metabolism(Simulation):
    ''' Metabolism surrogate '''

    def __init__(self, state):
        self.initial_time = state.get('time', 0.0)
        self.local_time = 0.0
        self.timestep = 1.0
        self.nAvogadro = constants.N_A
        self.volume = 1.0  # in fL
        self.division_time = 100
        self.division = []

        self.cell_density = 1100.0 * (units.g / units.L)
        self.dry_mass = 403.0 * units.fg
        self.cell_mass = 1339.0 * units.fg

        self.stoichiometry = stoichiometry
        self.biomass_stoichiometry = biomass_stoichiometry
        self.enzymes = enzymes
        self.concentrations = initial_concentrations

        self.exchange_molecules = transport_limits.keys()

        self.environment_change = {mol_id: 0.0 for mol_id in self.exchange_molecules}


    def update_state(self):

        # convert transport limits from mmol/gDCW/hr to mol/L
        coefficient = self.dry_mass / self.cell_mass * self.cell_density * (self.timestep * units.s)  # coeff is in g*s/L
        transport_limits_molar = {mol_id: (limit * coefficient).asNumber(units.mol/units.L)
                                  for mol_id, limit in transport_limits.iteritems()}

        fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=transport_limits.keys(),
            objective=self.biomass_stoichiometry['v_biomass'],
            objectiveType='standard',
            solver='glpk-linear',
        )

        fba.setExternalMoleculeLevels([transport_limits_molar[molID] for molID in self.exchange_molecules])

        for reactionID in self.stoichiometry:
            if reactionID in self.enzymes:
                enzymeID = self.enzymes[reactionID]
                if enzymeID in self.concentrations:
                    fba.setReactionFluxBounds(reactionID, upperBounds=KCAT_MAX * self.concentrations[enzymeID])

        exchange_fluxes = fba.getExternalExchangeFluxes()

        # update state based on internal and external concentrations
        countsToMolar = 1 / (self.nAvogadro * self.volume * 1e-15) # convert volume to L

        # update environmental nutrient counts
        # delta_counts = ((1 / countsToMolar) * exchange_fluxes).asNumber().astype(int)
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)

        # Get the deltas for environmental molecules
        environment_deltas = dict(zip(self.exchange_molecules, delta_exchange_counts))

        # Accumulate in environment_change
        self.accumulate_deltas(environment_deltas)

    def accumulate_deltas(self, environment_deltas):
        for molecule_id, count in environment_deltas.iteritems():
            self.environment_change[molecule_id] += count

    def check_division(self):
        # update division state based on time since initialization
        if self.local_time >= self.initial_time + self.division_time:
            self.division = [{'time': self.local_time}, {'time': self.local_time}]

        return self.division

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
            self.update_state()

        # self.check_division()
        self.local_time = run_until

        time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            # 'motile_force': self.motile_force,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': DEFAULT_COLOR,
            }
