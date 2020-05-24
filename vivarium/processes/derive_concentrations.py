from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.core.process import Deriver
from vivarium.utils.units import units



class DeriveConcentrations(Deriver):
    """
    Process for deriving concentrations from counts
    """

    defaults = {
        'concentration_keys': []}

    def __init__(self, initial_parameters={}):
        self.avogadro = AVOGADRO
        self.concentration_keys = self.or_default(
            initial_parameters, 'concentration_keys')

        ports = {
            'global': [
                'volume', 'mmol_to_counts'],
            'counts': self.concentration_keys,
            'concentrations': self.concentration_keys}

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveConcentrations, self).__init__(ports, parameters)

    def ports_schema(self):
        volume = 1.2
        mmol_to_counts = (self.avogadro * volume * units.fL).to('L/mmol').magnitude

        return {
            'global': {
                'volume': {
                    '_default': volume,
                    '_units': units.fL},
                'mmol_to_counts': {
                    '_default': mmol_to_counts,
                    '_units': units.L/units.mmol}},
            'counts': {
                concentration: {
                    '_divider': 'split'}
                for concentration in self.concentration_keys},
            'concentrations': {
                concentration: {
                    '_divider': 'set',
                    '_updater': 'set'}
                for concentration in self.concentration_keys}}

    def next_update(self, timestep, states):

        # states
        mmol_to_counts = states['global']['mmol_to_counts']
        counts = states['counts']

        # concentration update
        concentrations = {}
        for molecule, count in counts.items():
            concentrations[molecule] = count / mmol_to_counts

        for molecule, concentration in concentrations.items():
            assert concentration >= 0, 'derived {} concentration < 0'.format(molecule)

        return {
            'concentrations': concentrations}
