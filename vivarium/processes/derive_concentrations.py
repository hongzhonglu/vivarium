from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.compartment.process import Deriver
from vivarium.utils.units import units



class DeriveConcs(Deriver):
    """
    Process for deriving concentrations from counts
    """
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

        super(DeriveConcs, self).__init__(ports, parameters)

    def default_settings(self):
        volume = 1.2 * units.fL
        mmol_to_counts = (self.avogadro * volume).to('L/mmol')

        # default state
        default_state = {
            'global': {
                'volume': volume.magnitude,
                'mmol_to_counts': mmol_to_counts.magnitude}}

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'concentrations': {
                state_id : {
                    'updater': 'set',
                    'divide': 'set'}
                for state_id in self.ports['concentrations']}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def ports_schema(self):
        return {
            'global': {
                'volume': {
                    '_default': 0.0},
                'mmol_to_counts': {
                    '_default': 0.0}},
            'counts': {
                concentration: {
                    '_default': 0,
                    '_updater': 'set'}
                for concentration in self.concentration_keys},
            'concentrations': {
                concentration: {
                    '_default': 0.0,
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

        # counts = {port: state for port, state in states.items() if port not in ['concentrations', 'global']}

        # # concentration update
        # concentrations = {}
        # for port, states in counts.items():
        #     for state_id, count in states.items():
        #         concentrations[state_id] = count / mmol_to_counts

        for molecule, concentration in concentrations.items():
            assert concentration >= 0, 'derived {} concentration < 0'.format(molecule)

        return {
            'concentrations': concentrations}
