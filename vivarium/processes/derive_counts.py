from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.actor.process import Process
from vivarium.utils.units import units



class DeriveCounts(Process):
    """
    Process for deriving counts from concentrations
    """
    def __init__(self, initial_parameters={}):
        self.avogadro = AVOGADRO

        roles = initial_parameters.get('roles')
        roles.update({
            'global': ['volume', 'mmol_to_counts']})

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveCounts, self).__init__(roles, parameters)

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
            'counts': {
                state_id : {
                    'updater': 'set'}
                for state_id in self.roles['counts']}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        mmol_to_counts = states['global']['mmol_to_counts']
        concentrations = {role: state for role, state in states.items() if role not in ['counts', 'global']}

        counts = {}
        for role, states in concentrations.items():
            for state_id, conc in states.items():
                counts[state_id] = int(conc * mmol_to_counts)

        return {
            'counts': counts}