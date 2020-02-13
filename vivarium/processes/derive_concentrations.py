from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.actor.process import Process
from vivarium.utils.units import units



class DeriveConcs(Process):
    """
    Process for deriving concentrations from counts
    """
    def __init__(self, initial_parameters={}):
        self.avogadro = AVOGADRO

        roles = initial_parameters.get('roles')
        roles.update({
            'global': ['volume', 'mmol_to_counts']})

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveConcs, self).__init__(roles, parameters)

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

        # default updaters
        default_updaters = {
            'concentrations': {state_id: 'set'
                for state_id in self.roles['concentrations']}}


        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):

        # states
        mmol_to_counts = states['global']['mmol_to_counts']
        counts = {role: state for role, state in states.items() if role not in ['concentrations', 'global']}

        # concentration update
        concentrations = {}
        for role, states in counts.items():
            for state_id, count in states.items():
                concentrations[state_id] = count / mmol_to_counts

        for mol_id, conc in concentrations.items():
            assert conc >= 0, 'derived {} concentration < 0'.format(mol_id)

        return {
            'concentrations': concentrations}