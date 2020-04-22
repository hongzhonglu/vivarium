from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.compartment.process import Process
from vivarium.utils.units import units



class DeriveConcs(Process):
    """
    Process for deriving concentrations from counts
    """
    def __init__(self, initial_parameters={}):
        self.avogadro = AVOGADRO

        source_ports = initial_parameters.get('source_ports')
        target_ports = initial_parameters.get('target_ports')

        assert len(target_ports) == 1, 'DeriveConcs too many target ports'
        assert list(target_ports.keys())[0] == 'concentrations', 'DeriveConcs requires target port named concentrations'

        ports = {'global': ['volume', 'mmol_to_counts']}
        ports.update(source_ports)
        ports.update(target_ports)

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

    def next_update(self, timestep, states):

        # states
        mmol_to_counts = states['global']['mmol_to_counts']
        counts = {port: state for port, state in states.items() if port not in ['concentrations', 'global']}

        # concentration update
        concentrations = {}
        for port, states in counts.items():
            for state_id, count in states.items():
                concentrations[state_id] = count / mmol_to_counts

        for mol_id, conc in concentrations.items():
            assert conc >= 0, 'derived {} concentration < 0'.format(mol_id)

        return {
            'concentrations': concentrations}
