from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Process


class Template(Process):
    '''
    Need to add a boot method for this process to vivarium/environment/boot.py for it to run on its own
    '''

    defaults = {
        'parameters': {}
    }

    def __init__(self, initial_parameters={}):

        ports = {
            'internal': ['states'],
            'external': ['states'],
        }

        parameters = initial_parameters.get(
            'parameters', self.defaults['parameters'])

        super(Template, self).__init__(ports, parameters)

    def default_settings(self):
        '''
        state is a dictionary with:
        default_state = {
            'external': states (dict) -- external states ids with default initial values
            'internal': states (dict) -- internal states ids with default initial values

        emitter_keys is a dictionary with:
        keys = {
            'internal': states (list), # a list of states to emit from internal
            'external': states (list), # a list of states to emit from external
        }

        updaters defines the updater type for each state in ports.
        The default updater is to pass a delta,
        which is accumulated and passed to the environment at every exchange step
        '''

        # default state
        internal_state = {}
        external_state = {}
        default_state = {
            'internal': internal_state,
            'external': external_state}

        # default emitter keys
        default_emitter_keys = {
            'internal': ['states'],
            'external': ['states']}

        # schema -- define how each state is updater, divided, and its units
        schema = {
            'internal': {
                state_id : {
                    'updater': 'accumulate'}
                for state_id in self.ports['internal']}}

        # default derivers -- create a new derived port for these ports: keys
        deriver_setting = [{
            'type': 'mmol_to_counts',
            'source_port': 'internal',
            'derived_port': 'counts',
            'keys': self.ports['internal']}]

        default_settings = {
            'process_id': 'template',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'deriver_setting': deriver_setting,
            'time_step': 1.0}

        return default_settings


    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        internal_updates = 0
        external_updates = 0

        update = {
            'internal': {'states': internal_updates},
            'external': {'states': external_updates},
        }
        return update