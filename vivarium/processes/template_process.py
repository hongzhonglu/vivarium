from __future__ import absolute_import, division, print_function

from vivarium.actor.process import Process
from vivarium.utils.units import units


class Template(Process):
    '''
    Need to add a boot method for this process to vivarium/environment/boot.py for it to run on its own
    '''
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['states'],
            'external': ['states'],
        }
        parameters = {}

        super(Template, self).__init__(roles, parameters)

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

        updaters defines the updater type for each state in roles.
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

        # default updaters
        default_updaters = {'state': 'accumulate'}

        default_settings = {
            'process_id': 'template',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
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