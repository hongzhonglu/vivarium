from __future__ import absolute_import, division, print_function

from lens.actor.process import Process
from lens.utils.units import units


class Template(Process):
    '''
    Need to add a boot method for this process to lens/environment/boot.py for it to run on its own
    '''
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['states'],
            'external': ['states'],
        }
        parameters = {}

        super(Template, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
        default_state = {
            'external': states (dict) -- external states ids with default initial values
            'internal': states (dict) -- internal states ids with default initial values
        '''
        internal_state = {}
        external_state = {}
        return {
            'internal': internal_state,
            'external': external_state,
        }

    def default_emitter_keys(self):
        '''
        returns dictionary with:
        keys = {
            'internal': states (list), # a list of states to emit from internal
            'external': states (list), # a list of states to emit from external
        }
        '''
        keys = {
            'internal': ['states'],
            'external': ['states'],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta,
        which is accumulated and passed to the environment at every exchange step'''
        keys = {'state': 'accumulate'}
        return keys

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