from __future__ import absolute_import, division, print_function

import random

from vivarium.actor.process import Process
from vivarium.utils.units import units


class ProteinExpression(Process):
    '''
    a minimal protein expression process
    '''
    def __init__(self, initial_parameters={}):

        internal = ['protein']
        roles = {'internal': internal}
        parameters = {
            'expression_probability': 1e-3,
        }
        parameters.update(initial_parameters)

        super(ProteinExpression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal_state = {'protein': 0}
        default_state = {'internal': internal_state}

        # default emitter keys
        default_emitter_keys = {'internal': ['protein']}

        # default updaters
        default_updaters = {}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        internal = states['internal']

        internal_update = {}
        for state in internal.keys():
            if random.random() < self.parameters['expression_probability']:
                internal_update[state] = 1

        return {'internal': internal_update}