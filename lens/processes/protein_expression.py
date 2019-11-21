from __future__ import absolute_import, division, print_function

import random

from lens.actor.process import Process
from lens.utils.units import units


class ProteinExpression(Process):
    '''
    a minimal protein expression process
    '''
    def __init__(self, initial_parameters={}):

        internal = ['protein']
        roles = {'internal': internal}
        parameters = {
            'expression_probability': 1e-1,  # 1e-2,
        }
        parameters.update(initial_parameters)

        super(ProteinExpression, self).__init__(roles, parameters)

    def default_state(self):
        internal_state = {'protein': 0}
        return {'internal': internal_state}

    def default_emitter_keys(self):
        keys = {'internal': ['protein']}
        return keys

    def next_update(self, timestep, states):
        internal = states['internal']

        internal_update = {}
        for state in internal.keys():
            if random.random() < self.parameters['expression_probability']:
                internal_update[state] = 1

        return {'internal': internal_update}