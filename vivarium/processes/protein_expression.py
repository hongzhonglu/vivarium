from __future__ import absolute_import, division, print_function

import random

from vivarium.actor.process import Process
from vivarium.utils.units import units


class ProteinExpression(Process):
    '''
    a minimal protein expression process
    '''
    def __init__(self, initial_parameters={}):

        expression_rates = {
            'protein1': 1e-3,
            'protein2': 1e-2,
        }
        parameters = {'expression_rates': expression_rates}
        self.internal_states = list(expression_rates.keys())
        roles = {'internal': self.internal_states}

        parameters.update(initial_parameters)

        super(ProteinExpression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0 for state_id in self.internal_states}
        default_state = {'internal': internal}

        # default emitter keys
        default_emitter_keys = {'internal': self.internal_states}

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
        for state_id in internal.keys():
            if random.random() < self.parameters['expression_rates'][state_id]:
                internal_update[state_id] = 1

        return {'internal': internal_update}