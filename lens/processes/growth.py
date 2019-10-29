from __future__ import absolute_import, division, print_function

import numpy as np

from lens.actor.process import Process


class Growth(Process):
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['mass', 'volume'],
        }
        parameters = {'growth_rate': 0.0006}
        parameters.update(initial_parameters)
        super(Growth, self).__init__(roles, parameters)

    def default_state(self):
        mass = 1339  # (wet mass in fg)
        default_state = {'mass': mass}
        return {
            'internal': default_state}

    def default_emitter_keys(self):
        keys = {'internal': ['mass']}
        return keys

    def default_updaters(self):
        updater_types = {'internal': {'mass': 'set'}}
        return updater_types

    def next_update(self, timestep, states):
        mass = states['internal']['mass']
        new_mass = mass * np.exp(self.parameters['growth_rate'] * timestep)
        return {'internal': {'mass': new_mass}}
