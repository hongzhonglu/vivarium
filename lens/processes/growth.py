from __future__ import absolute_import, division, print_function

from lens.actor.process import Process


class Growth(Process):
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['mass', 'volume'],
        }
        parameters = {'growth_rate': 0.005}
        parameters.update(initial_parameters)
        super(Growth, self).__init__(roles, parameters)

    def default_state(self):
        mass = 1.339e-12  # g  # 1339 (wet mass in fg)
        default_state = {'mass': mass}
        return {
            'environment_deltas': [],
            'environment_ids': [],
            'internal': default_state}

    def default_emitter_keys(self):
        keys = {'internal': ['mass']}
        return keys

    def next_update(self, timestep, states):
        mass = states['internal']['mass']
        growth_rate = mass * self.parameters['growth_rate'] * timestep
        update = {'mass': growth_rate}
        return {'internal': update}