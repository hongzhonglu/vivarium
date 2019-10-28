from __future__ import absolute_import, division, print_function

from lens.actor.process import Process


class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = False

        roles = {'internal': ['volume']}
        parameters = {'division_volume': 2.0}
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def default_state(self):
        default_state = {'division': False}
        return {
            'internal': default_state}

    def next_update(self, timestep, states):
        volume = states['internal']['volume']
        if volume >= self.parameters['division_volume']:
            self.division = True

        return {'internal': {'division': self.division}}
