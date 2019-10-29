from __future__ import absolute_import, division, print_function

from lens.actor.process import Process


class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = 0

        roles = {'internal': ['volume']}
        parameters = {'division_volume': 2.4}  # TODO -- make division at 2X initial_volume?
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def default_state(self):
        default_state = {'division': False}
        return {
            'internal': default_state}

    def default_updaters(self):
        updater_types = {
            'internal': {'division': 'set'}}
        return updater_types

    def next_update(self, timestep, states):
        volume = states['internal']['volume']

        if volume >= self.parameters['division_volume']:
            self.division = 1

        return {'internal': {'division': self.division}}
