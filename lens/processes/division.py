from __future__ import absolute_import, division, print_function

import uuid

from lens.actor.process import Process


class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = []

        roles = {'internal': ['volume']}
        parameters = {'division_volume': 2.0}
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def next_update(self, timestep, states):
        volume = states['internal']['volume']
        if volume >= self.parameters['division_volume']:
            self.division = self.daughter_config(volume)

        return {'division': self.division}  # TODO -- division is not a role....

    def daughter_config(self, volume):
        config1 = {
            'id': str(uuid.uuid4()),
            'time': self.time(),
            'volume': volume * 0.5}
        config2 = {
            'id': str(uuid.uuid4()),
            'time': self.time(),
            'volume': volume * 0.5}
        return [config1, config2]
