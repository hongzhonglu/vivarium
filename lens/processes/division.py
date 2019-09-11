from __future__ import absolute_import, division, print_function

import uuid

from lens.actor.process import Process


class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division_volume = 2.0
        self.division = []

    def generate_inner_update(self):
        # TODO (Eran) -- pass division to compartment
        return {
            'division': self.division,
        }

    def next_update(self, timestep, states):
        volume = states['internal']['volume']
        if volume >= self.division_volume:
            self.division = self.daughter_config(volume)

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