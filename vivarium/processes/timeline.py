from __future__ import absolute_import, division, print_function

import copy

from vivarium.utils.dict_utils import deep_merge
from vivarium.core.process import Process


class Timeline(Process):

    def __init__(self, initial_parameters={}):
        self.timeline = copy.deepcopy(initial_parameters['timeline'])

        # get ports
        ports = {'global': ['time']}
        for event in self.timeline:
            for state in list(event[1].keys()):
                port = {state[0]: [state[1]]}
                ports = deep_merge(ports, port)

        parameters = {
            'timeline': self.timeline}

        super(Timeline, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            'global': {
                'time': {
                    '_default': 0,
                    '_updater': 'accumulate'}}}

    def next_update(self, timestep, states):

        time = states['global']['time']

        update = {'global': {'time': timestep}}
        for (t, change_dict) in self.timeline:
            if time >= t:
                for state, value in change_dict.items():
                    port = state[0]
                    variable = state[1]
                    if port not in update:
                        update[port] = {}
                    update[port][variable] = {
                        '_value': value,
                        '_updater': 'set'}
                self.timeline.pop(0)

        return update