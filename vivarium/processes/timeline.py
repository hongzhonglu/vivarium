from __future__ import absolute_import, division, print_function

from vivarium.utils.dict_utils import deep_merge
from vivarium.compartment.process import Process

class Timeline(Process):

    def __init__(self, initial_parameters={}):
        self.timeline = initial_parameters['timeline']

        # get ports
        ports = {'global': ['time']}
        for event in self.timeline:
            # states = list(event[1].keys())
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

        for (t, change_dict) in timeline:
            if time >= t:
                for port_id, change in change_dict.items():
                    port = compartment.states.get(port_id)
                    port.assign_values(change)
                timeline.pop(0)



        import ipdb; ipdb.set_trace()
        # TODO -- need to set the new value!

        return {
            'global': {'time': timestep}
        }
