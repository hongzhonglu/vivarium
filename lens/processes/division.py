from __future__ import absolute_import, division, print_function

import random

from lens.actor.process import Process


def divide_condition(compartment):
    division = compartment.states['cell'].state_for(['division'])
    if division.get('division', 0) == 0:  # 0 is false
        divide = False
    else:
        divide = True
    return divide

def divide_state(compartment):
    divided = [{}, {}]
    for state_key, state in compartment.states.items():
        left = random.randint(0, 1)
        for index in range(2):
            divided[index][state_key] = {}
            for key, value in state.to_dict().items():
                if key == 'division':
                    divided[index][state_key][key] = 0
                else:
                    divided[index][state_key][key] = value // 2 + (value % 2 if index == left else 0)

    print('divided {}'.format(divided))
    return divided



class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = 0

        roles = {'internal': ['volume', 'division']}
        parameters = {'division_volume': 2.4}  # TODO -- make division at 2X initial_volume?  Pass this in from initial_parameters
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def default_state(self):
        default_state = {
            'volume': 1,
            'division': False}
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
