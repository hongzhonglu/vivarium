from __future__ import absolute_import, division, print_function

import random

from vivarium.actor.process import Process


def divide_condition(compartment):
    division_role = compartment.division_role
    division = compartment.states[division_role].state_for(['division'])
    if division.get('division', 0) == 0:  # 0 means false
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
                    # TODO -- this should not divide everything. if 'set' updater, just set value
                    divided[index][state_key][key] = value // 2 + (value % 2 if index == left else 0)

    print('divided {}'.format(divided))
    return divided



class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = 0

        initial_volume = initial_parameters.get('initial_volume', 1.2)
        division_volume = initial_volume * 2

        roles = {'internal': ['volume', 'division']}

        parameters = {'division_volume': division_volume}  # TODO -- make division at 2X initial_volume?  Pass this in from initial_parameters
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {
            'volume': 1,
            'division': False}
        default_state = {'internal': internal}

        # default emitter keys
        default_emitter_keys = {}

        # default updaters
        default_updaters = {
            'internal': {'division': 'set'}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        volume = states['internal']['volume']

        if volume >= self.parameters['division_volume']:
            self.division = 1

        return {'internal': {'division': self.division}}
