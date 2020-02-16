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



class Division(Process):
    def __init__(self, initial_parameters={}):
        self.division = 0

        initial_state = initial_parameters.get('initial_state', {})
        initial_volume = initial_state.get('global', {}).get('volume', 1.2)  # L
        division_volume = initial_volume * 2

        roles = {'global': ['volume', 'division']}

        parameters = {'division_volume': division_volume}  # TODO -- make division at 2X initial_volume?  Pass this in from initial_parameters
        parameters.update(initial_parameters)

        super(Division, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        globals = {
            'volume': 1.2,
            'division': False}
        default_state = {'global': globals}

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'global': {
                'division': {
                    'updater': 'set',
                    'divide': 'zero'}}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        volume = states['global']['volume']

        if volume >= self.parameters['division_volume']:
            self.division = 1

        return {'global': {'division': self.division}}
