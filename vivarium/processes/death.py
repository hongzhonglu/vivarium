'''Simulate cell death'''

from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.actor.process import (
    convert_to_timeseries,
    plot_simulation_output,
    Process,
)


class CheckerInterface(object):
    '''Interface that should be subclassed by all survivability checkers

    Each subclass should check for a condition that might kill the cell.
    '''

    def __init__(self):
        self.needed_state_keys = {}

    def check_can_survive(self, states):
        '''Check whether the current state is survivable by the cell

        Arguments:
            states: The states of each role in the cell, as a
                dictionary.

        Returns:
            True if the cell can survive, False if it cannot.
        '''
        raise NotImplementedError(
            'Checker should implement check_can_survive')


class AntibioticChecker(CheckerInterface):

    def __init__(
        self, antibiotic_threshold=0.9, antibiotic_key='antibiotic'
    ):
        super(AntibioticChecker, self).__init__()
        self.threshold = antibiotic_threshold
        self.key = antibiotic_key
        self.needed_state_keys.setdefault(
            'internal', set()).add('antibiotic')

    def check_can_survive(self, states):
        concentration = states['internal'][self.key]
        if concentration > self.threshold:
            return False
        return True


CHECKER_CLASSES = {
    'antibiotic': AntibioticChecker,
}


class DeathFreezeState(Process):

    def __init__(self, initial_parameters={}):
        self.checkers = [
            CHECKER_CLASSES[name](**config_dict)
            for name, config_dict in initial_parameters.get(
                'checkers', []).items()
        ]
        roles = {'internal': set(['death_freeze_state'])}
        for checker in self.checkers:
            needed_keys = checker.needed_state_keys
            for role in needed_keys:
                keys = roles.setdefault(role, set())
                keys |= needed_keys[role]
        for role, keys in roles.items():
            roles[role] = list(keys)
        super(DeathFreezeState, self).__init__(roles, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'internal': {
                    'death_freeze_state': False,
                },
            },
            'emitter_keys': {},
            'updaters': {
                'internal': {'death_freeze_state': 'set'},
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        for checker in self.checkers:
            if not checker.check_can_survive(states):
                return {
                    'internal': {
                        'death_freeze_state': True
                    }
                }
        return {}


def test_death_freeze_state(end_time=10, antibiotic_step=1, asserts=True):
    THRESHOLD = 5
    parameters = {
        'checkers': {
            'antibiotic': {
                'antibiotic_threshold': THRESHOLD,
            }
        }
    }
    process = DeathFreezeState(parameters)
    states = process.default_settings()['state']
    states['internal']['antibiotic'] = 0

    time = 0
    timestep = 1
    saved_states = {
        time: copy.deepcopy(states)
    }

    while time < end_time:
        update = process.next_update(timestep, states)
        # Update is this easy only because the only updater is 'set'
        for role, state in states.items():
            state.update(update.get(role, {}))
        time += timestep
        states['internal']['antibiotic'] += antibiotic_step
        saved_states[time] = copy.deepcopy(states)

    if asserts:
        expected_death = THRESHOLD / antibiotic_step
        expected_saved_states = {
            time: {
                'internal': {
                    'death_freeze_state': (
                        False if time <= expected_death else True),
                    'antibiotic': time * antibiotic_step,
                }
            }
            for time in range(end_time)
        }

    return saved_states

def plot_death_freeze_state_test():
    out_dir = os.path.join('out', 'tests', 'death_freeze_state')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    saved_data = test_death_freeze_state(asserts=False)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_settings = {}
    plot_simulation_output(timeseries, plot_settings, out_dir)


if __name__ == '__main__':
    plot_death_freeze_state_test()
