'''Simulate cell death'''

from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.composition import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_compartment,
)
from vivarium.actor.process import (
    Process,
    initialize_state,
    load_compartment,
)
from vivarium.actor.composition import COMPARTMENT_STATE, get_schema


TOY_ANTIBIOTIC_THRESHOLD = 5.0
TOY_INJECTION_RATE = 2.0


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
        roles = {
            'internal': set(),
            'compartment': ['processes'],
        }
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
                'internal': {},
            },
            'emitter_keys': {},
            'updaters': {
                'compartment': {'processes': 'set'},
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        for checker in self.checkers:
            if not checker.check_can_survive(states):
                return {'compartment': {'processes': []}}
        return {}


class ToyAntibioticInjector(Process):

    def __init__(self, initial_parameters={}):
        self.injection_rate = initial_parameters.get(
            'injection_rate', 1.0)
        roles = {'internal': ['antibiotic']}
        super(ToyAntibioticInjector, self).__init__(roles, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'internal': {
                    'antibiotic': 0.0
                }
            },
            'emitter_keys': {'antibiotic'},
        }
        return default_settings

    def next_update(self, timestep, states):
        delta = timestep * self.injection_rate
        return {'internal': {'antibiotic': delta}}


def compose_toy_death(config):
    death_parameters = {
        'checkers': {
            'antibiotic': {
                'antibiotic_threshold': TOY_ANTIBIOTIC_THRESHOLD,
            }
        }
    }
    death_process = DeathFreezeState(death_parameters)
    injector_parameters = {
        'injection_rate': TOY_INJECTION_RATE,
    }
    injector_process = ToyAntibioticInjector(injector_parameters)
    processes = [
        {
            'death': death_process,
            'injector': injector_process,
        },
    ]
    topology = {
        'death': {
            'internal': 'cell',
            'compartment': COMPARTMENT_STATE,
        },
        'injector': {
            'internal': 'cell',
        }
    }
    init_state = {
        'cell': {
            'antibiotic': 0.0
        }
    }
    schema = get_schema(processes, topology)
    states = initialize_state(processes, topology, schema, init_state)
    options = {
        'topology': topology,
        'schema': schema,
    }
    return {
        'processes': processes,
        'states': states,
        'options': options,
    }


def test_death_freeze_state(end_time=10, asserts=True):
    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_toy_death, boot_config)
    settings = {
        'timeline': [(end_time, {})]
    }
    saved_states = simulate_compartment(compartment, settings)
    if asserts:
        # Add 1 because dies when antibiotic strictly above threshold
        expected_death = 1 + TOY_ANTIBIOTIC_THRESHOLD // TOY_INJECTION_RATE
        expected_saved_states = {
            time: {
                'cell': {
                    'antibiotic': (
                        time * TOY_INJECTION_RATE
                        if time <= expected_death
                        # Add one because death will only be detected
                        # the iteration after antibiotic above
                        # threshold. This happens because death and
                        # injector run "concurrently" in the composite,
                        # so their updates are applied after both have
                        # finished.
                        else (expected_death + 1) * TOY_INJECTION_RATE
                    )
                }
            }
            for time in range(end_time + 1)
        }
        assert expected_saved_states == saved_states

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
