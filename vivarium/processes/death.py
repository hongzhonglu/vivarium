'''Simulate cell death'''

from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.composition import (
    plot_simulation_output,
    simulate_compartment,
    load_compartment)
from vivarium.compartment.process import (
    Process,
    initialize_state,
    COMPARTMENT_STATE,
)


TOY_ANTIBIOTIC_THRESHOLD = 5.0
TOY_INJECTION_RATE = 2.0


class DetectorInterface(object):
    '''Interface that should be subclassed by all death detectors

    Each subclass should check for a condition that might kill the cell.
    '''

    def __init__(self):
        self.needed_state_keys = {}

    def check_can_survive(self, states):
        '''Check whether the current state is survivable by the cell

        Arguments:
            states: The states of each port in the cell, as a
                dictionary.

        Returns:
            True if the cell can survive, False if it cannot.
        '''
        raise NotImplementedError(
            'Detector should implement check_can_survive')


class AntibioticDetector(DetectorInterface):

    def __init__(
        self, antibiotic_threshold=0.9, antibiotic_key='antibiotic'
    ):
        super(AntibioticDetector, self).__init__()
        self.threshold = antibiotic_threshold
        self.key = antibiotic_key
        self.needed_state_keys.setdefault(
            'internal', set()).add('antibiotic')

    def check_can_survive(self, states):
        concentration = states['internal'][self.key]
        if concentration > self.threshold:
            return False
        return True


DETECTOR_CLASSES = {
    'antibiotic': AntibioticDetector,
}


class DeathFreezeState(Process):

    def __init__(self, initial_parameters=None):
        if initial_parameters is None:
            initial_parameters = {}
        self.detectors = [
            DETECTOR_CLASSES[name](**config_dict)
            for name, config_dict in initial_parameters.get(
                'detectors', {}).items()
        ]
        # List of names of processes that will remain after death
        self.enduring_processes = initial_parameters.get(
            'enduring_processes', [])
        ports = {
            'internal': set(),
            'compartment': ['processes'],
            'global': ['dead'],
        }
        for detector in self.detectors:
            needed_keys = detector.needed_state_keys
            for port in needed_keys:
                keys = ports.setdefault(port, set())
                keys |= needed_keys[port]
        for port, keys in ports.items():
            ports[port] = list(keys)
        super(DeathFreezeState, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'internal': {},
                'global': {
                    'dead': 0,
                },
            },
            'emitter_keys': {'global': ['dead']},
            'updaters': {
                'compartment': {'processes': 'set'},
                'global': {'dead': 'set'},
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        for detector in self.detectors:
            if not detector.check_can_survive(states):
                cur_processes = states['compartment']['processes'][0]
                return {
                    'compartment': {
                        'processes': [
                            {
                                process_name: cur_processes[process_name]
                                for process_name in self.enduring_processes
                            }
                        ]
                    },
                    'global': {
                        'dead': 1,
                    },
                }
        return {}


class ToyAntibioticInjector(Process):

    def __init__(self, initial_parameters=None):
        if initial_parameters is None:
            initial_parameters = {}
        self.injection_rate = initial_parameters.get(
            'injection_rate', 1.0)
        self.antibiotic_name = initial_parameters.get(
            'antibiotic_name', 'antibiotic')
        ports = {'internal': [self.antibiotic_name]}
        super(ToyAntibioticInjector, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'internal': {
                    self.antibiotic_name: 0.0
                }
            },
            'emitter_keys': {'internal': [self.antibiotic_name]},
        }
        return default_settings

    def next_update(self, timestep, states):
        delta = timestep * self.injection_rate
        return {'internal': {self.antibiotic_name: delta}}


def compose_toy_death(config):
    death_parameters = {
        'detectors': {
            'antibiotic': {
                'antibiotic_threshold': TOY_ANTIBIOTIC_THRESHOLD,
            }
        },
        'enduring_processes': ['enduring_injector'],
    }
    death_process = DeathFreezeState(death_parameters)
    injector_parameters = {
        'injection_rate': TOY_INJECTION_RATE,
    }
    injector_process = ToyAntibioticInjector(injector_parameters)
    enduring_parameters = {
        'injection_rate': TOY_INJECTION_RATE,
        'antibiotic_name': 'enduring_antibiotic'
    }
    enduring_process = ToyAntibioticInjector(enduring_parameters)
    processes = [
        {
            'death': death_process,
            'injector': injector_process,
            'enduring_injector': enduring_process,
        },
    ]
    topology = {
        'death': {
            'internal': 'cell',
            'compartment': COMPARTMENT_STATE,
            'global': 'global',
        },
        'injector': {
            'internal': 'cell',
        },
        'enduring_injector': {
            'internal': 'cell',
        },
    }
    init_state = {
        'cell': {
            'antibiotic': 0.0,
            'enduring_antibiotic': 0.0,
        },
        'global': {
            'dead': 0,
        },
    }
    states = initialize_state(processes, topology, init_state)
    options = {
        'topology': topology,
    }
    return {
        'processes': processes,
        'states': states,
        'options': options,
    }


def test_death_freeze_state(end_time=10, asserts=True):
    compartment = load_compartment(compose_toy_death)
    settings = {
        'timeline': [(end_time, {})],
    }
    saved_states = simulate_compartment(compartment, settings)
    if asserts:
        # Add 1 because dies when antibiotic strictly above threshold
        expected_death = 1 + TOY_ANTIBIOTIC_THRESHOLD // TOY_INJECTION_RATE
        expected_saved_states = {
            'cell': {
                'antibiotic': [],
                'enduring_antibiotic': [],
            },
            'global': {
                'dead': [],
            },
            'time': [],
        }
        for i in range(end_time):
            time = i + 1
            expected_saved_states['cell']['antibiotic'].append(
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
            expected_saved_states['cell']['enduring_antibiotic'].append(
                time * TOY_INJECTION_RATE)
            expected_saved_states['global']['dead'].append(
                0 if time <= expected_death else 1)
            expected_saved_states['time'].append(time)
        assert expected_saved_states == saved_states

    return saved_states


def plot_death_freeze_state_test():
    out_dir = os.path.join('out', 'tests', 'death_freeze_state')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    timeseries = test_death_freeze_state(asserts=False)
    plot_settings = {}
    plot_simulation_output(timeseries, plot_settings, out_dir)


if __name__ == '__main__':
    plot_death_freeze_state_test()
