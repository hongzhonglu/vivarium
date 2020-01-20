from __future__ import absolute_import, division, print_function

import os
import random
import copy

from vivarium.actor.process import Process, convert_to_timeseries, plot_simulation_output



default_step_size = 1

class MinimalDegradation(Process):
    '''
    a minimal protein degradation process

    parameters:
        degradation_rates (dict) with {'mol_id': probability_of_degradation (1/sec)}
    '''
    def __init__(self, initial_parameters={}):

        degradation_rates = initial_parameters.get('degradation_rates')
        self.internal_states = list(degradation_rates.keys()) if degradation_rates else []

        roles = {'internal': self.internal_states}

        parameters = {
            'degradation_rates': degradation_rates,
            'step_size': initial_parameters.get('step_size', default_step_size)
        }
        parameters.update(initial_parameters)

        super(MinimalDegradation, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # TODO -- load in initial state, or have compartment set to 0
        internal = {state_id: 100 for state_id in self.internal_states}
        default_state = {'internal': internal}

        # default emitter keys
        default_emitter_keys = {'internal': self.internal_states}

        # default updaters
        default_updaters = {}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        internal = states['internal']
        step_size = self.parameters['step_size']
        n_steps = int(timestep / step_size)

        internal_update = {state_id: 0 for state_id in internal.keys()}
        for state_id in internal.keys():
            rate = self.parameters['degradation_rates'][state_id]
            for step in range(n_steps):
                if random.random() < rate:
                    internal_update[state_id] -= 1

        return {'internal': internal_update}



# test functions
def test_degradation(end_time=10):
    toy_degradation_rates = {
        'protein1': 1e-3,
        'protein2': 1e-2,
        'protein3': 1e-1}

    degradation_config = {
        'degradation_rates': toy_degradation_rates}

    # load process
    degradation = MinimalDegradation(degradation_config)

    # get initial state and parameters
    settings = degradation.default_settings()
    state = settings['state']
    skip_roles = ['exchange']

    # initialize saved data
    saved_state = {}

    # run the simulation
    time = 0
    timestep = 5

    saved_state[0] = state
    while time < end_time:
        time += timestep
        # get update
        update = degradation.next_update(timestep, state)

        # apply update
        for role_id, states_update in update.items():
            if role_id not in skip_roles:
                for state_id, change in states_update.items():
                    state[role_id][state_id] += change

        saved_state[time] = copy.deepcopy(state)

    return saved_state



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'minimal_degradation')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_degradation(1000)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)
