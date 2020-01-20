from __future__ import absolute_import, division, print_function

import os
import random
import copy

from vivarium.actor.process import Process, convert_to_timeseries, plot_simulation_output



default_step_size = 1
default_expression_rates = {
    'protein1': 1e-3,
    'protein2': 1e-2}

class MinimalExpression(Process):
    '''
    a minimal protein expression process

    parameters:
        expression_rates (dict) with {'mol_id': probability_of_expression (1/sec)}
    '''
    def __init__(self, initial_parameters={}):

        expression_rates = initial_parameters.get(
            'expression_rates', default_expression_rates)
        self.internal_states = list(expression_rates.keys())

        roles = {'internal': self.internal_states}

        parameters = {
            'expression_rates': expression_rates,
            'step_size': initial_parameters.get('step_size', default_step_size)
        }
        parameters.update(initial_parameters)

        super(MinimalExpression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = {state_id: 0 for state_id in self.internal_states}
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

        internal_update = {}
        for state_id in internal.keys():
            rate = self.parameters['expression_rates'][state_id]
            for step in range(n_steps):
                if random.random() < rate:
                    internal_update[state_id] = 1

        return {'internal': internal_update}



# test functions
def test_expression(end_time=10):
    toy_expression_rates = {
        'protein1': 1e-3,
        'protein2': 1e-2,
        'protein3': 1e-1}

    expression_config = {
        'expression_rates': toy_expression_rates}

    # load process
    expression = MinimalExpression(expression_config)

    # get initial state and parameters
    settings = expression.default_settings()
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
        update = expression.next_update(timestep, state)

        # apply update
        for role_id, states_update in update.items():
            if role_id not in skip_roles:
                for state_id, change in states_update.items():
                    state[role_id][state_id] += change

        saved_state[time] = copy.deepcopy(state)

    return saved_state



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'minimal_expression')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_expression(1000)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)
