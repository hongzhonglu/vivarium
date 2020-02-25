from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.actor.process import Process
from vivarium.utils.dict_utils import tuplify_role_dicts
from vivarium.utils.regulation_logic import build_rule
from vivarium.actor.composition import process_in_compartment, simulate_with_environment, convert_to_timeseries, \
    plot_simulation_output



default_step_size = 1

class MinimalExpression(Process):
    '''
    a minimal protein expression process.
    TO BE USED ONLY AS TRAINING WHEELS

    parameters:
        expression_rates (dict) with {'mol_id': probability_of_expression (1/sec)}
    '''
    def __init__(self, initial_parameters={}):

        expression_rates = initial_parameters.get('expression_rates')
        self.internal_states = list(expression_rates.keys()) if expression_rates else []
        regulation_logic = initial_parameters.get('regulation', {})
        self.regulation = {
            gene_id: build_rule(logic) for gene_id, logic in regulation_logic.items()}
        regulators = initial_parameters.get('regulators', [])
        internal_regulators = [state_id for role_id, state_id in regulators if role_id == 'internal']
        external_regulators = [state_id for role_id, state_id in regulators if role_id == 'external']

        roles = {
            'internal': self.internal_states + internal_regulators,
            'external': external_regulators,
            'concentrations': []}

        parameters = {
            'expression_rates': expression_rates,
            'step_size': initial_parameters.get('step_size', default_step_size)}

        parameters.update(initial_parameters)


        super(MinimalExpression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # TODO -- load in initial state, or have compartment set to 0
        internal = {state_id: 0 for state_id in self.internal_states}
        default_state = {'internal': internal}

        # default emitter keys
        default_emitter_keys = {'internal': self.internal_states}

        deriver_setting = [{
            'type': 'counts_to_mmol',
            'source_role': 'internal',
            'derived_role': 'concentrations',
            'keys': self.internal_states}]

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'deriver_setting': deriver_setting}

        return default_settings

    def next_update(self, timestep, states):
        internal = states['internal']
        step_size = self.parameters['step_size']
        n_steps = int(timestep / step_size)

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_role_dicts(states)
        regulation_state = {}
        for gene_id, reg_logic in self.regulation.items():
            regulation_state[gene_id] = reg_logic(flattened_states)

        internal_update = {state_id: 0 for state_id in internal.keys()}
        for state_id in internal.keys():
            if state_id in regulation_state and not regulation_state[state_id]:
                break
            rate = self.parameters['expression_rates'][state_id]
            for step in range(n_steps):
                if random.random() < rate:
                    internal_update[state_id] += 1

        return {
            'internal': internal_update}



# test functions
def test_expression(end_time=10):
    toy_expression_rates = {
        'protein1': 1e-2,
        'protein2': 1e-1,
        'protein3': 1e0}

    expression_config = {
        'expression_rates': toy_expression_rates}

    # load process
    expression = MinimalExpression(expression_config)


    settings = {
        'total_time': 100,
        # 'exchange_role': 'exchange',
        'environment_role': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(expression)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'minimal_expression')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_expression(1000)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)

