from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.compartment.process import Process
from vivarium.utils.dict_utils import tuplify_port_dicts
from vivarium.utils.regulation_logic import build_rule
from vivarium.compartment.composition import (
    simulate_process_in_experiment,
    plot_simulation_output
)



class MinimalExpression(Process):
    '''
    a minimal protein expression process.
    TO BE USED ONLY AS TRAINING WHEELS

    parameters:
        expression_rates (dict) with {'mol_id': probability_of_expression (1/sec)}
    '''

    defaults = {
        'step_size': 1,
        'regulation': {},
        'concentrations_deriver_key': 'expression_concentrations_deriver',
    }

    def __init__(self, initial_parameters={}):

        expression_rates = initial_parameters.get('expression_rates')
        self.internal_states = list(expression_rates.keys()) if expression_rates else []
        regulation_logic = initial_parameters.get('regulation', self.defaults['regulation'])
        self.regulation = {
            gene_id: build_rule(logic) for gene_id, logic in regulation_logic.items()}
        regulators = initial_parameters.get('regulators', [])
        internal_regulators = [state_id for port_id, state_id in regulators if port_id == 'internal']
        external_regulators = [state_id for port_id, state_id in regulators if port_id == 'external']

        self.concentration_keys = self.internal_states + internal_regulators

        ports = {
            'internal': self.internal_states + internal_regulators,
            'external': external_regulators,
            'concentrations': [],
            'global': []}

        parameters = {
            'expression_rates': expression_rates,
            'step_size': initial_parameters.get('step_size', self.defaults['step_size'])}

        parameters.update(initial_parameters)

        self.concentrations_deriver_key = self.or_default(initial_parameters, 'concentrations_deriver_key')

        super(MinimalExpression, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            'internal': {
                state : {
                    '_updater': 'accumulate',
                    '_default': 0.0,
                    '_emit': True}
                for state in self.ports['internal']}}

    def derivers(self):
        return {
            self.concentrations_deriver_key: {
                'deriver': 'counts_to_mmol',
                'port_mapping': {
                    'global': 'global',
                    'counts': 'internal',
                    'concentrations': 'concentrations'},
                'config': {
                    'concentration_keys': self.concentration_keys}}}

    def next_update(self, timestep, states):
        internal = states['internal']
        step_size = self.parameters['step_size']
        n_steps = int(timestep / step_size)

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
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
    settings = {'total_time': end_time}
    return simulate_process_in_experiment(expression, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'minimal_expression')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    timeseries = test_expression(1000)
    plot_simulation_output(timeseries, {}, out_dir)

