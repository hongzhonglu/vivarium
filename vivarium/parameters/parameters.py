from __future__ import absolute_import, division, print_function

import os
import copy
import itertools
import math

import numpy as np
import matplotlib.pyplot as plt

from vivarium.compartment.composition import (
    load_compartment,
    simulate_compartment,
    simulate_with_environment)

# composites
from vivarium.composites.master import compose_master



null_emitter = {'emitter': 'null'}



def get_nested(dict, keys):
    d = dict
    for key in keys[:-1]:
        if key in d:
            d = d[key]
    try:
        value = d[keys[-1]]
    except:
        value = None
        print('value not found for: {}'.format(keys))

    return value

def set_nested(dict, keys, value, create_missing=True):
    d = dict
    for key in keys[:-1]:
        if key in d:
            d = d[key]
        elif create_missing:
            d = d.setdefault(key, {})
        else:
            return dict
    if keys[-1] in d or create_missing:
        d[keys[-1]] = value
    return dict

def get_parameters_logspace(min, max, number):
    '''
    get list of n parameters logarithmically spaced between min and max
    '''
    range = np.logspace(np.log10(min), np.log10(max), number, endpoint=True)
    return list(range)

def run_sim_get_output(new_compartment, output_values, settings, total_time, simulate_environment=False):
    if simulate_environment:
        sim_out = simulate_with_environment(new_compartment, settings)
    else:
        sim_out = simulate_compartment(new_compartment, settings)

    last_state = sim_out[total_time]

    output = []
    for output_value in output_values:
        output.append(get_nested(last_state, output_value))

    return output

def parameter_scan(config):
    '''
    Pass in a config (dict) with:
        - composite (function) -- a function for the composite compartment
        - scan_parameters (dict) -- each parameter location (tuple) mapped to a list of values
        - output_values (list) -- a list of output values (tuple) with the (port, key)
        - conditions (list) -- a list of state values (dict) with {port: {variable: value}}
            for the default state the condition is and empty dict, [{}]
        - options (dict)

    Returns a list of all parameter combinations, and a dictionary with output values for those parameters
    '''

    composite = config['composite']
    scan_params = config['scan_parameters']
    metrics = config['metrics']
    options = config.get('options', {})
    conditions = config.get('conditions', [{}])
    n_conditions = len(conditions)

    ## Set up the parameter
    # how many parameter sets for scan?
    n_values = [len(v) for v in scan_params.values()]
    n_combinations = np.prod(np.array(n_values))
    print('parameter scan size: {}'.format(n_combinations))

    # get default parameters from compartment
    compartment = load_compartment(composite)
    default_params = compartment.current_parameters()

    # make all parameter sets for scan
    param_keys = list(scan_params.keys())
    param_values = list(scan_params.values())
    param_combinations = list(itertools.product(*param_values))  # a list of all parameter combinations
    param_sets = [dict(zip(param_keys, combo)) for combo in param_combinations]  # list of dicts with {param: value}

    ## Simulation settings for scan
    total_time = options.get('time', 10)
    timestep = options.get('timestep', 1)
    simulate_environment = options.get('simulate_with_environment', False)
    simulation_settings = options.get('simulation_settings', {})
    settings = {
        'timestep': timestep,
        'total_time': total_time,
        'return_raw_data': True}
    settings.update(simulation_settings)

    # run all parameters, and save results
    results = []
    for params_index, param_set in enumerate(param_sets):
        # set up the parameters
        parameters = copy.deepcopy(default_params)
        for param_key, param_value in param_set.items():
            parameters = set_nested(parameters, param_key, param_value)

        ## Run the parameter set for each condition's state
        for condition_index, condition_state in enumerate(conditions):
            print('running parameter set {}/{}, condition {}/{}'.format(
                params_index + 1,
                n_combinations,
                condition_index+1,
                n_conditions))

            # make the compartment
            new_compartment = load_compartment(composite, parameters)

            # update the states
            for store, values in condition_state.items():
                new_compartment.states[store].assign_values(values)

            # run the sim with these parameters and condition
            try:
                output = run_sim_get_output(
                    new_compartment,
                    metrics,
                    settings,
                    total_time,
                    simulate_environment)

                result = {
                    'parameter_index': params_index,
                    'condition_index': condition_index,
                    'output': output}
                results.append(result)

            except:
                print('failed simulation: parameter set {}, condition {}'.format(params_index, condition_index))

    # organize data by metric
    output_data = {
        'results': results,
        'metrics': metrics,
        'parameter_sets': {idx: param_set for idx, param_set in enumerate(param_sets)},
        'conditions': {idx: condition for idx, condition in enumerate(conditions)}}

    return organize_param_scan_results(output_data)

def organize_param_scan_results(data):
    results = data['results']
    metrics = data['metrics']
    parameter_sets = data['parameter_sets']
    conditions = data['conditions']

    param_indices = list(range(0, len(parameter_sets)))
    condition_indices = list(range(0, len(conditions)))

    # organize the results by metric
    metric_data = {
        metric: {
            condition_index: [None for param in param_indices]
            for condition_index in condition_indices}
        for metric in metrics}

    for result in results:
        param_index = result['parameter_index']
        condition_index = result['condition_index']
        output = result['output']

        for metric_index, datum in enumerate(output):
            metric = metrics[metric_index]
            metric_data[metric][condition_index][param_index] = datum

    return {
        'metric_data': metric_data,
        'parameter_sets': parameter_sets,
        'conditions': conditions}

def plot_scan_results(results, out_dir='out', filename='parameter_scan'):
    metric_data = results['metric_data']
    parameter_sets = results['parameter_sets']
    conditions = results['conditions']
    parameter_indices = [idx for idx, param in enumerate(parameter_sets)]

    ## make figure
    n_cols = 1
    lines_per_row = 8
    base_rows = len(metric_data)
    param_rows = math.ceil(len(parameter_indices)/lines_per_row)
    condition_rows = math.ceil(len(conditions)/lines_per_row)
    n_rows = base_rows + param_rows + condition_rows
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 2))
    grid = plt.GridSpec(n_rows, n_cols)
    font = {'size': 6}
    plt.rc('font', **font)

    row_idx = 0
    col_idx = 0
    for metric, data, in metric_data.items():
        ax = fig.add_subplot(grid[row_idx, col_idx])
        for condition, param_data in data.items():
            ax.scatter(parameter_indices, param_data, label=condition)

        ax.legend(title='condition', bbox_to_anchor=(1.2, 1.0))
        ax.set_ylabel(metric)
        ax.set_xticks(parameter_indices)
        ax.set_xlabel('parameter set #')

        row_idx += 1

    # prepare text
    param_text_row = 1 / lines_per_row / param_rows
    cond_text_row = 1 / lines_per_row / condition_rows
    parameter_text = [
        '{}: {}'.format(param_idx, param_set)
        for param_idx, param_set in parameter_sets.items()]
    condition_text = [
        '{}: {}'.format(condition_idx, condition)
        for condition_idx, condition in conditions.items()]

    ## plot text
    # parameters
    ax = fig.add_subplot(grid[base_rows:base_rows+param_rows, :])
    ax.text(0, 1.0, 'parameters')
    for text_idx, param in enumerate(parameter_text):
        ax.text(0, 0.9-text_idx*param_text_row, param)
    ax.axis('off')

    # conditions
    ax = fig.add_subplot(grid[base_rows+param_rows:, :])
    ax.text(0, 1.0, 'conditions')
    for text_idx, condition in enumerate(condition_text):
        ax.text(0, 0.9-text_idx*cond_text_row, condition)
    ax.axis('off')


    ## save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')

def scan_master():
    composite = compose_master

    # define scanned parameters, to replace defaults
    scan_params = {
        ('transport',
         'kinetic_parameters',
         'EX_glc__D_e',
         ('internal','PTSG'),
         'kcat_f'):
            get_parameters_logspace(1e-3, 1e0, 6)
    }

    output_values = [
        ('reactions', 'EX_glc__D_e'),
        ('reactions', 'GLCptspp'),
        ('global', 'growth_rate')]

    # set up simulation settings and scan options
    timeline = [(30, {})]
    sim_settings = {
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'environment_volume': 1e-6,  # L
        'timeline': timeline}

    scan_options = {
        'simulate_with_environment': True,
        'simulation_settings': sim_settings}

    scan_config = {
        'composite': composite,
        'scan_parameters': scan_params,
        'output_values': output_values,
        'options': scan_options}
    results = parameter_scan(scan_config)

    return results



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'master_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    results = scan_master()
    plot_scan_results(results, out_dir)
