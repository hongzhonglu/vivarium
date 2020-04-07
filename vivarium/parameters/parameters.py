from __future__ import absolute_import, division, print_function

import os
import copy
import itertools

import numpy as np
import matplotlib.pyplot as plt

from vivarium.compartment.composition import load_compartment, simulate_compartment

# composites
from vivarium.composites.master import compose_master


null_emitter = {'emitter': 'null'}

def get_nested(dict, keys):
    d = dict
    for key in keys[:-1]:
        if key in d:
            d = d[key]
    if keys[-1] in d:
        value = d[keys[-1]]
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

def parameter_scan(composite, scan_params, output_values):

    n_values = [len(v) for v in scan_params.values()]
    n_combinations = np.prod(np.array(n_values))

    # get initial parameters
    compartment = load_compartment(composite, null_emitter)
    default_params = compartment.current_parameters()

    # make list with combinations of all scanned parameters,
    # nested in the compartment's default parameters
    param_keys = list(scan_params.keys())
    param_values = list(scan_params.values())
    param_combinations = list(itertools.product(*param_values))
    param_sets = []
    for combo in param_combinations:
        new_params = copy.deepcopy(default_params)
        for param_key, value in zip(param_keys, combo):
            new_params = set_nested(new_params, param_key, value)
        param_sets.append(new_params)

    # simulation settings
    total_time = 10
    settings = {
        'timestep': 1,
        'total_time': total_time}

    # run all parameters, and save results
    results = []
    for params in param_sets:
        params.update(null_emitter)
        new_compartment = load_compartment(composite, params)
        sim_out = simulate_compartment(new_compartment, settings)
        last_state = sim_out[total_time]

        output = []
        for output_value in output_values:
            output.append(get_nested(last_state, output_value))
        results.append(output)


    # organize the results
    param_combo_ids = [dict(zip(param_keys, combo)) for combo in param_combinations]
    results_dict = {output_id: [result[out_idx] for result in results]
        for out_idx, output_id in enumerate(output_values)}

    return {
        'parameter combination': param_combo_ids,
        'output': results_dict}

def scan_master():
    composite = compose_master

    # define scanned parameters, to replace defaults
    scan_params = {
        ('metabolism', 'model_path'): ['models/e_coli_core.json'],
        ('transport', 'kinetic_parameters', 'EX_glc__D_e', ('internal','PTSG'), 'kcat_f'): [-1e5, -3e4, -6e3, -3e3, -3e2, -3e1]
    }

    output_values = [
        ('reactions', 'EX_glc__D_e'),
        ('reactions', 'GLCpts'),
        ('cell', 'growth_rate')
    ]

    results = parameter_scan(composite, scan_params, output_values)

    return results

def set_axes(ax, show_xaxis=False):
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)

def plot_scan_results(results, out_dir='out'):
    parameter_ids = results['parameter combination']
    outputs = results['output']
    param_indexes = list(range(0, len(parameter_ids)))

    # make figure
    n_cols = 1
    base_rows = len(outputs)
    n_rows = base_rows + 2

    fig = plt.figure(figsize=(n_cols * 6, n_rows * 2))
    grid = plt.GridSpec(n_rows, n_cols)
    font = {'size': 6}
    plt.rc('font', **font)

    row_idx = 0
    col_idx = 0
    for output_id, output in outputs.items():
        ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)
        ax.plot(param_indexes, output)
        ax.set_ylabel(output_id)

        row_idx += 1
        if row_idx == base_rows:
            # if last row of column
            set_axes(ax, True)
            ax.set_xticks(param_indexes)
            ax.set_xlabel('parameter set #')
        else:
            set_axes(ax)

    # parameter ids as text
    text_row = 0.08
    parameter_idxs = ['{}: {}'.format(idx, param_id)
        for idx, param_id in enumerate(parameter_ids)]

    ax = fig.add_subplot(grid[base_rows:, :])
    for text_idx, param_idx in enumerate(parameter_idxs):
        ax.text(0, 0.9-text_idx*text_row, param_idx)
    ax.axis('off')

    # save figure
    fig_path = os.path.join(out_dir, 'parameter_scan')
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')



if __name__ == '__main__':
    out_dir = os.path.join('out', 'parameters', 'scan_master')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    results = scan_master()
    plot_scan_results(results, out_dir)
