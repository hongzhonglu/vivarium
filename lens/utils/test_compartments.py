from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math

from lens.actor.process import Compartment
from lens.environment.lattice_compartment import LatticeCompartment

# TODO -- use minimal composite for pytest
toy_composite = {
    'processes': {},
    'states': {},
    'options': {}}

# TODO -- make a test_process, simulate_process

def test_compartment(composite=toy_composite):
    boot_config = {}
    composite_config = composite(boot_config)
    processes = composite_config['processes']
    states = composite_config['states']
    options = composite_config['options']

    compartment = Compartment(processes, states, options)

    return compartment

def test_lattice_compartment(composite=toy_composite):
    boot_config = {}
    composite_config = composite(boot_config)
    processes = composite_config['processes']
    states = composite_config['states']
    options = composite_config['options']

    lattice_compartment = LatticeCompartment(processes, states, options)

    return lattice_compartment

def simulate_compartment(cmprt):
    # TODO -- need to update according to process updaters

    timestep = 1
    saved_state = {}
    for step in np.arange(10):
        cmprt.update(timestep)
        saved_state[step] = cmprt.current_state()

    return saved_state

def convert_to_timeseries(saved_states):
    # TODO -- timeseries assumes state is 1 deep, make more general state update

    time_vec = list(saved_states.keys())
    initial_state = saved_states[time_vec[0]]
    # roles = list(initial_state.keys())
    timeseries = {role: {state: []
        for state, initial in states.items()}
        for role, states in initial_state.items()}
    timeseries['time'] = time_vec

    for time, all_states in saved_states.items():
        for role, states in all_states.items():
            for state_id, state in states.items():
                timeseries[role][state_id].append(state)

    return timeseries

def set_axes(ax, show_xaxis=False):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)

def plot_simulation_output(saved_states, out_dir='out'):
    timeseries = convert_to_timeseries(saved_states)
    skip_keys = ['time']
    roles = [role for role in timeseries.keys() if role not in skip_keys]
    time_vec = timeseries['time']

    # make figure, with grid for subplots
    n_data = [len(timeseries[key]) for key in roles]
    n_cols = len(n_data)
    n_rows = max(n_data)

    fig = plt.figure(figsize=(n_cols * 8, n_rows * 2.5))
    grid = plt.GridSpec(n_rows + 1, n_cols, wspace=0.4, hspace=1.5)

    # plot data
    row_idx = 0
    col_idx = 0
    for role_idx, role in enumerate(roles):
        for state_id, series in sorted(timeseries[role].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            ax.plot(time_vec, series)
            ax.title.set_text(str(role) + ': ' + str(state_id))
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

            row_idx += 1
            if row_idx > n_data[role_idx]:
                row_idx = 0
                col_idx += 1
                set_axes(ax, True)
                ax.set_xlabel('time')
            else:
                set_axes(ax)

    # save figure
    fig_path = os.path.join(out_dir, 'simulation')
    plt.subplots_adjust(wspace=0.5, hspace=0.2)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')
