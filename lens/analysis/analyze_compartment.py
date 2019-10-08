from __future__ import absolute_import, division, print_function

import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def compartment_analysis(data_dict, experiment_id, simulation_id, output_dir='out'):

    # remove series with all zeros
    zero_state = []
    for key1 in data_dict.iterkeys():
        if key1 is not 'time':
            for key2, series in data_dict[key1].iteritems():
                if all(v == 0 for v in series):
                    zero_state.append((key1,key2))
    for (key1, key2) in zero_state:
        del data_dict[key1][key2]

    data_keys = [key for key in data_dict.keys() if key is not 'time']

    time_vec = [t / 3600 for t in data_dict['time']]  # convert to hours
    n_data = [len(data_dict[key].keys()) for key in data_keys]
    n_rows = sum(n_data) + 2

    # make figure
    fig = plt.figure(figsize=(8, n_rows * 1.5))
    plot_idx = 1

    # plot data
    for key in data_keys:
        for mol_id, series in data_dict[key].iteritems():
            ax = fig.add_subplot(n_rows, 1, plot_idx)
            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            ax.set_xlabel('time (hrs)')
            plot_idx += 1

    # additional data as text
    zeros = ['{}[{}]'.format(state, role) for (role, state) in zero_state]
    ax = fig.add_subplot(n_rows, 1, plot_idx)
    ax.text(0.02, 0.1, 'states with all zeros: {}'.format(zeros), wrap=True)
    ax.axis('off')


    # make figure output directory and save figure
    fig_path = os.path.join(output_dir, experiment_id, simulation_id)
    if not os.path.isdir(fig_path):
        os.makedirs(fig_path)
    try:
        plt.tight_layout()
    except:
        pass
    plt.savefig(fig_path + '/compartment')
    plt.clf()
