from __future__ import absolute_import, division, print_function

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def compartment_analysis(data_dict, experiment_id, simulation_id, output_dir='out'):

    data_keys = [key for key in data_dict.keys() if key is not 'time']

    time_vec = [t / 3600 for t in data_dict['time']]  # convert to hours
    n_data = [len(data_dict[key].keys()) for key in data_keys]
    n_rows = sum(n_data)

    fig = plt.figure(figsize=(8, n_rows * 1.5))
    plot_idx = 1

    for key in data_keys:
        for mol_id, series in data_dict[key].iteritems():
            ax = fig.add_subplot(n_rows, 1, plot_idx)
            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            ax.set_xlabel('time (hrs)')
            plot_idx += 1

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
