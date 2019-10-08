from __future__ import absolute_import, division, print_function

import os
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


def location_trace(exp_data, history_data, experiment_id, output_dir='out'):

    agent_ids = [key for key in history_data.keys() if key is not 'time']
    time_vec = [t / 3600 for t in history_data['time']]  # convert to hours
    edge_length = exp_data['edge_length']

    plt.figure(figsize=(8, 8))
    for agent_id in agent_ids:
        # get locations and convert to 2D array
        locations = history_data[agent_id]['location']
        locations_array = np.array(locations)

        x_coord = locations_array[:, 0]
        y_coord = locations_array[:, 1]
        angle = locations_array[:, 2]

        for i in range(len(locations)):
            plt.plot(x_coord[i:i + 2], y_coord[i:i + 2], 'b-')

    plt.xlim((0, edge_length))
    plt.ylim((0, edge_length))


    # make figure output directory and save figure
    fig_path = os.path.join(output_dir, experiment_id)
    if not os.path.isdir(fig_path):
        os.makedirs(fig_path)
    try:
        plt.tight_layout()
    except:
        pass
    plt.savefig(fig_path + '/location_trace')
    plt.clf()
