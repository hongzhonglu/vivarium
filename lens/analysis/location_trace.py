from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from lens.analysis.analysis import Analysis, get_lattice

class LatticeTrace(Analysis):
    def __init__(self):
        super(LatticeTrace, self).__init__(analysis_type='environment')

    def get_data(self, client, query):
        query.update({'type': 'lattice'})
        history_data = client.find(query)
        history_data.sort('time')
        lattice_history = get_lattice(history_data)

        return lattice_history

    def analyze(self, experiment_config, history_data, output_dir):

        agent_ids = [key for key in history_data.keys() if key is not 'time']
        time_vec = [t / 3600 for t in history_data['time']]  # convert to hours
        edge_x = experiment_config['edge_length_x']
        edge_y = experiment_config['edge_length_y']
        scaling = 8 / min(edge_x, edge_y)

        # plot trajectories
        fig = plt.figure(figsize=(scaling*edge_x, scaling*edge_y))
        for agent_id in agent_ids:
            # get locations and convert to 2D array
            locations = history_data[agent_id]['location']
            locations_array = np.array(locations)
            x_coord = locations_array[:, 0]
            y_coord = locations_array[:, 1]

            plt.plot(x_coord, y_coord, 'b-')  # trajectory
            plt.plot(x_coord[0], y_coord[0], color=(0.0,0.8,0.0), marker='*')  # starting point
            plt.plot(x_coord[-1], y_coord[-1], color='r', marker='*')  #  ending point

        plt.xlim((0, edge_x))
        plt.ylim((0, edge_y))
        start = mlines.Line2D([], [], color=(0.0,0.8,0.0), marker='*', linestyle='None',
                                  markersize=10, label='start')
        end = mlines.Line2D([], [], color='r', marker='*', linestyle='None',
                                  markersize=10, label='end')
        plt.legend(handles=[start, end])

        plt.savefig(output_dir + '/location_trace')
        plt.close(fig)
