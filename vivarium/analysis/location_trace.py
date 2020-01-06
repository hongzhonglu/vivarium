from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from vivarium.analysis.analysis import Analysis, get_lattice

class LatticeTrace(Analysis):
    def __init__(self):
        super(LatticeTrace, self).__init__(analysis_type='environment')

    def get_data(self, client, query, options={}):
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

        markersize = 15

        # plot trajectories
        fig = plt.figure(figsize=(scaling*edge_x, scaling*edge_y))
        plt.rcParams.update({'font.size': 12, "font.family":"Times New Roman"})
        for agent_id in agent_ids:
            # get locations and convert to 2D array
            locations = history_data[agent_id]['location']
            locations_array = np.array(locations)
            x_coord = locations_array[:, 0]
            y_coord = locations_array[:, 1]

            plt.plot(x_coord, y_coord, label=agent_id)  # trajectory
            plt.plot(x_coord[0], y_coord[0], color=(0.0,0.8,0.0), marker='*', markersize=markersize)  # starting point
            plt.plot(x_coord[-1], y_coord[-1], color='r', marker='*', markersize=markersize)  #  ending point

        # set limits
        plt.xlim((0, edge_x))
        plt.ylim((0, edge_y))
        plt.xlabel(u'\u03bcm')
        plt.ylabel(u'\u03bcm')
        # specify the number of ticks
        plt.locator_params(axis='y', nbins=int(edge_y/10))
        plt.locator_params(axis='x', nbins=int(edge_x/10))

        # create legend for agent ids
        first_legend = plt.legend(title='agent ids', loc='center left', bbox_to_anchor=(1.01, 0.5), prop={'size': 12})
        ax = plt.gca().add_artist(first_legend)

        # create a legend for start/end markers
        start = mlines.Line2D([], [], color=(0.0,0.8,0.0), marker='*', linestyle='None',
                                  markersize=markersize, label='start')
        end = mlines.Line2D([], [], color='r', marker='*', linestyle='None',
                                  markersize=markersize, label='end')
        plt.legend(handles=[start, end], loc='upper left')


        plt.savefig(output_dir + '/location_trace', bbox_inches='tight')
        plt.close(fig)
