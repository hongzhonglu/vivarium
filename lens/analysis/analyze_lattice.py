from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis, get_lattice

class LatticeTrace(Analysis):
    def __init__(self):
        super(LatticeTrace, self).__init__(analysis_type='lattice')

    def get_data(self, client):
        query = {'type': 'lattice'}
        history_data = client.find(query)
        history_data.sort('time')
        lattice_history = get_lattice(history_data)

        return lattice_history

    def analyze(self, experiment_config, history_data, output_dir):

        agent_ids = [key for key in history_data.keys() if key is not 'time']
        time_vec = [t / 3600 for t in history_data['time']]  # convert to hours
        edge_length = experiment_config['edge_length']

        plt.figure(figsize=(8, 8))
        for agent_id in agent_ids:
            # get locations and convert to 2D array
            locations = history_data[agent_id]['location']
            locations_array = np.array(locations)

            x_coord = locations_array[:, 0]
            y_coord = locations_array[:, 1]
            angle = locations_array[:, 2]

            # mark starting point
            plt.plot(x_coord[0], y_coord[0], 'r*')
            for i in range(len(locations)):
                plt.plot(x_coord[i:i + 2], y_coord[i:i + 2], 'b-')

        plt.xlim((0, edge_length))
        plt.ylim((0, edge_length))

        plt.savefig(output_dir + '/location_trace')
        plt.clf()