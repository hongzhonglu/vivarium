from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis, get_lattice, get_compartment

class ChemotaxisTrace(Analysis):
    def __init__(self):
        super(ChemotaxisTrace, self).__init__(analysis_type='environment')

    def requirements(self):
        return ['motor']

    def get_data(self, client, query):
        query.update({'type': 'lattice'})
        history_data = client.find(query)
        history_data.sort('time')
        lattice_history = get_lattice(history_data)

        query = {'type': 'compartment'}
        history_data = client.find(query)
        history_data.sort('time')
        compartment_history = get_compartment(history_data)

        data = {
            'lattice_history': lattice_history,
            'compartment_history': compartment_history
        }
        return data

    def analyze(self, experiment_config, data, output_dir):
        lattice_history = data['lattice_history']
        # motor_state = data['compartment_history']['cell']['motor_state']

        agent_ids = [key for key in lattice_history.keys() if key is not 'time']
        time_vec = [t / 3600 for t in lattice_history['time']]  # convert to hours
        edge_length_x = experiment_config['edge_length_x']
        edge_length_y = experiment_config['edge_length_y']

        plt.figure(figsize=(8, 8))
        for agent_id in agent_ids:
            # get locations and convert to 2D array
            locations = lattice_history[agent_id]['location']
            locations_array = np.array(locations)

            x_coord = locations_array[:, 0]
            y_coord = locations_array[:, 1]
            angle = locations_array[:, 2]

            # mark starting point
            plt.plot(x_coord[0], y_coord[0], 'r*')
            for i in range(len(locations)):
                plt.plot(x_coord[i:i + 2], y_coord[i:i + 2], 'b-')

        plt.xlim((0, edge_length_x))
        plt.ylim((0, edge_length_y))

        plt.savefig(output_dir + '/chemotaxis_trace')
        plt.clf()