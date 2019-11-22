from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis, get_compartment

class Motor(Analysis):
    def __init__(self):
        super(Motor, self).__init__(analysis_type='compartment')

    def requirements(self):
        return ['motor']

    def get_data(self, client, query, options={}):

        sim_id = query['simulation_id']
        query.update({'type': 'compartment'})
        history_data = client.find(query)
        history_data.sort('time')
        compartment_history = get_compartment(history_data)
        # put sim_id back into data
        compartment_history['sim_id'] = sim_id

        return compartment_history


    def analyze(self, experiment_config, data, output_dir):
        # skip_keys = ['time', 'sim_id']
        time_vec = [t / 3600 for t in data['time']]  # convert to hours
        CheY_P_vec = data['cell']['CheY_P']
        ccw_motor_bias_vec = data['cell']['ccw_motor_bias']
        ccw_to_cw_vec = data['cell']['ccw_to_cw']
        motor_state_vec = data['cell']['motor_state']


        # plot results
        cols = 1
        rows = 4
        fig = plt.figure(figsize=(6 * cols, 1 * rows))

        ax1 = plt.subplot(rows, cols, 1)
        ax2 = plt.subplot(rows, cols, 2)
        ax4 = plt.subplot(rows, cols, 3)

        ax1.plot(CheY_P_vec, 'b')
        ax2.plot(ccw_motor_bias_vec, 'b', label='ccw_motor_bias')
        ax2.plot(ccw_to_cw_vec, 'g', label='ccw_to_cw')

        # get length of runs, tumbles
        run_lengths = []
        tumble_lengths = []
        prior_state = 0
        state_start_time = 0
        for state, time in zip(motor_state_vec, time_vec):
            if state == 0:  # run
                if prior_state != 0:
                    tumble_lengths.append(time - state_start_time)
                    state_start_time = time
            elif state == 1:  # tumble
                if prior_state != 1:
                    run_lengths.append(time - state_start_time)
                    state_start_time = time
            prior_state = state

        # plot run/tumble distributions
        max_length = max(run_lengths + tumble_lengths)
        bins = np.linspace(0, max_length, 20)
        ax4.hist([run_lengths, tumble_lengths], bins, label=['run_lengths', 'tumble_lengths'])

        # labels
        ax1.set_xticklabels([])
        ax1.set_ylabel("CheY_P", fontsize=10)
        ax2.set_xticklabels([])
        ax2.set_ylabel("motor bias", fontsize=10)
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax4.set_xlabel("motor state length (sec)", fontsize=10)

        plt.savefig(output_dir + '/motor', bbox_inches='tight')
        plt.close(fig)
