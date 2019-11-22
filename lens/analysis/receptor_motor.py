from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis, get_compartment

class ReceptorMotor(Analysis):
    def __init__(self):
        super(ReceptorMotor, self).__init__(analysis_type='compartment')

    def requirements(self):
        return ['receptor', 'motor']

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
        receptor_activities = data['cell']['chemoreceptor_activity']
        CheY_P_vec = data['cell']['CheY_P']
        ccw_motor_bias_vec = data['cell']['ccw_motor_bias']
        ccw_to_cw_vec = data['cell']['ccw_to_cw']
        motor_state_vec = data['cell']['motor_state']

        # plot results
        cols = 1
        rows = 4
        fig = plt.figure(figsize=(6 * cols, 1 * rows))

        plt.savefig(output_dir + '/receptor_motor', bbox_inches='tight')
        plt.close(fig)
