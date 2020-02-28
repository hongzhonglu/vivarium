from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from vivarium.analysis.analysis import Analysis, get_compartment, get_lattice

class Motor(Analysis):
    def __init__(self):
        super(Motor, self).__init__(analysis_type='compartment_with_env')

    def requirements(self):
        return ['motor']

    def get_data(self, client, query, options={}):

        type = options.get('type')
        if type is 'compartment':
            sim_id = query['simulation_id']
            query.update({'type': 'compartment'})
            history_data = client.find(query)
            history_data.sort('time')
            compartment_history = get_compartment(history_data)
            compartment_history['sim_id'] = sim_id  # put sim_id back into data

            return compartment_history

        elif type is 'environment':
            query.update({'type': 'lattice-agent'})
            history_data = client.find(query)
            history_data.sort('time')
            lattice_history = get_lattice(history_data)

            return lattice_history

    def analyze(self, experiment_config, data, output_dir):

        expected_speed = 14.2  # um/s (Berg)
        expected_run_chemotax = 0.86  # s (Berg) expected run length when chemotaxis
        expected_run = 0.42  # s (Berg) expected run length without chemotaxis
        expected_tumble = 0.14  # s (Berg)
        expected_angle_between_runs = 68  # degrees (Berg)

        # compartment data
        compartment_data = data['compartment']
        sim_id = compartment_data['sim_id']
        time_vec = compartment_data['time']  # convert to hours

        # get the internal port. assumes there is only one in all agents.
        internal_port = None
        for agent_id, specs in experiment_config['agents'].items():
            internal_port = specs['topology']['motor']['internal']

        # TODO -- why can internal_port be incorrect in topology?
        if internal_port not in compartment_data:
            internal_port = 'cell'

        CheY_P_vec = compartment_data[internal_port]['CheY_P']
        ccw_motor_bias_vec = compartment_data[internal_port]['ccw_motor_bias']
        ccw_to_cw_vec = compartment_data[internal_port]['ccw_to_cw']
        motor_state_vec = compartment_data[internal_port]['motor_state']

        # environment data for this sim
        env_time_vec = data['environment']['time']  # seconds
        environment_data = data['environment'][sim_id]
        location_vec = environment_data['location']

        # get speed
        speed_vec = [0]
        previous_time = env_time_vec[0]
        previous_loc = location_vec[0]
        previous_motor_state = motor_state_vec[0]
        run_angle = 0
        angle_between_runs = []
        for time, location, m_state in zip(env_time_vec[1:], location_vec[1:], motor_state_vec[1:]):
            dt = time - previous_time
            distance = ((location[0] - previous_loc[0])**2 + (location[1] - previous_loc[1])**2)**0.5
            speed_vec.append(distance/dt)  # um/sec

            # get change in angles between motor states
            if m_state == 0 and previous_motor_state == 1:
                angle_change = abs(location[2] - run_angle) / np.pi * 180 % 360  # convert to absolute degrees
                angle_between_runs.append(angle_change)
            elif m_state == 0:
                run_angle = location[2]

            # update previous states
            previous_time = time
            previous_loc = location
            previous_motor_state = m_state

        avg_speed = sum(speed_vec) / len(speed_vec)
        avg_angle_between_runs = sum(angle_between_runs) / len(angle_between_runs)

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

        # make figure
        n_cols = 1
        n_rows = 5
        fig = plt.figure(figsize=(6 * n_cols, 2 * n_rows))
        plt.rcParams.update({'font.size': 10})
        fig.suptitle('{}'.format(sim_id))

        # define subplots
        ax1 = plt.subplot(n_rows, n_cols, 1)
        ax2 = plt.subplot(n_rows, n_cols, 2)
        ax3 = plt.subplot(n_rows, n_cols, 3)
        ax4 = plt.subplot(n_rows, n_cols, 4)
        ax5 = plt.subplot(n_rows, n_cols, 5)

        # plot data
        ax1.plot(CheY_P_vec, 'b')
        ax2.plot(ccw_motor_bias_vec, 'b', label='ccw_motor_bias')
        ax2.plot(ccw_to_cw_vec, 'g', label='ccw_to_cw')
        ax3.plot(speed_vec)
        ax3.axhline(y=avg_speed, color='r', linestyle='dashed', label='mean')
        ax3.axhline(y=expected_speed, color='b', linestyle='dashed', label='expected mean')

        # plot change in direction between runs
        ax4.plot(angle_between_runs)
        ax4.axhline(y=avg_angle_between_runs, color='b', linestyle='dashed', label='mean angle between runs')
        ax4.axhline(y=expected_angle_between_runs, color='r', linestyle='dashed', label='exp. angle between runs')

        # plot run/tumble distributions
        max_length = max(run_lengths + tumble_lengths)
        bins = np.linspace(0, max_length, 10)
        logbins = np.logspace(0, np.log10(bins[-1]), len(bins))
        ax5.hist([run_lengths, tumble_lengths], bins=logbins, label=['run_lengths', 'tumble_lengths'], color=['b', 'm'])

        # plot expected values
        ax5.axvline(x=expected_tumble, color='m', linestyle='dashed', label='expected tumble')
        ax5.axvline(x=expected_run, color='b', linestyle='dashed', label='expected run')
        ax5.axvline(x=expected_run_chemotax, color='c', linestyle='dashed', label='expected chemotaxis run')


        # labels
        ax1.set_xticklabels([])
        ax1.set_ylabel("CheY_P")
        ax2.set_xticklabels([])
        ax2.set_ylabel("motor bias")
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax3.set_ylabel(u'speed (\u03bcm/sec)')
        ax3.set_xlabel('time')
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax4.set_ylabel('angle between runs')
        ax4.set_xlabel('run #')
        ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax5.set_xlabel("motor state length (sec)")
        ax5.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax5.set_xscale('log')


        plt.savefig(output_dir + '/motor', bbox_inches='tight')
        plt.close(fig)
