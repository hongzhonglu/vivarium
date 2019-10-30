from __future__ import absolute_import, division, print_function

import os
import numpy as np
import random

from lens.actor.process import Process, dict_merge


TUMBLE_JITTER = 1.0  # (radians)

# parameters
DEFAULT_PARAMETERS = {
    # 'k_A': 5.0,  #
    'k_y': 100.0,  # 1/uM/s
    'k_z': 30.0,  # / CheZ,
    'gamma_Y': 0.1,
    'k_s': 0.45,  # scaling coefficient
    'adaptPrecision': 3,
    # motor
    'mb_0': 0.65,  # steady state motor bias (Cluzel et al 2000)
    'n_motors': 5,
    'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
}

##initial state
INITIAL_STATE = {
    # response regulator proteins
    'CheY_tot': 9.7,  # (uM) #0.0097,  # (mM) 9.7 uM = 0.0097 mM
    'CheY_P': 0.0,
    'CheZ': 0.01*100,  # (uM) #phosphatase 100 uM = 0.1 mM (0.01 scaling from RapidCell1.4.2)
    'CheA': 0.01*100,  # (uM) #100 uM = 0.1 mM (0.01 scaling from RapidCell1.4.2)
    # sensor activity
    'chemoreceptor_activity': 1/3,
    # motor activity
    'motile_force': 0,
    'motile_torque': 0,
    'motor_state': 1,  # motor_state 1 for tumble, 0 for run
}

class MotorActivity(Process):
    '''
    Model of motor activity from:
        Vladimirov, N., Lovdok, L., Lebiedz, D., & Sourjik, V. (2008).
        Dependence of bacterial chemotaxis on gradient shape and adaptation rate.
    '''
    def __init__(self, initial_parameters={}):

        roles = {
            'internal': ['chemoreceptor_activity',
                         'CheA',
                         'CheZ',
                         'CheY_tot',
                         'CheY_P',
                         'ccw_motor_bias',
                         'ccw_to_cw',
                         'motile_force',
                         'motile_torque',
                         'motor_state'],
            'external': []
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(MotorActivity, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''
        internal = INITIAL_STATE
        return {
            'external': {},
            'internal': dict_merge(internal, {'volume': 1})}

    def default_emitter_keys(self):
        keys = {
            'internal': [
                'ccw_motor_bias',
                'ccw_to_cw',
                'motile_force',
                'motile_torque',
                'motor_state',
                'CheA',
                'CheY_P'],
            'external': [],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''
        set_states = [
            'ccw_motor_bias',
            'ccw_to_cw',
            'motile_force',
            'motile_torque',
            'motor_state',
            'CheA',
            'CheY_P']
        updater_types = {
            'internal': {state_id: 'set' for state_id in set_states},
            'external': {}}

        return updater_types

    def next_update(self, timestep, states):
        '''
        CheY phosphorylation model from:
            Kollmann, M., Lovdok, L., Bartholome, K., Timmer, J., & Sourjik, V. (2005).
            Design principles of a bacterial signalling network. Nature.
        Motor switching model from:
            Scharf, B. E., Fahrner, K. A., Turner, L., and Berg, H. C. (1998).
            Control of direction of flagellar rotation in bacterial chemotaxis. PNAS.

        An increase of attractant inhibits CheA activity (chemoreceptor_activity),
        but subsequent methylation returns CheA activity to its original level.
        TODO -- add CheB phosphorylation
        '''

        internal = states['internal']
        P_on = internal['chemoreceptor_activity']
        motor_state = internal['motor_state']

        # parameters
        adaptPrecision = self.parameters['adaptPrecision']
        k_y = self.parameters['k_y']
        k_s = self.parameters['k_s']
        k_z = self.parameters['k_z']
        gamma_Y  =self.parameters['gamma_Y']
        mb_0 = self.parameters['mb_0']
        cw_to_ccw = self.parameters['cw_to_ccw']

        ## Kinase activity
        # relative steady-state concentration of phosphorylated CheY. Assumes Che_P = P_on
        scaling = 1.  # 1.66889  # 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        CheY_P = adaptPrecision * scaling * k_y * k_s * P_on / (k_y * k_s * P_on + k_z + gamma_Y)  # CheZ cancels out of k_z

        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = mb_0 / (CheY_P * (1 - mb_0) + mb_0)
        ccw_to_cw = cw_to_ccw * (1 / ccw_motor_bias - 1)

        if motor_state == 0:  # 0 for run
            # switch to tumble?
            prob_switch = ccw_to_cw * timestep
            if np.random.random(1)[0] <= prob_switch:
                motor_state = 1
                force, torque = self.tumble()
            else:
                force, torque = self.run()

        elif motor_state == 1:  # 1 for tumble
            # switch to run?
            prob_switch = cw_to_ccw * timestep
            if np.random.random(1)[0] <= prob_switch:
                motor_state = 0
                force, torque = self.run()
            else:
                force, torque = self.tumble()

        # TODO -- should force/torque accumulate over exchange timestep?
        update = {
            'internal': {
                'ccw_motor_bias': ccw_motor_bias,
                'ccw_to_cw': ccw_to_cw,
                'motile_force': force,
                'motile_torque': torque,
                'motor_state': motor_state,
                'CheY_P': CheY_P
            }
        }
        return update

    def tumble(self):
        force = 2.0  # 5.0
        torque = random.normalvariate(0, TUMBLE_JITTER)
        return force, torque

    def run(self):
        force = 4.0  # 15.0
        torque = 0.0
        return force, torque


def test_motor_control():
    # TODO -- add asserts for test

    motor = MotorActivity()
    state = motor.default_state()
    receptor_activity = 1./3.
    state['internal']['chemoreceptor_activity'] = receptor_activity

    CheY_P_vec = []
    ccw_motor_bias_vec = []
    ccw_to_cw_vec = []
    motor_state_vec = []
    time_vec = []

    # run simulation
    time = 0
    timestep = 0.01  # sec
    end_time = 360  # secs
    while time < end_time:

        update = motor.next_update(timestep, state)
        CheY_P = update['internal']['CheY_P']
        ccw_motor_bias = update['internal']['ccw_motor_bias']
        ccw_to_cw = update['internal']['ccw_to_cw']
        motor_state = update['internal']['motor_state']

        # update motor state
        state['internal']['motor_state'] = motor_state

        CheY_P_vec.append(CheY_P)
        ccw_motor_bias_vec.append(ccw_motor_bias)
        ccw_to_cw_vec.append(ccw_to_cw)
        motor_state_vec.append(motor_state)
        time_vec.append(time)

        time += timestep

    return {
        'CheY_P_vec': CheY_P_vec,
        'ccw_motor_bias_vec': ccw_motor_bias_vec,
        'ccw_to_cw_vec': ccw_to_cw_vec,
        'motor_state_vec': motor_state_vec,
        'time_vec': time_vec}

def test_variable_receptor():
    from numpy import linspace

    motor = MotorActivity()
    state = motor.default_state()
    timestep = 1
    receptor_activities = linspace(0.0, 1.0, 501).tolist()
    CheY_P_vec = []
    ccw_motor_bias_vec = []
    ccw_to_cw_vec = []
    motor_state_vec = []
    for activity in receptor_activities:
        state['internal']['chemoreceptor_activity'] = activity
        update = motor.next_update(timestep, state)
        CheY_P = update['internal']['CheY_P']
        ccw_motor_bias = update['internal']['ccw_motor_bias']
        ccw_to_cw = update['internal']['ccw_to_cw']
        motile_state = update['internal']['motor_state']

        CheY_P_vec.append(CheY_P)
        ccw_motor_bias_vec.append(ccw_motor_bias)
        ccw_to_cw_vec.append(ccw_to_cw)
        motor_state_vec.append(motile_state)

    # check ccw_to_cw bias is strictly increasing with increased receptor activity
    assert all(i < j for i, j in zip(ccw_to_cw_vec, ccw_to_cw_vec[1:]))

    return {
        'receptor_activities': receptor_activities,
        'CheY_P_vec': CheY_P_vec,
        'ccw_motor_bias_vec': ccw_motor_bias_vec,
        'ccw_to_cw_vec': ccw_to_cw_vec,
        'motor_state_vec': motor_state_vec}

def plot_motor_control(output, out_dir='out'):
    # TODO -- make this into an analysis figure
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    # receptor_activities = output['receptor_activities']
    CheY_P_vec = output['CheY_P_vec']
    ccw_motor_bias_vec = output['ccw_motor_bias_vec']
    ccw_to_cw_vec = output['ccw_to_cw_vec']
    motor_state_vec = output['motor_state_vec']
    time_vec = output['time_vec']

    # plot results
    cols = 1
    rows = 4
    plt.figure(figsize=(6 * cols, 1 * rows))

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
    for state, time in zip(motor_state_vec,time_vec):
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

    fig_path = os.path.join(out_dir, 'motor_control')
    plt.subplots_adjust(wspace=0.7, hspace=0.5)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


def plot_variable_receptor(output, out_dir='out'):
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    receptor_activities = output['receptor_activities']
    CheY_P_vec = output['CheY_P_vec']
    ccw_motor_bias_vec = output['ccw_motor_bias_vec']
    ccw_to_cw_vec = output['ccw_to_cw_vec']
    motor_state_vec = output['motor_state_vec']

    # plot results
    cols = 1
    rows = 4
    plt.figure(figsize=(6 * cols, 1 * rows))

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)
    ax4 = plt.subplot(rows, cols, 4)

    ax1.plot(receptor_activities, 'b')
    ax2.plot(CheY_P_vec, 'b')
    ax3.plot(ccw_motor_bias_vec, 'b', label='ccw_motor_bias')
    ax3.plot(ccw_to_cw_vec, 'g', label='ccw_to_cw')
    ax4.plot(motor_state_vec, '.b')

    ax1.set_xticklabels([])
    ax1.set_ylabel("receptor activity \n P(on) ", fontsize=10)
    ax2.set_xticklabels([])
    ax2.set_ylabel("CheY_P", fontsize=10)
    ax3.set_xticklabels([])
    ax3.set_ylabel("motor bias", fontsize=10)
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_yticks([0.0, 1.0])
    ax4.set_yticklabels(["run", "tumble"])

    fig_path = os.path.join(out_dir, 'motor_variable_receptor')
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'Vladimirov2008_motor')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output1 = test_motor_control()
    plot_motor_control(output1, out_dir)

    output2 = test_variable_receptor()
    plot_variable_receptor(output2, out_dir)
