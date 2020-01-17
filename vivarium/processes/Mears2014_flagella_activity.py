from __future__ import absolute_import, division, print_function

import os
import numpy as np
import random
import math
import uuid

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch

from vivarium.actor.process import Process, deep_merge


DEFAULT_N_FLAGELLA = 5
DEFAULT_PMF = 170  # PMF ~170mV at pH 7, ~140mV at pH 7.7 (Berg)

DEFAULT_PARAMETERS = {
    # parameters from Mears, Koirala, Rao, Golding, Chemla (2014)
    'ccw_to_cw': 0.26,  # (1/s) Motor switching rate from CCW->CW
    'cw_to_ccw': 1.7,  # (1/s) Motor switching rate from CW->CCW
    'CB': 0.13,  # average CW bias of wild-type motors
    'lambda': 0.68,  # (1/s) transition rate from semi-coiled to curly-w state

    # CheY-P flucutations
    'YP_ss': 2.59,  # (uM) steady state concentration of CheY-P
    'sigma2_Y': 1.0,  # (uM^2) variance in CheY-P
    'tau': 0.2,  # (s) characteristic time-scale fluctuations in [CheY-P]

    # CW bias
    'K_d': 3.1,  # (uM) midpoint of CW bias vs CheY-P response curve
    'H': 10.3,  # Hill coefficient for CW bias vs CheY-P response curve

    # rotational state of individual flagella
    # parameters from Sneddon, Pontius, and Emonet (2012)
    'omega': 1.3,  # (1/s) characteristic motor switch time
    'g_0': 40,  # (k_B/T) free energy barrier for CCW-->CW
    'g_1': 40,  # (k_B/T) free energy barrier for CW-->CCW
    'K_D': 3.06,  # binding constant of Chey-P to base of the motor
}

# initial state
INITIAL_STATE = {
    'CheY': 2.59,
    'CheY_P': 2.59,  # (uM) mean concentration of CheY-P
    'cw_bias': 0.5,  # (made up)
    'motile_state': 0, # 1 for tumble, -1 for run, 0 for none
    'motile_force': 0,
    'motile_torque': 0,
}



def tumble(PMF):
    tumble_jitter = 0.4
    force = 0.006 * PMF  # 1.0
    torque = random.normalvariate(0, tumble_jitter)
    return [force, torque]

def run(PMF):
    force = 0.0124 * PMF  # 2.1
    torque = 0.0
    return [force, torque]

class FlagellaActivity(Process):
    '''
    Model of multi-flagellar motor activity and CheY-P fluctuations from:
        Mears, P. J., Koirala, S., Rao, C. V., Golding, I., & Chemla, Y. R. (2014).
        Escherichia coli swimming is robust against variations in flagellar number.

    TODO:
        - Complete PMF control over motor thrust.
        - Mears et al. has flagella with 3 conformational states for flagella (normal (CCW), semi (CW), curly (CW)).
        - Flagella states should be nested dictionaries with multiple states (rotational state, # of motors engaged).
        - If flagella counts are modified outside of this process, there needs to be a function that detects it and adjusts the flagella states.
        - Flagella will need to be separated upon division, rather than having each daughter inherit all the flagella.
    '''

    def __init__(self, initial_parameters={}):

        self.n_flagella = initial_parameters.get('n_flagella', DEFAULT_N_FLAGELLA)
        self.flagella_ids = [str(uuid.uuid1()) for flagella in range(self.n_flagella)]

        roles = {
            'internal': [
                'chemoreceptor_activity',
                'CheY',
                'CheY_P',
                'cw_bias',
                'motile_state',
                'motile_force',
                'motile_torque',
                'n_flagella'],
            'membrane': ['PMF', 'protons_flux_accumulated'],
            'flagella': self.flagella_ids,
            'external': []
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # flagella motor state: -1 for CCW, 1 for CW
        # motile state: -1 for run, 1 for tumble, 0 for no state
        internal = INITIAL_STATE
        default_state = {
            'external': {},
            'membrane': {'PMF': DEFAULT_PMF, 'PROTONS': 0},
            'flagella': {flagella_id: random.choice([-1, 1]) for flagella_id in self.flagella_ids},
            'internal': deep_merge(internal, {'volume': 1, 'n_flagella': self.n_flagella})}

        # default emitter keys
        default_emitter_keys = {
            'internal': [
                'n_flagella',
                'motile_force',
                'motile_torque',
                'motile_state'],
            'flagella': [],
            'external': [],
        }

        # default updaters
        internal_set_states = [
            'CheY',
            'CheY_P',
            'cw_bias',
            'motile_state',
            'motile_force',
            'motile_torque']

        default_updaters = {
            'internal': {state_id: 'set' for state_id in internal_set_states},
            'membrane': {'PROTONS': 'accumulate'},
            'flagella': {flagella_id: 'set' for flagella_id in self.flagella_ids},
            'external': {}}

        default_settings = {
            'process_id': 'motor',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'time_step': 0.001}

        return default_settings

    def next_update(self, timestep, states):

        internal = states['internal']
        n_flagella = states['internal']['n_flagella']
        flagella = states['flagella']
        PMF = states['membrane']['PMF']

        # states
        # TODO -- chemoreceptor?
        CheY = internal['CheY']
        CheY_P = internal['CheY_P']

        # parameters
        tau = self.parameters['tau']
        YP_ss = self.parameters['YP_ss']
        sigma = self.parameters['sigma2_Y']**0.5
        K_d = self.parameters['K_d']
        H = self.parameters['H']

        ## update CheY-P
        # Mears et al. eqn S6
        dYP = -(1 / tau) * (CheY_P - YP_ss) * timestep + sigma * (2 * timestep / tau)**0.5 * random.normalvariate(0, 1)
        CheY_P = max(CheY_P + dYP, 0.0)  # keep value positive
        CheY = max(CheY - dYP, 0.0)  # keep value positive

        ## CW bias, Hill function
        # Cluzel, P., Surette, M., & Leibler, S. (2000).
        cw_bias = CheY_P**H / (K_d**H + CheY_P**H)

        ## update all flagella
        flagella_update = {}
        for flagella_id, motor_state in flagella.items():
            new_motor_state = self.update_flagellum(motor_state, cw_bias, CheY_P, timestep)
            flagella_update.update({flagella_id: new_motor_state})

        ## get cell motile state.
        # if any flagella is rotating CW, the cell tumbles.
        # flagella motor state: -1 for CCW, 1 for CW
        # motile state: -1 for run, 1 for tumble, 0 for no state
        if any(state == 1 for state in flagella_update.values()):
            motile_state = 1
            [force, torque] = tumble(PMF)
        elif len(flagella_update) > 0:
            motile_state = -1
            [force, torque] = run(PMF)
        else:
            motile_state = 0
            force = 0
            torque = 0

        return {
            'flagella': flagella_update,
            'internal' : {
                'CheY': CheY,
                'CheY_P': CheY_P,
                'cw_bias': cw_bias,
                'motile_state': motile_state,
                'motile_force': force,
                'motile_torque': torque,
            }}

    def update_flagellum(self, motor_state, cw_bias, CheY_P, timestep):
        '''
        Rotational state of an individual flagellum from:
            Sneddon, M. W., Pontius, W., & Emonet, T. (2012).
            Stochastic coordination of multiple actuators reduces
            latency and improves chemotactic response in bacteria.

        # TODO -- normal, semi, curly states from Sneddon
        '''
        g_0 = self.parameters['g_0']  # (k_B/T) free energy barrier for CCW-->CW
        g_1 = self.parameters['g_1']  # (k_B/T) free energy barrier for CW-->CCW
        K_D = self.parameters['K_D']  # binding constant of Chey-P to base of the motor
        omega = self.parameters['omega']  # (1/s) characteristic motor switch time

        # free energy barrier
        delta_g = g_0 / 4 - g_1 / 2 * (CheY_P / (CheY_P + K_D))

        # switching frequency
        CW_to_CCW = omega * math.exp(delta_g)
        CCW_to_CW = omega * math.exp(-delta_g)
        # switch_freq = CCW_to_CW * (1 - cw_bias) + CW_to_CCW * cw_bias

        # flagella motor state: -1 for CCW, 1 for CW
        if motor_state == -1:
            prob_switch = CCW_to_CW * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 1
            else:
                new_motor_state = -1

        elif motor_state == 1:
            prob_switch = CW_to_CCW * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = -1
            else:
                new_motor_state = 1

        return new_motor_state



# testing functions
def test_activity(parameters, total_time=10):
    # TODO -- add asserts for test

    initial_params = {}
    initial_params.update(parameters)

    motor = FlagellaActivity(initial_params)
    settings = motor.default_settings()
    state = settings['state']

    receptor_activity = 1./3.
    state['internal']['chemoreceptor_activity'] = receptor_activity

    saved_data = {
        'internal': {state_id: [value] for state_id, value in state['internal'].items()},
        'flagella': {state_id: [value] for state_id, value in state['flagella'].items()},
        'time': [0]}

    # run simulation
    time = 0
    timestep = 0.001  # sec
    while time < total_time:
        time += timestep

        update = motor.next_update(timestep, state)

        # state is set
        state['flagella'].update(update['flagella'])
        state['internal'].update(update['internal'])

        saved_data['time'].append(time)
        for role in ['internal', 'flagella',]:
            for state_id, value in state[role].items():
                saved_data[role][state_id].append(value)

    return saved_data

def test_motor_PMF():
    from numpy import linspace

    # range of PMF value for test
    PMF_values = linspace(50.0, 200.0, 501).tolist()
    timestep = 1

    # initialize process and state
    motor = FlagellaActivity()
    settings = motor.default_settings()
    state = settings['state']

    motor_state_vec = []
    motor_force_vec = []
    motor_torque_vec = []
    for PMF in PMF_values:
        state['membrane']['PMF'] = PMF
        update = motor.next_update(timestep, state)

        motile_state = update['internal']['motile_state']
        motile_force = update['internal']['motile_force']
        motile_torque = update['internal']['motile_torque']

        # save
        motor_state_vec.append(motile_state)
        motor_force_vec.append(motile_force)
        motor_torque_vec.append(motile_torque)

    return {
        'motile_state': motor_state_vec,
        'motile_force': motor_force_vec,
        'motile_torque': motor_torque_vec,
        'PMF': PMF_values,
    }

def plot_activity(output, out_dir='out', filename='motor_control'):
    # receptor_activities = output['receptor_activities']
    CheY_vec = output['internal']['CheY']
    CheY_P_vec = output['internal']['CheY_P']
    cw_bias_vec = output['internal']['cw_bias']
    motile_state_vec = output['internal']['motile_state']
    flagella = output['flagella']
    time_vec = output['time']

    # get all flagella states into an activity grid
    flagella_ids = list(flagella.keys())
    activity_grid = np.zeros((len(flagella_ids), len(time_vec)))
    total_CW = np.zeros((len(time_vec)))
    for flagella_id, rotation_states in flagella.items():
        flagella_index = flagella_ids.index(flagella_id)
        modified_rotation_state = [1 if x==-1 else 2 if x==1 else 0 for x in rotation_states]
        activity_grid[flagella_index, :] = modified_rotation_state

        CW_rotation_state = [1 if x == 1 else 0 for x in rotation_states]
        total_CW += np.array(CW_rotation_state)

    # grid for cell state
    motile_state_grid = np.zeros((1, len(time_vec)))
    motile_state_grid[0, :] = motile_state_vec

    # set up colormaps
    # cell motile state
    cmap1 = colors.ListedColormap(['black', 'red', 'white'])
    bounds1 = [-1, -1/3, 1/3, 1]
    norm1 = colors.BoundaryNorm(bounds1, cmap1.N)
    motile_legend_elements = [
        Patch(facecolor='k', edgecolor='k', label='Run'),
        Patch(facecolor='r', edgecolor='k', label='None'),
        Patch(facecolor='w', edgecolor='k', label='Tumble')]

    # rotational state
    cmap2 = colors.ListedColormap(['black', 'white', 'blue'])
    bounds2 = [0, 0.5, 1.5, 2]
    norm2 = colors.BoundaryNorm(bounds2, cmap2.N)
    rotational_legend_elements = [
        Patch(facecolor='b', edgecolor='k', label='CW'),
        Patch(facecolor='w', edgecolor='k', label='CCW')]


    # plot results
    cols = 1
    rows = 5
    plt.figure(figsize=(10 * cols, 2 * rows))

    # define subplots
    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)
    ax4 = plt.subplot(rows, cols, 4)
    ax5 = plt.subplot(rows, cols, 5)

    # plot Che-P state
    ax1.plot(time_vec, CheY_vec, label='CheY')
    ax1.plot(time_vec, CheY_P_vec, label='CheY_P')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_xticks([])
    ax1.set_xlim(time_vec[0], time_vec[-1])
    ax1.set_ylabel('concentration (uM)')

    # plot CW bias
    ax2.plot(time_vec, cw_bias_vec)
    ax2.set_xticks([])
    ax2.set_xlim(time_vec[0], time_vec[-1])
    ax2.set_ylabel('CW bias')

    # plot cell motile state
    im1 = ax3.imshow(motile_state_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap1,
               norm=norm1)
    # cbar = plt.colorbar(im1, cmap=cmap1, norm=norm1, boundaries=bounds1, ticks=[0,1])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_ylabel('cell motile state')

    # legend
    ax3.legend(
        handles=motile_legend_elements,
        loc='center left',
        bbox_to_anchor=(1, 0.5))

    # plot flagella states in a grid
    if len(activity_grid) > 0:
        im2 = ax4.imshow(activity_grid,
                   interpolation='nearest',
                   aspect='auto',
                   cmap=cmap2,
                   norm=norm2)
        # cbar = plt.colorbar(im2, cmap=cmap2, norm=norm2, boundaries=bounds2, ticks=[0,1,2])
        # cbar.set_ticklabels(['none', 'CCW', 'CW'])
        plt.locator_params(axis='y', nbins=len(flagella_ids))
        ax4.set_yticks(list(range(len(flagella_ids))))
        ax4.set_xticks([])
        ax4.set_ylabel('flagella #')

        # legend
        ax4.legend(
            handles=rotational_legend_elements,
            loc='center left',
            bbox_to_anchor=(1, 0.5))
    else:
        # no flagella
        ax4.set_axis_off()

    # plot number of flagella CW
    ax5.plot(time_vec, total_CW)
    ax5.set_xlim(time_vec[0], time_vec[-1])
    ax5.set_xlabel('time (sec)')
    ax5.set_ylabel('number of flagella CW')

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')

def plot_motor_PMF(output, out_dir='out'):
    motile_state = output['motile_state']
    motile_force = output['motile_force']
    motile_torque = output['motile_torque']
    PMF = output['PMF']

    # plot results
    cols = 1
    rows = 1
    plt.figure(figsize=(10 * cols, 2 * rows))

    # define subplots
    ax1 = plt.subplot(rows, cols, 1)

    # plot motile_state
    ax1.plot(PMF, motile_force)
    ax1.set_xlabel('PMF (mV)')
    ax1.set_ylabel('force')

    # save figure
    fig_path = os.path.join(out_dir, 'motor_PMF')
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'Mears2014_flagella_activity')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    zero_flagella = {'n_flagella': 0}
    output1 = test_activity(zero_flagella, 10)
    plot_activity(output1, out_dir, 'motor_control_zero_flagella')

    five_flagella = {'n_flagella': 5}
    output2 = test_activity(five_flagella, 10)
    plot_activity(output2, out_dir, 'motor_control')

    output3 = test_motor_PMF()
    plot_motor_PMF(output3, out_dir)
