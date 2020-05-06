from __future__ import absolute_import, division, print_function

import os
import random
import math
import uuid
import argparse

import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch

from vivarium.compartment.process import Process
from vivarium.compartment.composition import simulate_process_with_environment


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

    defaults = {
        'flagella': 5,
        'parameters': DEFAULT_PARAMETERS,
        'initial_state': INITIAL_STATE,

        # motile force parameters
        'PMF': 170,  # PMF ~170mV at pH 7, ~140mV at pH 7.7 (Berg H, E. coli in motion, 2004, pg 113)
        'flagellum_thrust': 25,  # (pN) (Berg H, E. coli in motion, 2004, pg 113)
        'tumble_jitter': 0.4,
    }

    def __init__(self, initial_parameters={}):

        self.n_flagella = initial_parameters.get('flagella', self.defaults['flagella'])
        self.flagellum_thrust = initial_parameters.get('flagellum_thrust', self.defaults['flagellum_thrust'])
        self.tumble_jitter = initial_parameters.get('tumble_jitter', self.defaults['tumble_jitter'])
        self.tumble_scaling = 0.5/self.defaults['PMF']
        self.run_scaling = 1 / self.defaults['PMF']


        self.flagella_ids = [str(uuid.uuid1()) for flagella in range(self.n_flagella)]

        ports = {
            'internal': [
                'chemoreceptor_activity',
                'CheY',
                'CheY_P',
                'cw_bias',
                'motile_state',
                'motile_force',
                'motile_torque'],
            'flagella_counts': [
                'flagella'],
            'membrane': [
                'PMF',
                'protons_flux_accumulated'],
            'flagella_activity': [
                'flagella'],
            'external': []}

        parameters = self.defaults['parameters']
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(ports, parameters)

    def default_settings(self):

        # default state
        # flagella motor state: -1 for CCW, 1 for CW
        # motile state: -1 for run, 1 for tumble, 0 for no state
        internal = self.defaults['initial_state']
        default_state = {
            'external': {},
            'membrane': {
                'PMF': self.defaults['PMF'],
                'PROTONS': 0},
            'flagella_counts': {
                'flagella': self.n_flagella},
            'flagella_activity': {
                'flagella': {
                    flagella_id: random.choice([-1, 1]) for flagella_id in self.flagella_ids}},
            'internal': internal}

        # default emitter keys
        default_emitter_keys = {
            'internal': ['motile_force', 'motile_torque', 'motile_state', 'CheY', 'CheY_P', 'cw_bias'],
            'flagella_counts': ['flagella'],
            'flagella_activity': ['flagella']}

        # default updaters
        internal_set_states = [
            'CheY',
            'CheY_P',
            'cw_bias',
            'motile_state',
            'motile_force',
            'motile_torque']

        # schema
        internal_schema = {
            state_id: {
                'updater': 'set',
                'divide': 'set'}
            for state_id in internal_set_states}
        schema = {
            'internal': internal_schema,
            'membrane': {'PROTONS': {'updater': 'accumulate'}},
            'flagella_activity': {'flagella': {
                'updater': 'set',
                'divide': 'split_dict'}}}

        default_settings = {
            'process_id': 'motor',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'time_step': 0.01}  # 0.001

        return default_settings

    def next_update(self, timestep, states):

        internal = states['internal']
        n_flagella = states['flagella_counts']['flagella']
        flagella = states['flagella_activity']['flagella']
        PMF = states['membrane']['PMF']

        # adjust number of flagella
        new_flagella = int(n_flagella) - len(flagella)
        if new_flagella < 0:
            # remove flagella
            remove = random.sample(self.flagella_ids, abs(new_flagella))
            for flg_id in remove:
                self.flagella_ids.remove(flg_id)
                del flagella[flg_id]

        elif new_flagella > 0:
            # add flagella
            new_flagella_ids = [str(uuid.uuid1())
                for flagella in range(new_flagella)]
            self.flagella_ids.extend(new_flagella_ids)
            new_flagella_states = {flg_id: random.choice([-1, 1])
                for flg_id in new_flagella_ids}
            flagella.update(new_flagella_states)

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
            [force, torque] = self.tumble(n_flagella, PMF)
        elif len(flagella_update) > 0:
            motile_state = -1
            [force, torque] = self.run(n_flagella, PMF)
        else:
            motile_state = 0
            force = 0
            torque = 0

        return {
            'flagella_activity': {
                'flagella': flagella_update},
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

    def tumble(self, n_flagella, PMF):
        force = self.tumble_scaling * PMF * self.flagellum_thrust * n_flagella
        torque = random.normalvariate(0, self.tumble_jitter)
        return [force, torque]

    def run(self, n_flagella, PMF):
        force = self.run_scaling * PMF * self.flagellum_thrust * n_flagella
        torque = 0.0
        return [force, torque]


# testing functions
default_params = {'flagella': 5}
default_timeline = [(10, {})]
def test_activity(parameters=default_params, timeline=default_timeline):
    motor = FlagellaActivity(parameters)

    settings = {
        'timeline': timeline,
        'environment_port': 'external',
    }

    return simulate_process_with_environment(motor, settings)

def test_motor_PMF():

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
    motile_force_vec = output['internal']['motile_force']
    flagella_activity = output['flagella_activity']['flagella']
    time_vec = output['time']

    # get flagella ids by order appearance
    flagella_ids = []
    for state in flagella_activity:
        flg_ids = list(state.keys())
        for flg_id in flg_ids:
            if flg_id not in flagella_ids:
                flagella_ids.append(flg_id)

    # make flagella activity grid
    activity_grid = np.zeros((len(flagella_ids), len(time_vec)))
    total_CW = np.zeros((len(time_vec)))
    for time_index, flagella_state in enumerate(flagella_activity):
        for flagella_id, rotation_states in flagella_state.items():

            # get this flagella's index
            flagella_index = flagella_ids.index(flagella_id)

            modified_rotation_state = 0
            CW_rotation_state = 0
            if rotation_states == -1:
                modified_rotation_state = 1
            elif rotation_states == 1:
                modified_rotation_state = 2
                CW_rotation_state = 1

            activity_grid[flagella_index, time_index] = modified_rotation_state
            total_CW += np.array(CW_rotation_state)

    # grid for cell state
    motile_state_grid = np.zeros((1, len(time_vec)))
    motile_state_grid[0, :] = motile_state_vec

    # set up colormaps
    # cell motile state
    cmap1 = colors.ListedColormap(['steelblue', 'lightgray', 'darkorange'])
    bounds1 = [-1, -1/3, 1/3, 1]
    norm1 = colors.BoundaryNorm(bounds1, cmap1.N)
    motile_legend_elements = [
        Patch(facecolor='steelblue', edgecolor='k', label='Run'),
        Patch(facecolor='darkorange', edgecolor='k', label='Tumble'),
        Patch(facecolor='lightgray', edgecolor='k', label='N/A')]

    # rotational state
    cmap2 = colors.ListedColormap(['lightgray', 'steelblue', 'darkorange'])
    bounds2 = [0, 0.5, 1.5, 2]
    norm2 = colors.BoundaryNorm(bounds2, cmap2.N)
    rotational_legend_elements = [
        Patch(facecolor='steelblue', edgecolor='k', label='CCW'),
        Patch(facecolor='darkorange', edgecolor='k', label='CW'),
        Patch(facecolor='lightgray', edgecolor='k', label='N/A')]

    # plot results
    cols = 1
    rows = 5
    plt.figure(figsize=(4 * cols, 1.5 * rows))

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
    ax1.set_ylabel('concentration \n (uM)')

    # plot CW bias
    ax2.plot(time_vec, cw_bias_vec)
    ax2.set_xticks([])
    ax2.set_xlim(time_vec[0], time_vec[-1])
    ax2.set_ylabel('CW bias')

    # plot flagella states in a grid
    if len(activity_grid) > 0:
        ax3.imshow(activity_grid,
                   interpolation='nearest',
                   aspect='auto',
                   cmap=cmap2,
                   norm=norm2,
                   # extent=[-1,1,-1,1]
                   extent=[time_vec[0], time_vec[-1], len(flagella_ids)+0.5, 0.5]
                   )
        plt.locator_params(axis='y', nbins=len(flagella_ids))
        ax3.set_yticks(list(range(1, len(flagella_ids) + 1)))
        ax3.set_xticks([])
        ax3.set_ylabel('flagella #')

        # legend
        ax3.legend(
            handles=rotational_legend_elements,
            loc='center left',
            bbox_to_anchor=(1, 0.5))
    else:
        # no flagella
        ax3.set_axis_off()

    # plot cell motile state
    ax4.imshow(motile_state_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap1,
               norm=norm1,
               extent=[time_vec[0], time_vec[-1], 0, 1])
    ax4.set_yticks([])
    ax4.set_xticks([])
    ax4.set_ylabel('cell motile state')

    # legend
    ax4.legend(
        handles=motile_legend_elements,
        loc='center left',
        bbox_to_anchor=(1, 0.5))

    # plot motor thrust
    ax5.plot(time_vec, motile_force_vec)
    ax5.set_xlim(time_vec[0], time_vec[-1])
    ax5.set_ylabel('total motor thrust (pN)')
    ax5.set_xlabel('time (sec)')


    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')

def plot_motor_PMF(output, out_dir='out', figname='motor_PMF'):
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
    fig_path = os.path.join(out_dir, figname)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')

def run_variable_flagella(out_dir):
    # variable flagella
    init_params = {'flagella': 5}
    timeline = [
        (0, {}),
        (60, {
            'flagella_counts': {
                'flagella': 6}}),
        (200, {
            'flagella_counts': {
                'flagella': 7}}),
        (240, {})]
    output3 = test_activity(init_params, timeline)
    plot_activity(output3, out_dir, 'variable_flagella')

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'Mears2014_flagella_activity')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='flagella expression')
    parser.add_argument('--variable', '-v', action='store_true', default=False,)
    args = parser.parse_args()

    if args.variable:
        run_variable_flagella(out_dir)
    else:
        zero_flagella = {'flagella': 0}
        timeline = [(10, {})]
        output1 = test_activity(zero_flagella, timeline)
        plot_activity(output1, out_dir, 'motor_control_zero_flagella')

        five_flagella = {'flagella': 5}
        timeline = [(10, {})]
        output2 = test_activity(five_flagella, timeline)
        plot_activity(output2, out_dir, 'motor_control')

        run_variable_flagella(out_dir)


