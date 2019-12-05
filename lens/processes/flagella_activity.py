from __future__ import absolute_import, division, print_function

import os
import copy
import numpy as np
import random
import uuid

from lens.actor.process import Process, deep_merge


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
    # 'n_motors': 5,
    'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
}

##initial state
DEFAULT_N_FLAGELLA = 10

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

class FlagellaActivity(Process):
    '''
    Model of motor activity from:
        Vladimirov, N., Lovdok, L., Lebiedz, D., & Sourjik, V. (2008).
        Dependence of bacterial chemotaxis on gradient shape and adaptation rate.
    '''
    def __init__(self, initial_parameters={}):

        self.n_flagella = initial_parameters.get('n_flagella', DEFAULT_N_FLAGELLA)
        self.flagella_ids = [str(uuid.uuid1()) for flagella in range(self.n_flagella)]

        roles = {
            'internal': ['n_flagella',
                         'chemoreceptor_activity',
                         'CheA',
                         'CheZ',
                         'CheY_tot',
                         'CheY_P',
                         # 'ccw_motor_bias',
                         # 'ccw_to_cw',
                         'motile_force',
                         'motile_torque',
                         'motor_state'],
            'membrane': ['PMF', 'PROTONS'],  # use PROTONS to count accumulated flux
            'flagella': self.flagella_ids,
            'external': []
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # flagella motor state: 0 for CCW, 1 for CW
        # PMF range 180-200 mV
        internal = INITIAL_STATE
        default_state = {
            'external': {},
            'membrane': {'PMF': 180, 'PROTONS': 0},
            'flagella': {flagella_id: random.choice([0, 1]) for flagella_id in self.flagella_ids},
            'internal': deep_merge(internal, {
                'volume': 1,
                'n_flagella': DEFAULT_N_FLAGELLA})}

        # default emitter keys
        default_emitter_keys = {
            'internal': [
                'motile_force',
                'motile_torque',
                'motor_state',
                'CheA',
                'CheY_P'],
            'flagella': [
                'ccw_motor_bias',
                'ccw_to_cw'],
            'external': [],
        }

        # default updaters
        internal_set_states = [
            'motile_force',
            'motile_torque',
            'motor_state',
            'CheA',
            'CheY_P']
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
            'time_step': 0.01}

        return default_settings

    def next_update(self, timestep, states):

        internal = states['internal']
        n_flagella = states['internal']['n_flagella']
        flagella = states['flagella']
        membrane = states['membrane']
        # TODO add new flagella if n_flagella > len(flagella)

        # determine behavior from motor states of all flagella
        flagella_update = {}
        for flagella_id, motor_state in flagella.items():
            new_motor_state = self.update_flagellum(internal, motor_state, timestep)
            flagella_update.update({flagella_id: new_motor_state})

        # TODO -- if all CCW corresponds then run.
        # TODO -- force and torque from motor state, linear with PMF
        fraction_CCW = flagella_update.values().count(0) / n_flagella

        if fraction_CCW > 0.8:
            force, torque = run()
        else:
            force, torque = tumble()


        CheY_P = 0.0
        # TODO -- should force/torque accumulate over exchange timestep?
        return {
            'flagella': flagella_update,
            'internal': {
                'motile_force': force,
                'motile_torque': torque,
                'CheY_P': CheY_P}
        }


    def update_flagellum(self, internal, motor_state, timestep):
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

        P_on = internal['chemoreceptor_activity']

        # parameters
        adaptPrecision = self.parameters['adaptPrecision']
        k_y = self.parameters['k_y']
        k_s = self.parameters['k_s']
        k_z = self.parameters['k_z']
        gamma_Y = self.parameters['gamma_Y']
        mb_0 = self.parameters['mb_0']
        cw_to_ccw = self.parameters['cw_to_ccw']

        ## Kinase activity
        # relative steady-state concentration of phosphorylated CheY.
        scaling = 1.  # 1.66889  # 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        CheY_P = adaptPrecision * scaling * k_y * k_s * P_on / (
                    k_y * k_s * P_on + k_z + gamma_Y)  # CheZ cancels out of k_z

        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = mb_0 / (CheY_P * (1 - mb_0) + mb_0)  # (1/s)
        ccw_to_cw = cw_to_ccw * (1 / ccw_motor_bias - 1)  # (1/s)

        new_motor_state = copy.copy(motor_state)
        if motor_state == 0:  # 0 for run
            # switch to tumble?
            prob_switch = ccw_to_cw * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 1

        elif motor_state == 1:  # 1 for tumble
            # switch to run?
            prob_switch = cw_to_ccw * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 0

        return new_motor_state


def tumble():
    tumble_jitter = 0.4  # (radians)  # TODO -- put in parameters
    force = 1.0
    torque = random.normalvariate(0, tumble_jitter)
    return [force, torque]

def run():
    force = 2.1
    torque = 0.0
    return [force, torque]


# testing functions
def test_motor_control(total_time=10):
    # TODO -- add asserts for test

    initial_params = {
        'mb_0': 0.65,  # steady state motor bias (Cluzel et al 2000)
        'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
    }

    motor = FlagellaActivity(initial_params)
    settings = motor.default_settings()
    state = settings['state']
    receptor_activity = 1./3.
    state['internal']['chemoreceptor_activity'] = receptor_activity

    flagella_vec = []
    CheY_P_vec = []
    time_vec = []

    # run simulation
    time = 0
    timestep = 0.01  # sec
    while time < total_time:
        time += timestep

        update = motor.next_update(timestep, state)
        CheY_P = update['internal']['CheY_P']
        motor_states = update['flagella']

        flagella_vec.append(motor_states)
        CheY_P_vec.append(CheY_P)
        time_vec.append(time)


    return {
        'flagella_vec': flagella_vec,
        'CheY_P_vec': CheY_P_vec,
        'time_vec': time_vec}

def plot_motor_control(output, out_dir='out'):
    # TODO -- make this into an analysis figure
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib import colors

    expected_run = 0.42  # s (Berg) expected run length without chemotaxis
    expected_tumble = 0.14  # s (Berg)


    # # receptor_activities = output['receptor_activities']
    CheY_P_vec = output['CheY_P_vec']
    flagella_vec = output['flagella_vec']
    time_vec = output['time_vec']

    # get all flagella states
    initial_flagella = flagella_vec[0].keys()
    activity_grid = np.zeros((len(initial_flagella), len(time_vec)))
    for time_index, flagella_states in enumerate(flagella_vec):
        for flagella_id, motor_state in flagella_states.items():
            flagella_index = initial_flagella.index(flagella_id)
            activity_grid[flagella_index, time_index] = motor_state + 1

    # plot results
    cols = 1
    rows = 2
    plt.figure(figsize=(10 * cols, 2 * rows))


    ax1 = plt.subplot(rows, cols, 1)

    ## plot all flagella states in grid
    # make a color map of fixed colors
    cmap = colors.ListedColormap(['white', 'red', 'blue'])
    bounds = [0, 0.5, 1.5, 2]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # plot all flagella states in grid
    im = ax1.imshow(activity_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap,
               norm=norm
               )
    cbar = plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,1,2])
    cbar.set_ticklabels(['none', 'run', 'tumble'])
    plt.locator_params(axis='y', nbins=len(initial_flagella))
    ax1.set_yticks([])
    ax1.set_ylabel('flagella #')
    ax1.set_xlabel('time')


    # save the figure
    fig_path = os.path.join(out_dir, 'motor_control')
    plt.subplots_adjust(wspace=0.7, hspace=0.5)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'flagella_activity')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output1 = test_motor_control(200)
    plot_motor_control(output1, out_dir)
