from __future__ import absolute_import, division, print_function

import os
import numpy as np
import random
import uuid

from lens.actor.process import Process, deep_merge


# parameters, from Mears et al.
DEFAULT_N_FLAGELLA = 3
DEFAULT_PMF = 180
DEFAULT_PARAMETERS = {
    'ccw_to_cw': 0.26,  # (1/s) Motor switching rate from CCW->CW
    'cw_to_ccw': 1.7,  # (1/s) Motor switching rate from CW->CCW
    'CB': 0.13,  # average CW bias of wild-type motors
    'omega': 0.5,  # (s) characteristic motor switch time
    'lambda': 0.68,  # (1/s) transition rate from semi-coiled to curly-w state
    # 'x': DEFAULT_N_FLAGELLA,  # number of flagella that must be normal for a run to occur

    # CheY-P flucutations
    'YP_ss': 2.59,  # (uM) steady state concentration of CheY-P
    'sigma2_Y': 1.0,  # (uM^2) variance in CheY-P
    'tau': 0.2,  # (s) characteristic time-scale fluctuations in [CheY-P]

    # CW bias
    'K_d': 3.1,  # (uM) midpoint of CW bias vs CheY-P response curve
    'H': 10.3,  # Hill coefficient for CW bias vs CheY-P response curve
}

##initial state
INITIAL_STATE = {
    'CheY': 2.59,
    'CheY_P': 2.59,  # (uM) mean concentration of CheY-P
    'cw_bias': 0.5,  # (made up)
}

class FlagellaActivity(Process):
    '''
    Model of multi-flagellar motor activity and CheY-P fluctuations from:
        Mears, P. J., Koirala, S., Rao, C. V., Golding, I., & Chemla, Y. R. (2014).
        Escherichia coli swimming is robust against variations in flagellar number.
    '''
    def __init__(self, initial_parameters={}):

        self.n_flagella = initial_parameters.get('n_flagella', DEFAULT_N_FLAGELLA)
        self.flagella_ids = [str(uuid.uuid1()) for flagella in range(self.n_flagella)]

        roles = {
            'internal': [
                'chemoreceptor_activity',
                'CheY',
                'CheY_P',
                'cw_bias'
            ],
            'membrane': ['PMF', 'protons_flux_accumulated'],
            'flagella': self.flagella_ids,
            'external': []
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # flagella motor state: 0 for CCW, 1 for CW
        # PMF ranges 180-200 mV (Berg)
        internal = INITIAL_STATE
        default_state = {
            'external': {},
            'membrane': {'PMF': DEFAULT_PMF, 'PROTONS': 0},
            'flagella': {flagella_id: random.choice([0, 1]) for flagella_id in self.flagella_ids},
            'internal': deep_merge(internal, {'volume': 1, 'n_flagella': DEFAULT_N_FLAGELLA})}

        # default emitter keys
        default_emitter_keys = {
            'internal': [],
            'flagella': [],
            'external': [],
        }

        # default updaters
        internal_set_states = ['cw_bias']
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
        YP = internal['CheY_P']

        # parameters
        tau = self.parameters['tau']
        YP_ss = self.parameters['YP_ss']
        sigma = self.parameters['sigma2_Y']**0.5

        K_d = self.parameters['K_d']
        H = self.parameters['H']

        ## update CheY-P
        dYP = -(1 / tau) * (YP - YP_ss) * timestep + sigma * (2 * timestep / tau)**0.5 * random.normalvariate(0, 1)

        ## CW bias
        # Hill function from Cluzel, P., Surette, M., & Leibler, S. (2000).
        # An ultrasensitive bacterial motor revealed by monitoring signaling proteins in single cells.
        cw_bias = YP ** H / (K_d ** H + YP ** H)

        ## update flagella
        # determine behavior from motor states of all flagella
        flagella_update = {}
        for flagella_id, motor_state in flagella.items():
            new_motor_state = self.update_flagellum(motor_state, cw_bias, timestep)
            flagella_update.update({flagella_id: new_motor_state})

        return {
            'flagella': flagella_update,
            'internal' : {
                'CheY_P': dYP,
                'cw_bias': cw_bias,
            }}

    def update_flagellum(self, motor_state, cw_bias, timestep):
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

        ccw_to_cw = self.parameters['ccw_to_cw']
        cw_to_ccw = self.parameters['cw_to_ccw']

        # import ipdb; ipdb.set_trace()

        # # CCW corresponds to run. CW corresponds to tumble
        # ccw_motor_bias = mb_0 / (CheY_P * (1 - mb_0) + mb_0)  # (1/s)
        # ccw_to_cw = cw_to_ccw * (1 / ccw_motor_bias - 1)  # (1/s)
        if motor_state == 0:  # 0 for run
            # switch to tumble?
            prob_switch = ccw_to_cw * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 1
            else:
                new_motor_state = 0

        elif motor_state == 1:  # 1 for tumble
            # switch to run?
            prob_switch = cw_to_ccw * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 0
            else:
                new_motor_state = 1

        return new_motor_state



# testing functions
def test_motor_control(total_time=10):
    # TODO -- add asserts for test

    initial_params = {}

    motor = FlagellaActivity(initial_params)
    settings = motor.default_settings()
    state = settings['state']

    receptor_activity = 1./3.
    state['internal']['chemoreceptor_activity'] = receptor_activity

    saved_data = {
        'internal': {state_id: [value] for state_id, value in state['internal'].items()},
        'flagella': {state_id: [value] for state_id, value in state['flagella'].items()},
        'time': [0]}


    accumulate_internal_states = ['CheY_P']

    # run simulation
    time = 0
    timestep = 0.001  # sec
    while time < total_time:
        time += timestep

        update = motor.next_update(timestep, state)

        # flagella state is set
        state['flagella'] = update['flagella']

        # apply updates to internal state
        for state_id, value in update['internal'].items():
            if state_id in accumulate_internal_states:
                # accumulate
                state['internal'][state_id] += value

                if state['internal'][state_id] < 0:  # TODO -- why does CheY-P go below 0?
                    state['internal'][state_id] = 0
            else:
                # set
                state['internal'][state_id] = value

        saved_data['time'].append(time)
        for role in ['internal', 'flagella',]:
            for state_id, value in state[role].items():
                saved_data[role][state_id].append(value)

    return saved_data

def plot_motor_control(output, out_dir='out'):
    # TODO -- make this into an analysis figure
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib import colors

    # receptor_activities = output['receptor_activities']
    CheY_P_vec = output['internal']['CheY_P']
    cw_bias_vec = output['internal']['cw_bias']
    flagella = output['flagella']
    time_vec = output['time']

    # get all flagella states into an activity grid
    flagella_ids = flagella.keys()
    activity_grid = np.zeros((len(flagella_ids), len(time_vec)))
    total_CW = np.zeros((len(time_vec)))
    for flagella_id, motor_states in flagella.items():
        flagella_index = flagella_ids.index(flagella_id)
        activity_grid[flagella_index, :] = [x + 1 for x in motor_states]

        total_CW += np.array(motor_states)

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
    ax1.plot(time_vec, CheY_P_vec)
    ax1.set_ylabel('[CheY-P] (uM)')

    # plot CW bias
    ax2.plot(time_vec, cw_bias_vec)
    ax2.set_ylabel('CW bias')

    # plot cell state
    ax3.set_ylabel('cell')

    # plot flagella states in a grid
    cmap = colors.ListedColormap(['black', 'white', 'blue'])
    bounds = [0, 0.5, 1.5, 2]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    im = ax4.imshow(activity_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap,
               norm=norm)
    # cbar = plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,1,2])
    # cbar.set_ticklabels(['none', 'CCW', 'CW'])
    plt.locator_params(axis='y', nbins=len(flagella_ids))
    ax4.set_yticks(list(range(len(flagella_ids))))
    ax4.set_ylabel('flagella #')
    ax4.set_xlabel('time')


    # plot number of flagella CW
    ax5.plot(time_vec, total_CW)
    ax5.set_ylabel('number of flagella CW')


    # save figure
    fig_path = os.path.join(out_dir, 'motor_control')
    plt.subplots_adjust(wspace=0.7, hspace=0.5)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'Mears2014_flagella_activity')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output1 = test_motor_control(20)
    plot_motor_control(output1, out_dir)
