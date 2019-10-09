from __future__ import absolute_import, division, print_function

import numpy as np
import random

from lens.actor.process import Process, dict_merge


TUMBLE_JITTER = 0.5  # (radians)

# parameters
DEFAULT_PARAMETERS = {
    'k_A': 5.0,  #
    'k_y': 100.0,  # 1/uM/s
    'k_z': 30.0,  # / self.CheZ,
    'gamma_Y': 0.1,
    'k_s': 0.45,  # scaling coefficient
    'adaptPrecision': 0.3,  # set to 1.0 for perfect adaptation
    # motor
    'mb_0': 0.65,
    'n_motors': 5,
    'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
}

##initial state
INITIAL_STATE = {
    # response regulator proteins
    'CheY_tot': 0.0097,  # (mM) 9.7 uM = 0.0097 mM
    'CheY_P': 0.0,
    'CheZ': 0.0,  # phosphatase
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
        PLoS computational biology.
    '''
    def __init__(self, initial_parameters={}):

        roles = {
            'internal': ['chemoreceptor_activity',
                         'CheZ',
                         'CheY_tot',
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
            'internal': ['motile_force', 'motile_torque', 'motor_state'],
            'external': [],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''
        set_states = ['motile_force', 'motile_torque', 'motor_state']
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
        '''

        internal = states['internal']
        # external = states['external']
        P_on = internal['chemoreceptor_activity'] # this comes from chemoreceptor
        motor_state = internal['motor_state']
        CheZ = internal['CheZ']
        CheY_tot = internal['CheY_tot']

        # print('receptor activity: {}'.format(P_on))


        ## CheY phosphorylation, scaled to 1.0 at rest
        # self.CheA_P = self.CheA_tot * self.P_on * k_A / (self.P_on * k_A + k_y * self.CheY_tot)  # amount of phosphorylated CheA
        # scaling = 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        # self.CheY_P = adaptPrecision * scaling * self.CheY_tot * k_y * self.CheA_P / (k_y * self.CheA_P + k_z * self.CheZ + gamma_Y)
        CheY_P = self.parameters['adaptPrecision'] * 3 * self.parameters['k_y'] * CheY_tot * P_on * self.parameters['k_s'] / \
                 (self.parameters['k_y'] * self.parameters['k_s'] * P_on + self.parameters['k_z'] * CheZ + self.parameters['gamma_Y'])

        # CheB phosphorylation
        # TODO!
        # CheB =


        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = self.parameters['mb_0'] / (CheY_P/CheY_tot * (1 - self.parameters['mb_0']) + self.parameters['mb_0'])
        cww_to_cw = self.parameters['cw_to_ccw'] * (1 / ccw_motor_bias - 1)

        if motor_state == 0:  # 0 for run
            # switch to tumble?
            prob_switch = cww_to_cw * timestep
            if prob_switch <= np.random.random(1)[0]:
                motor_state = 1
                force, torque = self.tumble()
            else:
                force, torque = self.run()

        elif motor_state == 1:  # 1 for tumble
            # switch to run?
            prob_switch = self.parameters['cw_to_ccw'] * timestep
            if prob_switch <= np.random.random(1)[0]:
                motor_state = 0
                force, torque = self.run()
            else:
                force, torque = self.tumble()

        # TODO -- should force/torque be accumulated over exchange timestep?
        update = {
            'internal': {
                'motile_force': force,
                'motile_torque': torque,
                'motor_state': motor_state}
        }

        return update

    def tumble(self):
        # print('motor state: TUMBLE!')
        force = 5.0
        torque = random.normalvariate(0, TUMBLE_JITTER)
        return force, torque

    def run(self):
        # print('motor state: RUN!')
        force = 15.0
        torque = 0.0
        return force, torque
