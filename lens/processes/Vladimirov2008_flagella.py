from __future__ import absolute_import, division, print_function

import numpy as np
import random

from lens.actor.process import Process


TUMBLE_JITTER = 0.5  # (radians)

## Parameters
DEFAULT_PARAMETERS = {
    'k_A': 5.0,  #
    'k_y': 100.0,  # 1/uM/s
    'k_z': 30.0,  # / self.CheZ,
    'gamma_Y': 0.1,
    'k_s': 0.45,  # scaling coefficient

    'receptor_adapt_rate': 1000,
    'adaptPrecision': 0.3,  # set to 1.0 for perfect adaptation

    # motor
    'mb_0': 0.65,
    'n_motors': 5,
    'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
}

# initial concentrations
initial = {
    'CheR_WT': 0.00016,  # (mM) wild type concentration. 0.16 uM = 0.00016 mM
    'CheB_WT': 0.00028,  # (mM) wild type concentration. 0.28 uM = 0.00028 mM
}

class FlagellaActivity(Process):
    def __init__(self, initial_parameters={}):

        roles = {
            'internal': [],
            'external': []
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(roles, parameters, deriver=True)

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''

        internal = {'chemoreceptor_P_on': 0.5,
                    'motor_state': 'tumble'}

        return {
            'external': {},
            'internal': internal}

    def default_emitter_keys(self):
        keys = {
            'internal': [],
            'external': [],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''

        updater_types = {
            'internal': {},
            'external': {}}

        return updater_types

    def next_update(self, timestep, states):

        internal = states['internal']
        external = states['external']
        P_on = internal['chemoreceptor_P_on']
        motor_state = internal['motor_state']
        CheZ = internal['CheZ']
        CheY_tot = internal['CheY_tot']

        # CheY phosphorylation, scaled to 1.0 at rest
        # self.CheA_P = self.CheA_tot * self.P_on * k_A / (self.P_on * k_A + k_y * self.CheY_tot)  # amount of phosphorylated CheA
        # scaling = 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        # self.CheY_P = adaptPrecision * scaling * self.CheY_tot * k_y * self.CheA_P / (k_y * self.CheA_P + k_z * self.CheZ + gamma_Y)
        CheY_P = self.parameters['adaptPrecision'] * 3 * self.parameters['k_y'] * CheY_tot * P_on * self.parameters['k_s'] / \
                 (self.parameters['k_y'] * self.parameters['k_s'] * P_on + self.parameters['k_z'] * CheZ + self.parameters['gamma_Y'])

        # CheB phosphorylation
        # TODO!


        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = self.parameters['mb_0'] / (CheY_P/CheY_tot * (1 - self.parameters['mb_0']) + self.parameters['mb_0'])
        cww_to_cw = self.parameters['cw_to_ccw'] * (1 / ccw_motor_bias - 1)


        if motor_state is 'run':
            # switch to tumble?
            prob_switch = cww_to_cw * timestep
            print('prob_switch ' + str(prob_switch))
            if prob_switch <= np.random.random(1)[0]:
                motor_state = 'tumble'
                motile_force = self.tumble()
            else:
                motile_force = self.run()

        elif motor_state is 'tumble':
            # switch to run?
            prob_switch = self.parameters['cw_to_ccw'] * timestep
            print('prob_switch ' + str(prob_switch))
            if prob_switch <= np.random.random(1)[0]:
                motor_state = 'run'
                motile_force = self.run()
            else:
                motile_force = self.tumble()


        update = {
            'motile_force': motile_force,
            'motor_state': motor_state,
            'CheY_P': CheY_P}

        return update

    def tumble(self):
        force = 5.0
        torque = random.normalvariate(0, TUMBLE_JITTER)
        return [force, torque]

    def run(self):
        force = 15.0
        torque = 0.0
        return [force, torque]
