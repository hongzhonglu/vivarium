from __future__ import absolute_import, division, print_function

import numpy as np
import random

from lens.actor.process import Process


TUMBLE_JITTER = 0.5  # (radians)

## Parameters
k_A = 5.0  #
k_y = 100.0  # 1/uM/s
k_z = 30.0  # / self.CheZ,
gamma_Y = 0.1
k_s = 0.45  # scaling coefficient

receptor_adapt_rate = 1000
adaptPrecision = 0.3  # set to 1.0 for perfect adaptation

# initial concentrations
CheR_WT = 0.00016  # (mM) wild type concentration. 0.16 uM = 0.00016 mM
CheB_WT = 0.00028  # (mM) wild type concentration. 0.28 uM = 0.00028 mM

# motor
mb_0 = 0.65
n_motors = 5
cw_to_ccw = 0.83  # 1/s (Block1983) motor bias, assumed to be constant


class FlagellaActivity(Process):
    def __init__(self, initial_parameters={}):

        roles = {
            'internal': [],
            'external': []
        }
        parameters = {}
        parameters.update(initial_parameters)

        super(FlagellaActivity, self).__init__(roles, parameters, deriver=True)

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''

        return {
            'external': {},
            'internal': {}}

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

        # update receptor clusters and get their probability of being active
        # for cluster in self.receptor_clusters:
        #     cluster.update(ligand_conc, self.CheR, self.CheB)
        clusters_p_on = [cluster.P_on for cluster in self.receptor_clusters]
        P_on = sum(clusters_p_on) / self.n_receptors_clusters

        # CheY phosphorylation, scaled to 1.0 at rest
        # self.CheA_P = self.CheA_tot * self.P_on * k_A / (self.P_on * k_A + k_y * self.CheY_tot)  # amount of phosphorylated CheA
        # scaling = 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        # self.CheY_P = adaptPrecision * scaling * self.CheY_tot * k_y * self.CheA_P / (k_y * self.CheA_P + k_z * self.CheZ + gamma_Y)
        self.CheY_P = adaptPrecision * 3 * k_y * self.CheY_tot * P_on * k_s / (k_y * k_s * P_on + k_z * self.CheZ + gamma_Y)

        # CheB phosphorylation
        # TODO!



        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = mb_0 / (self.CheY_P/self.CheY_tot * (1 - mb_0) + mb_0)
        cww_to_cw = cw_to_ccw * (1 / ccw_motor_bias - 1)

        # print("current state: " + str(self.motor_state))
        # print("ccw_motor_bias: " + str(ccw_motor_bias))
        # print("cw_to_ccw: " + str(cw_to_ccw))
        # print("cww_to_cw: " + str(cww_to_cw))

        if self.motor_state is 'run':
            # switch to tumble?
            prob_switch = cww_to_cw * self.timestep
            print('prob_switch ' + str(prob_switch))
            if prob_switch <= np.random.random(1)[0]:
                self.motor_state = 'tumble'
                self.tumble()
            else:
                self.run()

        elif self.motor_state is 'tumble':
            # switch to run?
            prob_switch = cw_to_ccw * self.timestep
            print('prob_switch ' + str(prob_switch))
            if prob_switch <= np.random.random(1)[0]:
                self.motor_state = 'run'
                self.run()
            else:
                self.tumble()


        update = {}

        return update

    def tumble(self):
        force = 5.0
        torque = random.normalvariate(0, TUMBLE_JITTER)
        self.motile_force = [force, torque]
        # print('TUMBLE!')

    def run(self):
        force = 15.0
        torque = 0.0
        self.motile_force = [force, torque]
        # print('RUN!')