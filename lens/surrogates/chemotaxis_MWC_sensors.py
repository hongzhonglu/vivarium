
from __future__ import absolute_import, division, print_function

import time
import numpy as np
import math
import random

from lens.actor.inner import Simulation


TUMBLE_JITTER = 0.5  # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 0 , 127]]

# MeAsp is an attractant
# TODO (Eran) -- add NiCl2, a repellent
LIGAND = 'GLC'  # in the original model the ligand is 'MeAsp'

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


class Chemotaxis(Simulation):
    '''
    Simple chemotaxis surrogate that can move up glucose gradients. It can take on two states 'run' and 'tumble'.
    State is a function of the current glucose concentrations, and internal CheY concentrations -- the 'memory' of
    glucose concentrations from the previous time step.
    Based on the model reported in:
    - Vladimirov, Nikita. Multiscale modeling of bacterial chemotaxis. Diss. 2009.
    - http://www.rapidcell.net
    '''

    def __init__(self, state):
        self.initial_time = state.get('time', 0.0)
        self.local_time = 0.0
        self.timestep = 1.0
        self.environment_change = {}
        self.volume = 1.0
        self.division_time = 100
        self.color = DEFAULT_COLOR

        # initial state
        self.motor_state = 'tumble'
        self.external = {'GLC': 0.0}

        ## Initial state
        # receptor-activated kinase
        # self.CheA_tot = 1000.0
        # self.CheA_P = 0.0

        # methylation regulators of receptors
        self.CheR = CheR_WT
        self.CheB = CheB_WT
        self.CheB_P = 0.0	 # (mM)

        # response regulator proteins
        self.CheY_tot = 0.0097  # (mM) 9.7 uM = 0.0097 mM
        self.CheY_P = 0.0
        self.CheZ = 0.0		# phosphatase

        # make cluster objects
        self.n_receptors_clusters = 1
        self.receptor_clusters = []
        for n in xrange(self.n_receptors_clusters):
            self.receptor_clusters.append(ReceptorCluster(receptor_adapt_rate, self.timestep))

        # other states
        self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
        self.division = []


    def update_state(self):

        ligand_conc = self.external[LIGAND]  # (mM)

        # update receptor clusters and get their probability of being active
        for cluster in self.receptor_clusters:
            cluster.update(ligand_conc, self.CheR, self.CheB)
        clusters_p_on = [cluster.P_on for cluster in self.receptor_clusters]
        P_on = sum(clusters_p_on) / self.n_receptors_clusters

        # CheY phosphorylation, scaled to 1.0 at rest
        # self.CheA_P = self.CheA_tot * self.P_on * k_A / (self.P_on * k_A + k_y * self.CheY_tot)  # amount of phosphorylated CheA
        # scaling = 19.3610  # scales CheY_P linearly so that CheY_P=1 at rest (P_on=1/3)
        # self.CheY_P = adaptPrecision * scaling * self.CheY_tot * k_y * self.CheA_P / (k_y * self.CheA_P + k_z * self.CheZ + gamma_Y)
        self.CheY_P = adaptPrecision * 3 * k_y * self.CheY_tot * P_on * k_s / (k_y * k_s * P_on + k_z * self.CheZ + gamma_Y)

        # CheB phosphorylation
        # TODO!

        # print('P_on: ' + str(P_on))
        # print('CheY_P: ' + str(self.CheY_P))
        # print('CheY_tot: ' + str(self.CheY_tot))

    def motor_switching(self):

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

    def check_division(self):
        # update division state based on time since initialization
        if self.local_time >= self.initial_time + self.division_time:
            self.division = [{'time': self.local_time}, {'time': self.local_time}]

        return self.division

    def time(self):
        return self.local_time

    def apply_outer_update(self, update):
        self.external = update['concentrations']

        self.environment_change = {}
        for molecule in self.external.iterkeys():
            self.environment_change[molecule] = 0

    def run_incremental(self, run_until):
        # update state once per message exchange
        while self.local_time < run_until:
            self.update_state()
            self.motor_switching()
            # self.check_division()
            self.local_time += self.timestep

        time.sleep(0.2)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            'motile_force': self.motile_force,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': self.color,
            }

    def synchronize_state(self, state):
        if 'time' in state:
            self.initial_time = state['time']


class ReceptorCluster(object):
    def __init__(self, adaptRate, timestep):

        self.timestep = timestep

        # initial number of methyl groups on receptor cluster (0 to 8)
        self.n_methyl = 3.0

        # initial probability of receptor cluster being on
        self.P_on = 1/3

        ## Parameters
        self.n_Tar = 6  # number of Tar receptors in a cluster
        self.n_Tsr = 12  # number of Tsr receptors in a cluster

        # dissociation constants (mM)
        # K_Tar_on = 12e-3  # Tar to Asp (Emonet05)
        # K_Tar_off = 1.7e-3  # Tar to Asp (Emonet05)
        # (Endres & Wingreen, 2006) has dissociation constants for serine binding, NiCl2 binding
        self.K_Tar_off = 0.02    # (mM) MeAsp binding by Tar (Endres06)
        self.K_Tar_on = 0.5    # (mM) MeAsp binding by Tar (Endres06)
        self.K_Tsr_off = 100.0  # (mM) MeAsp binding by Tsr (Endres06)
        self.K_Tsr_on = 10e6    # (mM) MeAsp binding by Tsr (Endres06)

        # self.k_CheR = 0.0182  # effective catalytic rate of CheR
        # self.k_CheB = 0.0364  # effective catalytic rate of CheB

        self.k_meth = 0.0625    # Catalytic rate of methylation
        self.k_demeth = 0.0714  # Catalytic rate of demethylation

        self.adaptRate = adaptRate

    def update(self, ligand_conc, CheR, CheB):
        '''
        Monod-Wyman-Changeux model for mixed cluster activity
        from Endres & Wingreen. (2006). Precise adaptation in bacterial chemotaxis through "assistance neighborhoods"
        '''
        if self.n_methyl < 0:
            self.n_methyl = 0
        elif self.n_methyl > 8:
            self.n_methyl = 8
        else:
            # d_methyl = self.adaptRate * (self.k_CheR * CheR * (1.0 - self.P_on) - self.k_CheB * CheB * self.P_on) * self.timestep
            d_methyl = self.adaptRate * (self.k_meth * CheR * (1.0 - self.P_on) - self.k_demeth * CheB * self.P_on) * self.timestep
            self.n_methyl += d_methyl

        # get free-energy offsets from methylation
        # piece-wise linear model. Assumes same offset energy (epsilon) for both Tar and Tsr
        if self.n_methyl < 0:
            offset_energy = 1.0
        elif self.n_methyl < 2:
            offset_energy = 1.0 - 0.5 * self.n_methyl
        elif self.n_methyl < 4:
            offset_energy = -0.3 * (self.n_methyl - 2.0)
        elif self.n_methyl < 6:
            offset_energy = -0.6 - 0.25 * (self.n_methyl - 4.0)
        elif self.n_methyl < 7:
            offset_energy = -1.1 - 0.9 * (self.n_methyl - 6.0)
        elif self.n_methyl < 8:
            offset_energy = -2.0 - (self.n_methyl - 7.0)
        else:
            offset_energy = -3.0

        # free energy of the receptors
        # TODO -- generalize to multiple types of ligands. See eqn 9 in supporting info for Endres & Wingreen 2006
        Tar_free_energy = self.n_Tar * (offset_energy + math.log((1+ligand_conc/self.K_Tar_off)/(1+ligand_conc/self.K_Tar_on)))
        Tsr_free_energy = self.n_Tsr * (offset_energy + math.log((1+ligand_conc/self.K_Tsr_off)/(1+ligand_conc/self.K_Tsr_on)))

        # free energy of receptor clusters
        cluster_free_energy = Tar_free_energy + Tsr_free_energy  #  free energy of the cluster
        self.P_on = 1.0/(1.0 + math.exp(cluster_free_energy))  # probability that receptor cluster is ON (CheA is phosphorylated)
