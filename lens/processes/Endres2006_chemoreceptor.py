from __future__ import absolute_import, division, print_function

import math

from lens.actor.process import Process


LIGAND_ID = 'MeAsp'

INITIAL_STATE = {
    'n_methyl': 3.0,  # initial number of methyl groups on receptor cluster (0 to 8)
    'chemoreceptor_P_on': 1 / 3,  # initial probability of receptor cluster being on
}

DEFAULT_PARAMETERS = {
    'timestep': 1.0,
    'n_Tar': 6,  # number of Tar receptors in a cluster
    'n_Tsr': 12,  # number of Tsr receptors in a cluster
    # dissociation constants (mM)
    # K_Tar_on = 12e-3  # Tar to Asp (Emonet05)
    # K_Tar_off = 1.7e-3  # Tar to Asp (Emonet05)
    # (Endres & Wingreen, 2006) has dissociation constants for serine binding, NiCl2 binding
    'K_Tar_off': 0.02,  # (mM) MeAsp binding by Tar (Endres06)
    'K_Tar_on': 0.5,  # (mM) MeAsp binding by Tar (Endres06)
    'K_Tsr_off': 100.0,  # (mM) MeAsp binding by Tsr (Endres06)
    'K_Tsr_on': 10e6,  # (mM) MeAsp binding by Tsr (Endres06)
    # self.k_CheR = 0.0182  # effective catalytic rate of CheR
    # self.k_CheB = 0.0364  # effective catalytic rate of CheB
    'k_meth': 0.0625,  # Catalytic rate of methylation
    'k_demeth': 0.0714,  # Catalytic rate of demethylation
    'adaptRate': 1000,
}



class ReceptorCluster(Process):
    def __init__(self, initial_parameters={}):

        self.ligand_id = LIGAND_ID

        roles = {
            'internal': ['n_methyl', 'chemoreceptor_P_on', 'CheR', 'CheB'],
            'external': [self.ligand_id]
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(ReceptorCluster, self).__init__(roles, parameters, deriver=True)

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''

        internal = INITIAL_STATE
        internal.update({'volume': 1})
        external = {self.ligand_id: 0.0}

        return {
            'external': external,
            'internal': internal}

    def default_emitter_keys(self):
        keys = {
            'internal': ['chemoreceptor_P_on', 'n_methyl'],
            'external': [self.ligand_id],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''

        updater_types = {
            'internal': {state_id: 'set' for state_id in ['chemoreceptor_P_on', 'n_methyl']},
            'external': {}}

        return updater_types

    def next_update(self, timestep, states):
        '''
        Monod-Wyman-Changeux model for mixed cluster activity
        from Endres & Wingreen. (2006). Precise adaptation in bacterial chemotaxis through "assistance neighborhoods"
        '''

        n_methyl = states['internal']['n_methyl']
        P_on = states['internal']['chemoreceptor_P_on']
        CheR = states['internal']['CheR']
        CheB = states['internal']['CheB']

        ligand_conc = states['external'][self.ligand_id]

        if n_methyl < 0:
            n_methyl = 0
        elif n_methyl > 8:
            n_methyl = 8
        else:
            d_methyl = self.parameters['adaptRate'] * (self.parameters['k_meth'] * CheR * (1.0 - P_on) - self.parameters['k_demeth'] * CheB * P_on) * self.parameters['timestep']
            n_methyl += d_methyl

        # get free-energy offsets from methylation
        # piece-wise linear model. Assumes same offset energy (epsilon) for both Tar and Tsr
        if n_methyl < 0:
            offset_energy = 1.0
        elif n_methyl < 2:
            offset_energy = 1.0 - 0.5 * n_methyl
        elif n_methyl < 4:
            offset_energy = -0.3 * (n_methyl - 2.0)
        elif n_methyl < 6:
            offset_energy = -0.6 - 0.25 * (n_methyl - 4.0)
        elif n_methyl < 7:
            offset_energy = -1.1 - 0.9 * (n_methyl - 6.0)
        elif n_methyl < 8:
            offset_energy = -2.0 - (n_methyl - 7.0)
        else:
            offset_energy = -3.0

        # free energy of the receptors
        # TODO -- generalize to multiple types of ligands. See eqn 9 in supporting info for Endres & Wingreen 2006
        Tar_free_energy = self.parameters['n_Tar'] * \
                          (offset_energy + math.log((1+ligand_conc/self.parameters['K_Tar_off'])/
                          (1+ligand_conc/self.parameters['K_Tar_on'])))
        Tsr_free_energy = self.parameters['n_Tsr'] * \
                          (offset_energy + math.log((1+ligand_conc/self.parameters['K_Tsr_off'])/
                          (1+ligand_conc/self.parameters['K_Tsr_on'])))

        # free energy of receptor clusters
        cluster_free_energy = Tar_free_energy + Tsr_free_energy  #  free energy of the cluster
        P_on = 1.0/(1.0 + math.exp(cluster_free_energy))  # probability that receptor cluster is ON (CheA is phosphorylated)

        update = {
            'internal': {
                'chemoreceptor_P_on': P_on,
                'n_methyl': n_methyl,
            }
        }

        return update
