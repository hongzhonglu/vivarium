from __future__ import absolute_import, division, print_function

import os
import math

from vivarium.actor.process import Process
from vivarium.utils.units import units


LIGAND_ID = 'MeAsp'

INITIAL_STATE = {
    'n_methyl': 2.0,  # initial number of methyl groups on receptor cluster (0 to 8)
    'chemoreceptor_activity': 1./3.,  # initial probability of receptor cluster being on
    'CheR': 0.00016,  # (mM) wild type concentration. 0.16 uM = 0.00016 mM
    'CheB': 0.00028,  # (mM) wild type concentration. 0.28 uM = 0.00028 mM. [CheR]:[CheB]=0.16:0.28
    'CheB_P': 0.0,  # phosphorylated CheB
}

# Parameters from Endres and Wingreen 2006. See paper for serine binding, NiCl2 binding rates.
DEFAULT_PARAMETERS = {
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
    # k_CheR = 0.0182  # effective catalytic rate of CheR
    # k_CheB = 0.0364  # effective catalytic rate of CheB
    'k_meth': 0.0625,  # Catalytic rate of methylation
    'k_demeth': 0.0714,  # Catalytic rate of demethylation
    'adaptRate': 2,  # adaptation rate relative to wild-type. cell-to-cell variation cause by variability in [CheR, CheB]
}


class ReceptorCluster(Process):
    def __init__(self, initial_parameters={}):

        self.ligand_id = initial_parameters.get('ligand', LIGAND_ID)
        roles = {
            'internal': ['n_methyl', 'chemoreceptor_activity', 'CheR', 'CheB'],
            'external': [self.ligand_id]
        }
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(ReceptorCluster, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        internal = INITIAL_STATE
        internal.update({'volume': 1})
        external = {self.ligand_id: 0.1}
        default_state = {
            'external': external,
            'internal': internal}

        # default emitter keys
        default_emitter_keys = {
            'internal': ['chemoreceptor_activity', 'n_methyl'],
            'external': [self.ligand_id]}

        # default updaters
        default_updaters = {
            'internal': {state_id: 'set' for state_id in ['chemoreceptor_activity', 'n_methyl']},
            'external': {}}

        default_settings = {
            'process_id': 'receptor',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'time_step': 1.0}

        return default_settings

    def next_update(self, timestep, states):
        '''
        Monod-Wyman-Changeux model for mixed cluster activity from:
            Endres & Wingreen. (2006). Precise adaptation in bacterial chemotaxis through "assistance neighborhoods"
        '''
        # states
        n_methyl = states['internal']['n_methyl']
        P_on = states['internal']['chemoreceptor_activity']
        CheR = states['internal']['CheR'] * (units.mmol / units.L)
        CheB = states['internal']['CheB'] * (units.mmol / units.L)
        ligand_conc = states['external'][self.ligand_id]   # mmol/L

        # convert to umol / L
        CheR = CheR.to('umol/L').magnitude
        CheB = CheB.to('umol/L').magnitude

        # parameters
        n_Tar = self.parameters['n_Tar']
        n_Tsr = self.parameters['n_Tsr']
        K_Tar_off = self.parameters['K_Tar_off']
        K_Tar_on = self.parameters['K_Tar_on']
        K_Tsr_off = self.parameters['K_Tsr_off']
        K_Tsr_on = self.parameters['K_Tsr_on']
        adaptRate = self.parameters['adaptRate']
        k_meth = self.parameters['k_meth']
        k_demeth = self.parameters['k_demeth']

        if n_methyl < 0:
            n_methyl = 0
        elif n_methyl > 8:
            n_methyl = 8
        else:
            d_methyl = adaptRate * (k_meth * CheR * (1.0 - P_on) - k_demeth * CheB * P_on) * timestep
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

        # free energy of the receptors.
        # TODO -- generalize to other ligands (repellents). See Endres & Wingreen 2006 eqn 9 in supporting info
        Tar_free_energy = n_Tar * (offset_energy + math.log((1+ligand_conc/K_Tar_off) / (1+ligand_conc/K_Tar_on)))
        Tsr_free_energy = n_Tsr * (offset_energy + math.log((1+ligand_conc/K_Tsr_off) / (1+ligand_conc/K_Tsr_on)))

        # free energy of receptor clusters
        cluster_free_energy = Tar_free_energy + Tsr_free_energy  #  free energy of the cluster
        P_on = 1.0/(1.0 + math.exp(cluster_free_energy))  # probability that receptor cluster is ON (CheA is phosphorylated). Higher free energy --> activity less probably

        update = {
            'internal': {
                'chemoreceptor_activity': P_on,
                'n_methyl': n_methyl,
            }}
        return update


# tests and analyses of process
def test_receptor():
    # TODO -- add asserts for test
    # define timeline with (time (s), ligand concentration (mmol/L))
    timeline = [
        (0, 0.0),
        (200, 0.01),
        (600, 0.0),
        (1000, 0.1),
        (1400, 0.0),
        (1800, 1.0),
        (2200, 0.0),
        (2600, 10.0),
        (3000, 0.0),
    ]
    time = 0
    timestep = 1

    # initialize process
    receptor = ReceptorCluster()
    settings = receptor.default_settings()
    state = settings['state']
    ligand_id = receptor.ligand_id

    # run simulation
    ligand_vec = []
    receptor_activity_vec = []
    n_methyl_vec = []
    while time < timeline[-1][0]:
        time += timestep
        for (t, conc) in timeline:
            if time < t:
                pass
            else:
                ligand_conc = conc
                state['external'][ligand_id] = ligand_conc
                continue

        # run step
        update = receptor.next_update(timestep, state)
        P_on = update['internal']['chemoreceptor_activity']
        n_methyl = update['internal']['n_methyl']

        # update state
        state['internal']['chemoreceptor_activity'] = P_on
        state['internal']['n_methyl'] = n_methyl

        # save state
        ligand_vec.append(ligand_conc)
        receptor_activity_vec.append(P_on)
        n_methyl_vec.append(n_methyl)

    return {
        'ligand_vec': ligand_vec,
        'receptor_activity_vec': receptor_activity_vec,
        'n_methyl_vec': n_methyl_vec}

def plot_output(output, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    ligand_vec = output['ligand_vec']
    receptor_activity_vec = output['receptor_activity_vec']
    n_methyl_vec = output['n_methyl_vec']

    # plot results
    cols = 1
    rows = 3
    plt.figure(figsize=(6, 6))

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)

    ax1.plot(ligand_vec, 'b')
    ax2.plot(receptor_activity_vec, 'b')
    ax3.plot(n_methyl_vec, 'b')

    ax1.set_xticklabels([])
    ax1.set_ylabel("ligand \n log(mM) ", fontsize=10)
    ax1.set_yscale('log')
    ax2.set_xticklabels([])
    ax2.set_ylabel("activity \n P(on)", fontsize=10)
    ax3.set_xlabel("time (s)", fontsize=12)
    ax3.set_ylabel("average \n methylation", fontsize=10)

    fig_path = os.path.join(out_dir, 'response')
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    output = test_receptor()
    out_dir = os.path.join('out', 'tests', 'Endres2006_chemoreceptor')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_output(output, out_dir)
