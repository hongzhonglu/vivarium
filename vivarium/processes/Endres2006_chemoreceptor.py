from __future__ import absolute_import, division, print_function

import os
import math
import copy
import random

import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.utils.units import units


LIGAND_ID = 'MeAsp'
STEADY_STATE_DELTA = 1e-6

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
    'adapt_rate': 4,  # adaptation rate relative to wild-type. cell-to-cell variation cause by variability in [CheR, CheB]
}

def run_step(receptor, state, timestep):
    update = receptor.next_update(timestep, state)
    state['internal']['chemoreceptor_activity'] = update['internal']['chemoreceptor_activity']
    state['internal']['n_methyl'] = update['internal']['n_methyl']

def run_to_steady_state(receptor, state, timestep):
    P_on = state['internal']['chemoreceptor_activity']
    n_methyl = state['internal']['n_methyl']
    delta = 1
    while delta > STEADY_STATE_DELTA:
        run_step(receptor, state, timestep)
        d_P_on = P_on - state['internal']['chemoreceptor_activity']
        d_n_methyl = n_methyl - state['internal']['n_methyl']
        delta = (d_P_on**2 + d_n_methyl**2)**0.5
        P_on = state['internal']['chemoreceptor_activity']
        n_methyl = state['internal']['n_methyl']



class ReceptorCluster(Process):

    defaults = {
        'parameters': DEFAULT_PARAMETERS
    }

    def __init__(self, initial_parameters={}):

        self.ligand_id = initial_parameters.get('ligand', LIGAND_ID)
        self.initial_ligand = initial_parameters.get('initial_ligand', 5.0)
        ports = {
            'internal': ['n_methyl', 'chemoreceptor_activity', 'CheR', 'CheB'],
            'external': [self.ligand_id]}

        parameters = self.defaults['parameters']
        parameters.update(initial_parameters)

        super(ReceptorCluster, self).__init__(ports, parameters)

    def default_settings(self):

        set_keys = ['chemoreceptor_activity', 'n_methyl']

        # get default state by running to steady state
        internal = INITIAL_STATE
        external = {self.ligand_id: self.initial_ligand}
        state = {
            'external': external,
            'internal': internal}
        run_to_steady_state(self, state, 1.0)

        # default emitter keys
        default_emitter_keys = {
            'internal': set_keys,
            'external': [self.ligand_id]}

        # schema
        schema = {
            'internal': {
                state_id : {
                    'updater': 'set',
                    'divide': 'set'}
                for state_id in set_keys}}

        default_settings = {
            'process_id': 'receptor',
            'state': state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'time_step': 1.0}

        return default_settings

    def next_update(self, timestep, states):
        '''
        Monod-Wyman-Changeux model for mixed cluster activity from:
            Endres & Wingreen. (2006). Precise adaptation in bacterial chemotaxis through "assistance neighborhoods"
        '''
        # states
        n_methyl = copy.copy(states['internal']['n_methyl'])
        P_on = copy.copy(states['internal']['chemoreceptor_activity'])
        CheR = copy.copy(states['internal']['CheR'] * (units.mmol / units.L))
        CheB = copy.copy(states['internal']['CheB'] * (units.mmol / units.L))
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
        adapt_rate = self.parameters['adapt_rate']
        k_meth = self.parameters['k_meth']
        k_demeth = self.parameters['k_demeth']

        if n_methyl < 0:
            n_methyl = 0
        elif n_methyl > 8:
            n_methyl = 8
        else:
            d_methyl = adapt_rate * (k_meth * CheR * (1.0 - P_on) - k_demeth * CheB * P_on) * timestep
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
def get_pulse_timeline():
    timeline = [
        (0, 0.0),
        (100, 0.01),
        (200, 0.0),
        (300, 0.1),
        (400, 0.0),
        (500, 1.0),
        (600, 0.0),
        (700, 0.0)]
    return timeline

def get_linear_step_timeline(config):
    time = config.get('time', 100)
    slope = config.get('slope', 2e-3)  # mM/um
    speed = config.get('speed', 14)     # um/s
    conc_0 = config.get('initial_conc', 0)  # mM
    timeline = [(t, conc_0 + slope*t*speed) for t in range(time)]
    return timeline

def get_exponential_step_timeline(config):
    time = config.get('time', 100)
    base = config.get('base', 1+1e-4)  # mM/um
    speed = config.get('speed', 14)     # um/s
    conc_0 = config.get('initial_conc', 0)  # mM
    timeline = [(t, conc_0 + base**(t*speed) - 1) for t in range(time)]
    return timeline

def get_exponential_random_timeline(config):
    # exponential space with random direction changes
    time = config.get('time', 100)
    base = config.get('base', 1+1e-4)  # mM/um
    speed = config.get('speed', 14)     # um/s
    conc_0 = config.get('initial_conc', 0)  # mM

    conc = conc_0
    timeline = [(0, conc)]
    for t in range(time):
        conc += base**(random.choice((-1, 1)) * speed) - 1
        if conc<0:
            conc = 0
        timeline.append((t, conc))

    return timeline

def test_receptor(timeline=get_pulse_timeline(), timestep = 1):

    end_time = timeline[-1][0]
    time = 0

    # initialize process
    config = {
        'initial_ligand': timeline[0][1]}
    receptor = ReceptorCluster(config)
    settings = receptor.default_settings()
    state = settings['state']
    ligand_id = receptor.ligand_id
    state['external'][ligand_id] = timeline[0][1]

    # run simulation
    ligand_vec = []
    receptor_activity_vec = []
    n_methyl_vec = []
    time_vec = []
    ligand_conc=timeline[0][1]
    while time < end_time:
        time += timestep

        # update ligand from timeline
        timeline_steps = [
            index for index, (t, conc) in enumerate(timeline) if t < time]
        if len(timeline_steps) > 0:
            step = timeline_steps[-1]
            ligand_conc = timeline[step][1]
            state['external'][ligand_id] = ligand_conc
            del timeline[0:step+1]

        # run step
        run_step(receptor, state, timestep)
        P_on = state['internal']['chemoreceptor_activity']
        n_methyl = state['internal']['n_methyl']

        # save state
        ligand_vec.append(ligand_conc)
        receptor_activity_vec.append(P_on)
        n_methyl_vec.append(n_methyl)
        time_vec.append(time)

    return {
        'ligand_vec': ligand_vec,
        'receptor_activity_vec': receptor_activity_vec,
        'n_methyl_vec': n_methyl_vec,
        'time_vec': time_vec}

def plot_output(output, out_dir='out', filename='response'):
    ligand_vec = output['ligand_vec']
    receptor_activity_vec = output['receptor_activity_vec']
    n_methyl_vec = output['n_methyl_vec']
    time_vec = output['time_vec']

    # plot results
    cols = 1
    rows = 3
    plt.figure(figsize=(cols * 4, rows * 1.5))

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)

    ax1.plot(time_vec, ligand_vec, 'b')
    ax2.plot(time_vec, receptor_activity_vec, 'b')
    ax3.plot(time_vec, n_methyl_vec, 'b')

    ax1.set_xticklabels([])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(right=False, top=False)
    ax1.set_ylabel("external ligand \n (mM) ", fontsize=10)
    # ax1.set_yscale('log')

    ax2.set_xticklabels([])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(right=False, top=False)
    ax2.set_ylabel("cluster activity \n P(on)", fontsize=10)

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.tick_params(right=False, top=False)
    ax3.set_xlabel("time (s)", fontsize=12)
    ax3.set_ylabel("average \n methylation", fontsize=10)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.2)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'Endres2006_chemoreceptor')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    timeline = get_pulse_timeline()
    output = test_receptor(timeline)
    plot_output(output, out_dir, 'pulse')

    linear_config = {
        'time': 10,
        'slope': 1e-3,
        'speed': 14}
    timeline2 = get_linear_step_timeline(linear_config)
    output2 = test_receptor(timeline2)
    plot_output(output2, out_dir, 'linear')

    exponential_config = {
        'time': 10,
        'base': 1+4e-4,
        'speed': 14}
    timeline3 = get_exponential_step_timeline(exponential_config)
    output3 = test_receptor(timeline3)
    plot_output(output3, out_dir, 'exponential_4e-4')

    exponential_random_config = {
        'time': 60,
        'base': 1+4e-4,
        'speed': 14}
    timeline4 = get_exponential_random_timeline(exponential_random_config)
    output4 = test_receptor(timeline4, 0.1)
    plot_output(output4, out_dir, 'exponential_random')
