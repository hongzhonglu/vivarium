from __future__ import absolute_import, division, print_function

import os

import numpy as np
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process_with_environment,
    convert_to_timeseries)

from vivarium.processes.diffusion_field import LAPLACIAN_2D


PI = np.pi
DIFFUSION_CONSTANT = 1e-2
init_in_half = True


class MinSystem(Process):
    ''''''
    def __init__(self, initial_parameters={}):

        # cell size parameters
        self.length = initial_parameters.get('length', 10.0)  # micrometers
        self.radius = initial_parameters.get('radius', 0.5)  # micrometers
        self.diameter = 2 * self.radius

        # chemical parameters
        self.diffusion = initial_parameters.get('diffusion', 2.5)  # micrometer**2/sec
        self.k_ADP_ATP = initial_parameters.get('k_ADP_ATP', 1)  # sec**-1. Conversion rate of MinD:ADP to MinD:ATP
        self.k_D = initial_parameters.get('k_D', 0.025)  # micrometer/sec. Spontaneous attachment rate of MinD:ATP to membrane
        self.k_dD = initial_parameters.get('k_dD', 0.0015)  # micrometer**3/sec. Recruitment of MinD:ATP to membrane by attached MinD
        self.k_de = initial_parameters.get('k_de', 0.7 * 3)  # sec**-1. Rate of ATP hydrolysis in MinE:MinD:ATP complex, breaks apart complex and releases phosphate
        self.k_E = initial_parameters.get('k_E', 0.093 * 3)  # micrometer**3/sec. Rate of MinE attachment to membrane-associated MinD:ATP complex

        # grid parameters
        bin_size = initial_parameters.get('bin_size', 0.1)  # micrometers
        self.bins_y = int(self.length / bin_size) 
        self.bins_x = int(self.diameter / bin_size)
        self.bins_mem = self.bins_x + self.bins_y - 2  # membrane goes along x and wraps around to each cap's midpoint. subtract 2 for the corners
        self.dx = self.length / self.bins_x

        # initial concentrations
        self.MinD_conc = initial_parameters.get('MinD_conc', 1000 / (PI * self.radius ** 2))
        self.MinE_conc = initial_parameters.get('MinE_conc', 350 / (PI * self.radius ** 2))

        # make templates
        self.edges_template = np.pad(np.zeros((self.bins_y - 2, self.bins_x - 2)), (1, 1), 'constant', constant_values=(1, 1))

        # indices for the membrane, specifying contact sites along the body of the cylinder, and the caps.
        self.cap_pos1 = self.bins_y / 2 - 1
        self.cap_pos2 = self.bins_y / 2 + self.bins_x - 2


        half_bins_y = int(self.bins_y / 2)
        self.top_membrane = [(i, 0) for i in range(1, half_bins_y)]
        self.top_membrane.extend([(0, i) for i in range(self.bins_x)])
        self.top_membrane.extend([(i, -1) for i in range(1, half_bins_y)])

        self.bottom_membrane = [(i, 0) for i in range(half_bins_y, self.bins_y - 1)]
        self.bottom_membrane.extend([(-1, i) for i in range(self.bins_x)])
        self.bottom_membrane.extend([(i, -1) for i in range(half_bins_y, self.bins_y - 1)])

        # molecule indices
        index = {
            'MinD-ADP[c]': 0,
            'MinD-ATP[c]': 1,
            'MinE[c]': 2,
            'MinD-ATP[m]': 0,
            'MinE-MinD-ATP[m]': 1,
        }
        # molecule_ids = list(index.keys())


        

        # make ports
        ports = {
            'cytoplasm': ['MinD-ADP[c]', 'MinD-ATP[c]', 'MinE[c]'],
            'membrane': ['MinD-ATP[m]', 'MinE-MinD-ATP[m]']
        }

        parameters = {}
        parameters.update(initial_parameters)

        super(MinSystem, self).__init__(ports, parameters)


    def default_settings(self):
        ## Initialize molecular fields

        # cytoplasm
        cytoplasm = {}
        cytoplasm['MinD-ADP[c]'] = np.random.normal(self.MinD_conc, 1, (self.bins_y, self.bins_x))
        cytoplasm['MinD-ATP[c]'] = np.random.normal(self.MinD_conc, 1, (self.bins_y, self.bins_x))
        cytoplasm['MinE[c]'] = np.random.normal(self.MinE_conc, 1, (self.bins_y, self.bins_x))

        if init_in_half:
            half_bins_x = int(self.bins_x / 2)
            # put all MinD in one half of the cytoplasm to speed up time to oscillations
            cytoplasm['MinD-ADP[c]'][:, 0:half_bins_x] *= 2.0
            cytoplasm['MinD-ADP[c]'][:, half_bins_x:self.bins_x] *= 0.0
            cytoplasm['MinD-ATP[c]'][:, 0:half_bins_x] *= 2.0
            cytoplasm['MinD-ATP[c]'][:, half_bins_x:self.bins_x] *= 0.0

        # membrane
        membrane = {}
        membrane['MinD-ATP[m]'] = np.random.normal(self.MinD_conc, 1, (self.bins_mem,))
        membrane['MinE-MinD-ATP[m]'] = np.random.normal(self.MinE_conc, 1, (self.bins_mem,))

        initial_state = {
            'membrane': membrane,
            'cytoplasm': cytoplasm}

        return {
            'state': initial_state}

    def next_update(self, timestep, states):
        cytoplasm = states['cytoplasm']
        membrane = states['membrane']

        # grid
        top_membrane = self.top_membrane
        bottom_membrane = self.bottom_membrane
        
        # parameters
        diffusion = self.parameters['diffusion']
        k_D = self.parameters['k_D']
        k_dD = self.parameters['k_dD']
        k_E = self.parameters['k_E']
        k_de = self.parameters['k_de']
        k_ADP_ATP = self.parameters['k_ADP_ATP']

        bins_x = self.bins_x
        bins_y = self.bins_y
        dx = self.dx

        ## Initialize deltas and diffusion arrays
        d_cytoplasm = {}
        d_membrane = {}
        diffusion_c = {}


        ## Get diffusion
        diffusion_c['MinD-ADP[c]'] = convolve(
            cytoplasm['MinD-ADP[c]'],
            LAPLACIAN_2D,
            mode='reflect') / dx ** 2 * diffusion

        diffusion_c['MinD-ATP[c]'] = convolve(
            cytoplasm['MinD-ATP[c]'],
            LAPLACIAN_2D,
            mode='reflect') / dx ** 2 * diffusion

        diffusion_c['MinE[c]'] = convolve(
            cytoplasm['MinE[c]'],
            LAPLACIAN_2D,
            mode='reflect') / dx ** 2 * diffusion
        
        ## Get values at cytoplasm-membrane contact sites
        # cytoplasm to membrane contacts
        DT_c_top = np.array([cytoplasm['MinD-ATP[c]'][idx] for idx in top_membrane])
        DT_c_bottom = np.array([cytoplasm['MinD-ATP[c]'][idx] for idx in bottom_membrane])
        DT_c_exchange = (DT_c_top + DT_c_bottom) / 2

        E_c_top = np.array([cytoplasm['MinE[c]'][idx] for idx in top_membrane])
        E_c_bottom = np.array([cytoplasm['MinE[c]'][idx] for idx in bottom_membrane])
        exchange_E_c = (E_c_top + E_c_bottom) / 2

        # membrane to cytoplasm contacts
        exchange_DT_m = np.zeros((bins_y, bins_x))
        for mem_idx, cyto_idx in enumerate(top_membrane):
            exchange_DT_m[cyto_idx] = membrane['MinD-ATP[m]'][mem_idx]
        for mem_idx, cyto_idx in enumerate(bottom_membrane):
            exchange_DT_m[cyto_idx] = membrane['MinD-ATP[m]'][mem_idx]

        exchange_EDT_m = np.zeros((bins_y, bins_x))
        for mem_idx, cyto_idx in enumerate(top_membrane):
            exchange_EDT_m[cyto_idx] = membrane['MinE-MinD-ATP[m]'][mem_idx]
        for mem_idx, cyto_idx in enumerate(bottom_membrane):
            exchange_EDT_m[cyto_idx] = membrane['MinE-MinD-ATP[m]'][mem_idx]

        ## Calculate reaction rates
        # rates for cytoplasm
        rxn_1_c = (k_D + k_dD * (exchange_DT_m + exchange_EDT_m)) * self.edges_template * cytoplasm['MinD-ATP[c]']
        rxn_2_c = k_E * exchange_DT_m * cytoplasm['MinE[c]']
        rxn_3_c = k_de * exchange_EDT_m
        rxn_4_c = k_ADP_ATP * cytoplasm['MinD-ADP[c]']

        # rates for membrane
        rxn_1_m = (k_D + k_dD * (membrane['MinD-ATP[m]'] + membrane['MinE-MinD-ATP[m]'])) * DT_c_exchange
        rxn_2_m = k_E * membrane['MinD-ATP[m]'] * exchange_E_c
        rxn_3_m = k_de * membrane['MinE-MinD-ATP[m]']

        ## Get derivatives
        d_cytoplasm['MinD-ADP[c]'] = (diffusion_c['MinD-ADP[c]'] - rxn_4_c + rxn_3_c)
        d_cytoplasm['MinD-ATP[c]'] = (diffusion_c['MinD-ATP[c]'] + rxn_4_c - rxn_1_c)
        d_cytoplasm['MinE[c]'] = (diffusion_c['MinE[c]'] + rxn_3_c - rxn_2_c)
        d_membrane['MinD-ATP[m]'] = (-rxn_2_m + rxn_1_m)
        d_membrane['MinE-MinD-ATP[m]'] = (-rxn_3_m + rxn_2_m)

        return {
            'cytoplasm': d_cytoplasm,
            'membrane': d_membrane}



# testing functions
def get_min_config():

    length = 10.0  # um
    radius = 0.5  # um

    # initial concentrations
    MinD_conc = 1000 / (PI * radius ** 2)  # reported: 1000/um
    MinE_conc = 350 / (PI * radius ** 2)  # reported: 350/um


    return {
        'MinD_conc': MinD_conc,
        'MinE_conc': MinE_conc,
        'length': length,
        'radius': radius,
        'diffusion' : 2.5,  # micrometer**2/sec
        'k_ADP_ATP' : 1,  # sec**-1. Conversion rate of MinD:ADP to MinD:ATP
        'k_D' : 0.025,  # micrometer/sec. Spontaneous attachment rate of MinD:ATP to membrane
        'k_dD' : 0.0015,  # micrometer**3/sec. Recruitment of MinD:ATP to membrane by attached MinD
        'k_de' : 0.7*3,  # sec**-1. Rate of ATP hydrolysis in MinE:MinD:ATP complex, breaks apart complex and releases phosphate
        'k_E' : 0.093*3,  # micrometer**3/sec. Rate of MinE attachment to membrane-associated MinD:ATP complex
        }

def test_min_system(config = get_min_config(), time=10):
    min = MinSystem(config)
    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-2}

    return simulate_process_with_environment(min, settings)

def plot_min_system_output(data, config, out_dir='out', filename='output'):
    n_snapshots = 20
    plot_step_size = 10

    # get saved state
    cytoplasm = data['cytoplasm']
    membrane = data['membrane']


    # import ipdb;
    # ipdb.set_trace()

    fig, axes = plt.subplots(n_snapshots, 6, figsize=(8, 0.5 * n_snapshots))

    for slice in range(n_snapshots):
        # axes[slice, 0].text(0.5, 0.5, str(slice * plot_step_size * DT) + 's')

        # plot cytoplasm
        axes[slice, 1].imshow(cytoplasm['MinD-ATP[c]'][:, :, slice], cmap='YlGnBu', aspect="auto")
        axes[slice, 2].imshow(cytoplasm['MinD-ADP[c]'][:, :, slice], cmap='YlGnBu', aspect="auto")
        axes[slice, 3].imshow(cytoplasm['MinE[c]'][:, :, slice], cmap='YlOrRd', aspect="auto")

        # plot membrane
        axes[slice, 4].plot(membrane['MinD-ATP[m]'][:, slice].T)
        axes[slice, 5].plot(membrane['MinE-MinD-ATP[m]'][:, slice].T)

        # # add vertical lines for cap location
        # axes[slice, 4].axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
        # axes[slice, 4].axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')
        # axes[slice, 5].axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
        # axes[slice, 5].axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')

        axes[slice, 0].axis('off')

        for x in range(1, 6):
            axes[slice, x].set_xticks([])
            axes[slice, x].set_yticks([])

    # axes[0, 1].set_title('[MinD:ATP] cytoplasm', fontsize=6)
    # axes[0, 2].set_title('[MinD:ADP] cytoplasm', fontsize=6)
    # axes[0, 3].set_title('[MinE] cytoplasm', fontsize=6)
    # axes[0, 4].set_title('[MinD:ATP] membrane', fontsize=6)
    # axes[0, 5].set_title('[MinE:MinD:ATP] membrane', fontsize=6)










    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)






if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'min_system')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_min_config()
    saved_data = test_min_system(config, 30)
    timeseries = convert_to_timeseries(saved_data)
    plot_min_system_output(timeseries, config, out_dir)
