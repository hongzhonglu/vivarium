from __future__ import absolute_import, division, print_function

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import convolve

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    convert_to_timeseries)



# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


def make_fields(molecules, n_bins):
    # Create lattice and fill each site with concentrations dictionary
    # Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
    fields = {}
    for molecule_id, conc in molecules.items():
        fields[molecule_id] = np.full((n_bins), conc, dtype=np.float64)
    return fields

def fields_to_locations(fields, n_bins):
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    locations = {}
    for x in range(bins_x):
        for y in range(bins_y):
            concs = {
                mol_id: field[x][y]
                for mol_id, field in fields.items()}
            locations[(x,y)]= concs
    return locations

def locations_to_fields(locations, molecule_ids, n_bins):
    fields = {}
    for molecule_id in molecule_ids:
        field = np.empty((n_bins), dtype=np.float64)
        for x in range(n_bins[0]):
            for y in range(n_bins[1]):
                field[x,y] = locations[(x,y)][molecule_id]
        fields[molecule_id] = field
    return fields

def field_from_locations_series(locations_series, molecule_ids, n_bins, times):
    n_times = len(times)
    field_series = {
        mol_id: np.zeros((n_times, n_bins[0], n_bins[1]), dtype=np.float64)
        for mol_id in molecule_ids}

    for molecule_id in molecule_ids:
        for time_index in range(n_times):
            field = np.empty((n_bins), dtype=np.float64)
            for x in range(n_bins[0]):
                for y in range(n_bins[1]):
                    field[x, y] = locations_series[(x, y)][molecule_id][time_index]
            field_series[molecule_id][time_index] = field
    return field_series



class Diffusion(Process):
    '''
    '''

    def __init__(self, initial_parameters={}):

        # initial state
        initial_state = initial_parameters.get('initial_state', {})
        self.initial_sites = initial_state.get('sites', {})
        self.initial_membrane = initial_state.get('membrane_composition', {})

        # locations
        self.molecules = initial_parameters.get('molecules', {'glc': 1})
        self.molecule_ids = list(self.molecules.keys())

        # membranes
        self.membrane_locations = initial_parameters.get('membranes', [])
        self.channels = initial_parameters.get('channels', {})

        # parameters
        self.length_x = initial_parameters.get('length_x', 2)
        self.length_y = initial_parameters.get('length_y', 1)
        bins_x = initial_parameters.get('bins_x', 2)
        bins_y = initial_parameters.get('bins_y', 1)
        self.n_bins = (bins_x, bins_y)
        self.diffusion = initial_parameters.get('diffusion', 1e-1)

        self.dx = self.length_x / bins_x
        self.dy = self.length_y / bins_y
        self.dx2 = self.dx * self.dy
        self.diffusion_dt = 0.5 * self.dx ** 2 * self.dy ** 2 / (2 * self.diffusion * (self.dx ** 2 + self.dy ** 2))

        # make fields
        fields = make_fields(self.molecules, self.n_bins)
        ports = fields_to_locations(fields, self.n_bins)
        ports.update(
            {'membrane_composition': list(self.channels.keys())})

        parameters = {}
        parameters.update(initial_parameters)

        super(Diffusion, self).__init__(ports, parameters)


    def default_settings(self):
        initial_state = {
            'membrane_composition': self.initial_membrane}
        initial_state.update(self.initial_sites)

        return {
            'state': initial_state}

    def next_update(self, timestep, states):
        locations = {
            state_id: concs
            for state_id, concs in states.items() if state_id is not 'membrane_composition'}
        fields = locations_to_fields(locations, self.molecule_ids, self.n_bins)

        # run diffusion
        field_update = self.run_diffusion(fields, timestep)
        locations_update = fields_to_locations(field_update, self.n_bins)

        return locations_update

    # diffusion functions
    def diffusion_timestep(self, field):
        change_field = self.diffusion * self.diffusion_dt * convolve(field, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_field

    def run_diffusion(self, fields, timestep):
        fields_change = {}
        for mol_id, field in fields.items():
            delta = field.copy()
            t = 0.0
            while t < timestep:
                field += self.diffusion_timestep(field)
                t += self.diffusion_dt
            delta -= field
            fields_change[mol_id] = delta
        return fields_change



# testing functions
def get_two_compartment_config():
    initial_state = {
        'membrane_composition': {
            'porin': 1},
        'sites': {
            (0, 0): {
                'glc': 20.0},
            (1, 0): {
                'glc': 0.0}
        }}

    return {
        'initial_state': initial_state,
        'molecules': {
            'glc': 1
        },
        'membranes': [
            ((0,0),(1,0))
        ],
        'channels':{
            'porin': 1e-1  # diffusion rate through porin
        },
        'length_x': 2,
        'bins_x': 2,
        'length_y': 1,
        'bins_y': 1,
        'diffusion': 1e-1}

def get_cell_config():
    initial_state = {
        'membrane_composition': {
            'porin': 1},
        'sites': {
            (2, 2): {
                'glc': 20.0},
            (6, 2): {
                'glc': 20.0}
        }}

    return {
        'initial_state': initial_state,
        'molecules': {
            'glc': 1},
        'membranes': [],
        'length_x': 10,
        'bins_x': 10,
        'length_y': 4,
        'bins_y': 4,
        'diffusion': 1e-1}

def test_diffusion(config = get_two_compartment_config(), time=10):
    # load process
    diffusion = Diffusion(config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(diffusion)
    return simulate_with_environment(compartment, settings)

def plot_diffusion_output(data, config, out_dir='out', filename='field'):

    n_snapshots = 6

    # parameters
    molecules = config.get('molecules')
    molecule_ids = list(molecules)
    n_fields = len(molecule_ids)
    length_x = config.get('length_x')
    length_y = config.get('length_y')
    bins_x = config.get('bins_x')
    bins_y = config.get('bins_y')
    n_bins = (bins_x, bins_y)

    # data
    times = data.get('time')
    locations_series = {(x,y): data[(x,y)]
        for x in range(bins_x)
        for y in range(bins_y)}

    field_series = field_from_locations_series(locations_series, molecule_ids, n_bins, times)

     # plot fields
    time_vec = times
    plot_steps = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int).tolist()
    snapshot_times = [time_vec[i] for i in plot_steps]

    # make figure
    fig = plt.figure(figsize=(20 * n_snapshots, 10*n_fields))
    grid = plt.GridSpec(n_fields, n_snapshots, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    for mol_idx, mol_id in enumerate(molecule_ids):
        field_data = field_series[mol_id]

        for time_index, time_step in enumerate(plot_steps, 0):
            this_field = field_data[time_index]

            ax = fig.add_subplot(grid[mol_idx, time_index])  # grid is (row, column)
            # ax.title.set_text('time: {:.4f} hr | field: {}'.format(float(time) / 60. / 60., field_id))
            ax.set(xlim=[0, length_x], ylim=[0, length_y], aspect=1)
            ax.set_yticklabels([])
            ax.set_xticklabels([])

            plt.imshow(this_field,
                       # vmin=0,
                       # vmax=15.0,
                       origin='lower',
                       # extent=[0, length_y, 0, length_x],
                       interpolation='nearest',
                       cmap='YlGn')

    plt.colorbar()
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # config = get_two_compartment_config()
    config = get_cell_config()
    saved_data = test_diffusion(config, 10)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_diffusion_output(timeseries, config, out_dir)
