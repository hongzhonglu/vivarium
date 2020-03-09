from __future__ import absolute_import, division, print_function

import os

import numpy as np
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    convert_to_timeseries)



# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])

DIFFUSION_CONSTANT = 5e-1

def gaussian(deviation, distance):
    return np.exp(-np.power(distance, 2.) / (2 * np.power(deviation, 2.)))

def make_gradient(gradient, n_bins, size):
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    length_x = size[0]
    length_y = size[1]
    fields = {}
    
    if gradient.get('type') == 'gaussian':
        """
        gaussian gradient multiplies the base concentration of the given molecule
        by a gaussian function of distance from center and deviation

        'gradient': {
            'type': 'gradient',
            'molecules': {
                'mol_id1':{
                    'center': [0.25, 0.5],
                    'deviation': 30},
                'mol_id2': {
                    'center': [0.75, 0.5],
                    'deviation': 30}
            }},
        """

        for molecule_id, specs in gradient['molecules'].items():
            field = np.ones((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            deviation = specs['deviation']

            for x_patch in range(bins_x):
                for y_patch in range(bins_y):
                    # distance from middle of patch to center coordinates
                    dx = (x_patch + 0.5) * length_x / bins_x - center[0]
                    dy = (y_patch + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    scale = gaussian(deviation, distance)
                    # multiply gradient by scale
                    field[x_patch][y_patch] *= scale
            fields[molecule_id] = field

    elif gradient.get('type') == 'linear':
        """
        linear gradient adds to the base concentration of the given molecule
        as a function of distance from center and slope.

        'gradient': {
            'type': 'linear',
            'molecules': {
                'mol_id1':{
                    'center': [0.0, 0.0],
                    'slope': -10},
                'mol_id2': {
                    'center': [1.0, 1.0],
                    'slope': -5}
            }},
        """

        for molecule_id, specs in gradient['molecules'].items():
            field = np.ones((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            slope = specs['slope']

            for x_patch in range(bins_x):
                for y_patch in range(bins_y):
                    # distance from middle of patch to center coordinates
                    dx = (x_patch + 0.5) * length_x / bins_x - center[0]
                    dy = (y_patch + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    added = distance * slope
                    # add gradient to basal concentration
                    field[x_patch][y_patch] += added
            fields[fields <= 0.0] = 0.0
            fields[molecule_id] = field

    elif gradient.get('type') == 'exponential':
        """
        exponential gradient adds a delta (d) to the base concentration (c)
        of the given molecule as a function of distance  (x) from center and base (b),
        with d=c+x^d.

        'gradient': {
            'type': 'exponential',
            'molecules': {
                'mol_id1':{
                    'center': [0.0, 0.0],
                    'base': 1+2e-4},
                'mol_id2': {
                    'center': [1.0, 1.0],
                    'base': 1+2e-4}
            }},
        """

        for molecule_id, specs in gradient['molecules'].items():
            field = np.ones((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            base = specs['base']

            for x_patch in range(bins_x):
                for y_patch in range(bins_y):
                    dx = (x_patch + 0.5) * length_x / bins_x - center[0]
                    dy = (y_patch + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    added = base ** distance - 1

                    # add to base concentration
                    field[x_patch][y_patch] += added
            fields[fields <= 0.0] = 0.0
            fields[molecule_id] = field

    return fields



class DiffusionField(Process):
    ''''''
    def __init__(self, initial_parameters={}):

        # initial state
        molecule_ids = initial_parameters.get('molecules', ['glc'])
        self.initial_state = initial_parameters.get('initial_state', {})

        # parameters
        n_bins = initial_parameters.get('n_bins', (2,1))
        size = initial_parameters.get('size', (2, 1))

        # diffusion settings
        diffusion = initial_parameters.get('diffusion', DIFFUSION_CONSTANT)
        bins_x = n_bins[0]
        bins_y = n_bins[1]
        length_x = size[0]
        length_y = size[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        self.diffusion = diffusion / dx2
        self.diffusion_dt = 0.01
        # self.diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * self.diffusion * (dx ** 2 + dy ** 2))

        # initialize gradient fields
        gradient = initial_parameters.get('gradient', {})
        if gradient:
            gradient_fields = make_gradient(gradient, n_bins, size)
            self.initial_state.update(gradient_fields)

        # make ports
        ports = {
            'fields': molecule_ids}

        parameters = {}
        parameters.update(initial_parameters)

        super(DiffusionField, self).__init__(ports, parameters)

    def default_settings(self):
        return {
            'state': {
                'fields': self.initial_state}}

    def next_update(self, timestep, states):
        fields = states['fields']
        delta_fields = self.diffuse(fields, timestep)
        return {
            'fields': delta_fields}

    # diffusion functions
    def diffusion_delta(self, field, timestep):
        ''' calculate concentration changes cause by diffusion'''
        field_new = field.copy()
        t = 0.0
        dt = min(timestep, self.diffusion_dt)
        while t < timestep:
            field_new += self.diffusion * dt * convolve(field_new, LAPLACIAN_2D, mode='reflect')
            t += dt

        return field_new - field

    def diffuse(self, fields, timestep):
        delta_fields = {}
        for mol_id, field in fields.items():

            # run diffusion if molecule field is not uniform
            if len(set(field.flatten())) != 1:
                delta = self.diffusion_delta(field, timestep)
            else:
                delta = np.zeros_like(field)
            delta_fields[mol_id] = delta

        return delta_fields



# testing
def plot_field_output(data, config, out_dir='out', filename='field'):
    n_snapshots = 6

    # parameters
    molecules = config.get('molecules', {})
    molecule_ids = list(molecules)
    n_fields = len(molecule_ids)

    # n_bins = config.get('n_bins')
    size = config.get('size')
    length_x = size[0]
    length_y = size[1]

    # data
    times = data.get('time')
    field_series = data.get('fields')

     # plot fields
    time_vec = times
    plot_steps = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int).tolist()

    # make figure
    fig = plt.figure(figsize=(20 * n_snapshots, 10*n_fields))
    grid = plt.GridSpec(n_fields, n_snapshots, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    for mol_idx, mol_id in enumerate(molecule_ids):
        field_data = field_series[mol_id]
        vmin = np.amin(field_data)
        vmax = np.amax(field_data)

        for time_index, time_step in enumerate(plot_steps, 0):
            this_field = field_data[time_index]

            ax = fig.add_subplot(grid[mol_idx, time_index])  # grid is (row, column)
            ax.set(xlim=[0, length_x], ylim=[0, length_y], aspect=1)
            ax.set_yticklabels([])
            ax.set_xticklabels([])

            # plot field
            plt.imshow(np.transpose(this_field),
                       vmin=vmin,
                       vmax=vmax,
                       origin='lower',
                       extent=[0, length_x, 0, length_y],
                       interpolation='nearest',
                       cmap='YlGn')

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def get_field_config():
    n_bins = (10, 10)
    return {
        'molecules': ['glc'],
        'initial_state': {
            'glc': np.random.rand(n_bins[0], n_bins[1])},
        'n_bins': n_bins,
        'size': n_bins}

def get_gaussian_config():
    n_bins = (10, 10)
    return {
        'molecules': ['glc'],
        'n_bins': n_bins,
        'size': n_bins,
        'gradient': {
            'type': 'gaussian',
            'molecules': {
                'glc': {
                    'center': [0.5, 0.5],
                    'deviation': 1},
            }},
    }


def test_diffusion(config, time=10):
    diffusion = DiffusionField(config)
    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(diffusion)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion_field')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_field_config()
    saved_data = test_diffusion(config, 10)
    timeseries = convert_to_timeseries(saved_data)
    plot_field_output(timeseries, config, out_dir, 'field')

    gaussian_config = get_gaussian_config()
    gaussian_data = test_diffusion(gaussian_config, 10)
    gaussian_timeseries = convert_to_timeseries(gaussian_data)
    plot_field_output(gaussian_timeseries, gaussian_config, out_dir, 'gaussian_field')
