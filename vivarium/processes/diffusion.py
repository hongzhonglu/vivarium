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


DIFFUSION_CONSTANT = 1e-2

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])



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

def check_in_set(set_list, set):
    in_set = False
    for edge in set_list:
        if set[0] in edge and set[1] in edge:
            in_set = True
    return in_set

def make_location_network(molecules, n_bins):
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    locations = []
    adjacent_edges = set()

    # make location dict
    for x in range(bins_x):
        for y in range(bins_y):
            locations.append((x,y))

    # make adjacency list
    for x in range(bins_x):
        for y in range(bins_y):
            if y > 0:
                south = ((x, y), (x, y-1))
                if not check_in_set(adjacent_edges, south):
                    adjacent_edges.add(south)
            if y < bins_y-1:
                north = ((x, y), (x, y+1))
                if not check_in_set(adjacent_edges, north):
                    adjacent_edges.add(north)
            if x > 0:
                west = ((x, y), (x-1, y))
                if not check_in_set(adjacent_edges, west):
                    adjacent_edges.add(west)
            if x < bins_x-1:
                east = ((x, y), (x+1, y))
                if not check_in_set(adjacent_edges, east):
                    adjacent_edges.add(east)

    return locations, adjacent_edges


class Diffusion(Process):
    '''
    '''

    def __init__(self, initial_parameters={}):

        # initial state
        initial_state = initial_parameters.get('initial_state', {})
        self.initial_membrane = initial_state.get('membrane_composition', {})
        self.initial_sites = initial_state.get('sites', {})

        # membranes
        self.membrane_locations = initial_parameters.get('membrane_locations', [])
        self.channels = initial_parameters.get('channels', {})

        # parameters
        n_bins = initial_parameters.get('n_bins', (2,1))
        size = initial_parameters.get('size', (2, 1))
        bins_x = n_bins[0]
        bins_y = n_bins[1]
        length_x = size[0]
        length_y = size[1]

        ## make diffusion network
        molecule_ids = initial_parameters.get('molecules', ['glc'])
        diffusion = initial_parameters.get('diffusion', DIFFUSION_CONSTANT)
        locations, edges = make_location_network(molecule_ids, n_bins)
        self.diffusion_network = DiffusionNetwork(locations, edges, molecule_ids, diffusion)

        self.dx = length_x / bins_x
        self.dy = length_y / bins_y
        self.dx2 = self.dx * self.dy
        diffusion_dt = 0.5 * self.dx ** 2 * self.dy ** 2 / (2 * diffusion * (self.dx ** 2 + self.dy ** 2))
        self.diffusion_dt = min(diffusion_dt, 1)

        # make ports from locations and membrane channels
        ports = {
            site: molecule_ids
            for site in locations}
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
        concentrations = {
            state_id: concs
            for state_id, concs in states.items() if state_id is not 'membrane_composition'}

        diffusion_delta = self.diffusion_network.diffusion_delta(concentrations, timestep)

        return diffusion_delta


class DiffusionNetwork(object):

    def __init__(self, nodes, edges, molecule_ids, diffusion):
        self.nodes = nodes
        self.edges = edges
        self.molecule_ids = molecule_ids
        self.diffusion = diffusion

    def diffusion_delta(self, concentrations, timestep):
        diffusion_delta = {
            site: {mol_id: 0 for mol_id in self.molecule_ids}
            for site in self.nodes}

        for edge in self.edges:
            node1 = edge[0]
            node2 = edge[1]
            concs1 = concentrations[node1]
            concs2 = concentrations[node2]
            for mol_id in self.molecule_ids:
                delta1 = self.diffusion * timestep * (concs2[mol_id] - concs1[mol_id])
                delta2 = -delta1
                diffusion_delta[node1][mol_id] += delta1
                diffusion_delta[node2][mol_id] += delta2

        return diffusion_delta


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
        'molecules': ['glc'],
        'membrane_locations': [((0,0),(1,0))],
        'channels':{
            'porin': 1e-1  # diffusion rate through porin
        },
        'n_bins': (2, 1),
        'size': (2e-2, 1e-2),
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
        'molecules': ['glc'],
        'membrane_locations': [],
        'n_bins': (10, 4),
        'size': (10e-1, 4e-1),
        'diffusion': 1e-1}

def test_diffusion(config = get_two_compartment_config(), time=10):
    # load process
    diffusion = Diffusion(config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-2}

    compartment = process_in_compartment(diffusion)
    return simulate_with_environment(compartment, settings)

def plot_diffusion_output(data, config, out_dir='out', filename='field'):
    n_snapshots = 6

    # parameters
    molecules = config.get('molecules')
    molecule_ids = list(molecules)
    n_fields = len(molecule_ids)

    n_bins = config.get('n_bins')
    size = config.get('size')
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    length_x = size[0]
    length_y = size[1]

    # data
    times = data.get('time')
    locations_series = {(x,y): data[(x,y)]
        for x in range(bins_x)
        for y in range(bins_y)}

    field_series = field_from_locations_series(locations_series, molecule_ids, n_bins, times)

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
            # ax.title.set_text('time: {:.4f} hr | field: {}'.format(float(time) / 60. / 60., field_id))
            ax.set(xlim=[0, length_x], ylim=[0, length_y], aspect=1)
            ax.set_yticklabels([])
            ax.set_xticklabels([])

            plt.imshow(np.rot90(this_field),
                       vmin=vmin,
                       vmax=vmax,
                       origin='lower',
                       extent=[0, length_x, 0, length_y],
                       interpolation='nearest',
                       cmap='YlGn')

    # plt.colorbar()
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
    timeseries = convert_to_timeseries(saved_data)
    plot_diffusion_output(timeseries, config, out_dir)
