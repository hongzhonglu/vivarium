from __future__ import absolute_import, division, print_function

import os

import matplotlib.pyplot as plt
import numpy as np

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    convert_to_timeseries)



DIFFUSION_CONSTANT = 1e-2

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

def make_location_network(n_bins):
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

        # parameters
        n_bins = initial_parameters.get('n_bins', (2,1))
        size = initial_parameters.get('size', (2, 1))
        bins_x = n_bins[0]
        bins_y = n_bins[1]
        length_x = size[0]
        length_y = size[1]

        # make diffusion network
        molecule_ids = initial_parameters.get('molecules', ['glc'])
        diffusion = initial_parameters.get('diffusion', DIFFUSION_CONSTANT)
        locations, edges = make_location_network(n_bins)

        # membranes
        membrane_locations = initial_parameters.get('membrane_locations', [])
        channels = initial_parameters.get('channels', {})

        # diffusion settings
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        # diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * diffusion * (dx ** 2 + dy ** 2))

        diffusion_config = {
            'nodes': locations,
            'edges': edges,
            'molecule_ids': molecule_ids,
            'diffusion': diffusion/dx2,
            'membrane_locations': membrane_locations,
            'channels': channels}

        self.diffusion_network = DiffusionNetwork(diffusion_config)

        # make ports from locations and membrane channels
        ports = {
            site: molecule_ids
            for site in locations}
        ports.update(
            {'membrane_composition': list(channels.keys())})

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
        membrane = states['membrane_composition']

        diffusion_delta = self.diffusion_network.diffusion_delta(concentrations, membrane, timestep)

        return diffusion_delta


class DiffusionNetwork(object):

    def __init__(self, config):
        self.nodes = config.get('nodes')
        self.edges = config.get('edges')
        self.molecule_ids = config.get('molecule_ids')
        self.diffusion = config.get('diffusion')
        self.membrane_locations = config.get('membrane_locations')
        self.channel_diffusion = config.get('channels')

    def diffusion_delta(self, locations, membrane, timestep):
        diffusion_delta = {
            site: {mol_id: 0 for mol_id in self.molecule_ids}
            for site in self.nodes}

        for edge in self.edges:
            node1 = edge[0]
            node2 = edge[1]
            concs1 = locations[node1]
            concs2 = locations[node2]

            if edge in self.membrane_locations or (edge[1], edge[0]) in self.membrane_locations:
                diffusion = 0
                for channel_id, channel_conc in membrane.items():
                    diffusion += channel_conc * self.channel_diffusion[channel_id]
                diffusion = min(diffusion, self.diffusion)
            else:
                diffusion = self.diffusion

            for mol_id in self.molecule_ids:
                delta1 = diffusion * timestep * (concs2[mol_id] - concs1[mol_id])
                delta2 = -delta1
                diffusion_delta[node1][mol_id] += delta1
                diffusion_delta[node2][mol_id] += delta2

        return diffusion_delta



# testing functions
def get_two_compartment_config():
    initial_state = {
        'membrane_composition': {
            'porin': 5},
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
            'porin': 5e-2  # diffusion rate through porin
        },
        'n_bins': (2, 1),
        'size': (2e-2, 1e-2),
        'diffusion': 1e-1}

def get_grid_config():
    initial_state = {
        'membrane_composition': {
            'porin': 1e2},
        'sites': {
            (2, 1): {
                'glc': 20.0},
            (2, 2): {
                'glc': 20.0},
            (3, 1): {
                'glc': 20.0},
            (3, 2): {
                'glc': 20.0},
            (6, 1): {
                'glc': 20.0}
        }}

    return {
        'initial_state': initial_state,
        'molecules': ['glc'],
        'membrane_locations': [
            ((1, 0), (2, 0)),
            ((1, 1), (2, 1)),
            ((1, 2), (2, 2)),
            ((2, 3), (2, 2)),
            ((3, 3), (3, 2)),
            ((4, 3), (4, 2)),
            ((5, 3), (5, 2)),
            ((6, 2), (6, 3)),
        ],
        'channels': {
            'porin': 5e-4  # diffusion rate through porin
        },
        'n_bins': (10, 4),
        'size': (10, 4),
        'diffusion': 2e-1}

def test_diffusion(config = get_two_compartment_config(), time=10):
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
    membrane_locations = config.get('membrane_locations')

    n_bins = config.get('n_bins')
    size = config.get('size')
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    length_x = size[0]
    length_y = size[1]
    bin_size_x = length_x / bins_x
    bin_size_y = length_y / bins_y

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

            # plot membrane
            for membrane in membrane_locations:
                site1 = membrane[0]
                site2 = membrane[1]
                middle = (
                    (site1[0]+site2[0])/2 + 0.5,
                    (site1[1]+site2[1])/2 + 0.5)
                delta = (
                    (site1[1] - site2[1]),
                    (site1[0] - site2[0]))  # swap x and y
                point_0 = (
                    middle[0] + delta[0]/2,
                    middle[1] + delta[1]/2)
                point_1 = (
                    (middle[0] - delta[0]/2),
                    (middle[1] - delta[1]/2))
                x_0 = point_0[0] * bin_size_x
                y_0 = point_0[1] * bin_size_y
                x_1 = point_1[0] * bin_size_x
                y_1 = point_1[1] * bin_size_y

                plt.plot([x_0, x_1], [y_0, y_1], linewidth=5, color='r')

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_two_compartment_config()
    saved_data = test_diffusion(config, 20)
    timeseries = convert_to_timeseries(saved_data)
    plot_diffusion_output(timeseries, config, out_dir, '2_sites')

    config = get_grid_config()
    saved_data = test_diffusion(config, 20)
    timeseries = convert_to_timeseries(saved_data)
    plot_diffusion_output(timeseries, config, out_dir, 'grid')
