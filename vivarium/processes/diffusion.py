from __future__ import absolute_import, division, print_function

import os

import numpy as np
from scipy.ndimage import convolve

from vivarium.compartment.process import Process
from vivarium.utils.dict_utils import deep_merge, tuplify_port_dicts
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    convert_to_timeseries,
    plot_simulation_output)

# constants
DIFFUSION_CONSTANT = 1e3

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


class Diffusion(Process):
    '''
    '''

    def __init__(self, initial_parameters={}):

        # locations
        self.fields = initial_parameters.get('fields', [])
        self.membranes = initial_parameters.get('membranes', [])

        # diffusion parameters
        self.edge_length_x = 2
        self.patches_per_edge_x = 2
        self.edge_length_y = 1
        self.patches_per_edge_y = 1

        self.diffusion = initial_parameters.get('diffusion', DIFFUSION_CONSTANT)

        self.dx = self.edge_length_x / self.patches_per_edge_x
        self.dy = self.edge_length_y / self.patches_per_edge_y
        self.dx2 = self.dx * self.dy
        self.diffusion_dt = 0.5 * self.dx ** 2 * self.dy ** 2 / (2 * self.diffusion * (self.dx ** 2 + self.dy ** 2))

        # # get initial state
        # states = list(self.transcription.keys()) + list(self.translation.keys())
        # null_states = {'internal': {
        #     state_id: 0 for state_id in states}}
        # initialized_states = initial_parameters.get('initial_state', {})
        # self.initial_state = deep_merge(null_states, initialized_states)
        # internal = list(self.initial_state.get('internal', {}).keys())
        # external = list(self.initial_state.get('external', {}).keys())

        import ipdb; ipdb.set_trace()

        ports = {}

        parameters = {}
        parameters.update(initial_parameters)

        super(Diffusion, self).__init__(ports, parameters)

    def default_settings(self):
        default_settings = {}
        return default_settings

    def next_update(self, timestep, states):
        internal_state = states['internal']

        update = {}
        return update

    # diffusion ufnctions
    def diffusion_timestep(self, field, dt):
        ''' calculate concentration changes cause by diffusion'''
        change_field = self.diffusion * dt * convolve(field, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_field

    def run_diffusion(self, timestep):
        for index in range(len(self.fields)):
            molecule = self.fields[index]
            # run diffusion if molecule field is not uniform
            if len(set(molecule.flatten())) != 1:
                t = 0.0
                while t < timestep:
                    molecule += self.diffusion_timestep(molecule, self.diffusion_dt)
                    t += self.diffusion_dt


# testing
def get_cell_config():
    return {}

def test_diffusion(time=10):
    config = get_cell_config()

    # load process
    diffusion = Diffusion(config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(diffusion)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_diffusion(20)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)
