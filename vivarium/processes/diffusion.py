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

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


class Diffusion(Process):
    '''
    '''

    def __init__(self, initial_parameters={}):

        # locations
        self.molecules = initial_parameters.get('molecules', {'glc': 1})
        self.molecule_ids = list(self.molecules.keys())
        self.membranes = initial_parameters.get('membranes', [])

        # diffusion parameters
        self.length_x = initial_parameters.get('length_x', 2)
        self.bins_x = initial_parameters.get('bins_x', 2)
        self.length_y = initial_parameters.get('length_y', 1)
        self.bins_y = initial_parameters.get('bins_y', 1)
        self.diffusion = initial_parameters.get('diffusion', 1e3)

        self.dx = self.length_x / self.bins_x
        self.dy = self.length_y / self.bins_y
        self.dx2 = self.dx * self.dy
        self.diffusion_dt = 0.5 * self.dx ** 2 * self.dy ** 2 / (2 * self.diffusion * (self.dx ** 2 + self.dy ** 2))

        # make fields
        fields = self.fill_fields(self.molecules)

        # make ports
        ports = self.make_ports(fields)

        parameters = {}
        parameters.update(initial_parameters)

        super(Diffusion, self).__init__(ports, parameters)

    def make_ports(self, fields):
        ports = {}
        for x in range(self.bins_x):
            for y in range(self.bins_y):
                concs = {
                    mol_id: fields[mol_index][x][y]
                    for mol_index, mol_id in enumerate(self.molecule_ids)}
                ports[(x,y)]= concs
        return ports

    def make_fields(self, ports):
        fields = np.empty((len(self.molecule_ids), self.bins_x, self.bins_y), dtype=np.float64)
        for (x,y), conc_dict in ports.items():
            concs = [conc_dict[mol_id] for mol_id in self.molecule_ids]
            fields[:,x,y] = concs
        return fields

    def fill_fields(self, molecules={}):
        # Create lattice and fill each site with concentrations dictionary
        # Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
        fields = np.empty((len(self.molecule_ids), self.bins_x, self.bins_y), dtype=np.float64)
        for index, molecule_id in enumerate(self.molecule_ids):
            fields[index].fill(molecules[molecule_id])
        return fields

    def default_settings(self):
        default_settings = {
            'state': {},
        }
        return default_settings

    def next_update(self, timestep, states):
        fields = self.make_fields(states)
        self.run_diffusion(fields, timestep)

        import ipdb;
        ipdb.set_trace()

        update = {}
        return update


    # diffusion functions
    def diffusion_timestep(self, field, dt):
        ''' calculate concentration changes cause by diffusion'''
        change_field = self.diffusion * dt * convolve(field, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_field

    def run_diffusion(self, fields, timestep):
        for index in range(len(fields)):
            molecule = fields[index]
            # run diffusion if molecule field is not uniform
            if len(set(molecule.flatten())) != 1:
                t = 0.0
                while t < timestep:
                    molecule += self.diffusion_timestep(molecule, self.diffusion_dt)
                    t += self.diffusion_dt


# testing
def get_cell_config():
    return {
    'molecules': {
        'glc': 1},
    'membranes': [],
    'length_x': 10,
    'bins_x': 10,
    'length_y': 4,
    'bins_y': 4,
    'diffusion': 1e3}

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
