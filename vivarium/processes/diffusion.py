from __future__ import absolute_import, division, print_function

import os

import numpy as np
from scipy.ndimage import convolve

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    convert_to_timeseries)



# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])



class Diffusion(Process):
    '''
    '''

    def __init__(self, initial_parameters={}):

        # initial state
        initial_state = initial_parameters.get('initial_state', {})
        self.initial_sites = initial_state.get('sites', {})
        self.initial_membrane = initial_state.get('membrane', {})

        # locations
        self.molecules = initial_parameters.get('molecules', {'glc': 1})
        self.molecule_ids = list(self.molecules.keys())

        # membranes
        self.membranes = initial_parameters.get('membranes', [])
        self.channels = initial_parameters.get('channels', {})

        # parameters
        self.length_x = initial_parameters.get('length_x', 2)
        self.bins_x = initial_parameters.get('bins_x', 2)
        self.length_y = initial_parameters.get('length_y', 1)
        self.bins_y = initial_parameters.get('bins_y', 1)
        self.diffusion = initial_parameters.get('diffusion', 1e-1)

        self.dx = self.length_x / self.bins_x
        self.dy = self.length_y / self.bins_y
        self.dx2 = self.dx * self.dy
        self.diffusion_dt = 0.5 * self.dx ** 2 * self.dy ** 2 / (2 * self.diffusion * (self.dx ** 2 + self.dy ** 2))

        # make fields
        fields = self.fill_fields(self.molecules)

        # make ports from fields
        ports = self.make_ports(fields)
        ports.update(
            {'membrane': list(self.channels.keys())})

        parameters = {}
        parameters.update(initial_parameters)

        super(Diffusion, self).__init__(ports, parameters)


    def default_settings(self):
        initial_state = {
            'membrane': self.initial_membrane}
        initial_state.update(self.initial_sites)

        return {
            'state': initial_state}

    def next_update(self, timestep, states):
        sites = {
            state_id: concs
            for state_id, concs in states.items() if state_id is not 'membrane'}
        fields = self.make_fields(sites)
        field_update = self.run_diffusion(fields, timestep)
        ports_update = self.make_ports(field_update)

        return ports_update


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



    # diffusion functions
    def diffusion_timestep(self, field, dt):
        ''' calculate concentration changes cause by diffusion'''
        change_field = self.diffusion * dt * convolve(field, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_field

    def run_diffusion(self, fields, timestep):
        change_field = np.zeros((len(self.molecule_ids), self.bins_x, self.bins_y), dtype=np.float64)
        for index in range(len(fields)):
            field = fields[index]
            # run diffusion if field is not uniform
            if len(set(field.flatten())) != 1:
                t = 0.0
                while t < timestep:
                    change_field[index] += self.diffusion_timestep(field, self.diffusion_dt)
                    t += self.diffusion_dt

        return change_field


# testing
def get_two_compartment_config():
    initial_state = {
        'membrane': {
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
    return {
        'molecules': {
            'glc': 1},
        'membranes': [],
        'length_x': 10,
        'bins_x': 10,
        'length_y': 4,
        'bins_y': 4,
        'diffusion': 1e-1}

def test_diffusion(time=10):
    # config = get_cell_config()
    config = get_two_compartment_config()

    # load process
    diffusion = Diffusion(config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(diffusion)
    return simulate_with_environment(compartment, settings)


def plot_diffusion_output(timeseries, settings, out_dir):
    pass

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_diffusion(20)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)

    import ipdb; ipdb.set_trace()

    plot_diffusion_output(timeseries, {}, out_dir)
