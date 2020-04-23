from __future__ import absolute_import, division, print_function

import os

import numpy as np
from scipy import constants
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import simulate_process


# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])
AVOGADRO = constants.N_A

AGENT_KEYS = ['location', 'exchange', 'local_environment']

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

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    # distance from middle of bin to center coordinates
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    scale = gaussian(deviation, distance)
                    # multiply gradient by scale
                    field[x_bin][y_bin] *= scale
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

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    # distance from middle of bin to center coordinates
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    added = distance * slope
                    # add gradient to basal concentration
                    field[x_bin][y_bin] += added
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

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    added = base ** distance - 1

                    # add to base concentration
                    field[x_bin][y_bin] += added
            fields[fields <= 0.0] = 0.0
            fields[molecule_id] = field

    return fields



class DiffusionField(Process):
    '''
    Diffusion in 2-dimensional fields of molecules, with agent locations for uptake and secretion.

    Notes:
    - Diffusion constant of glucose in 0.5 and 1.5 percent agarose gel = ~6 * 10^-10 m^2/s
        (Weng et al. 2005. Transport of glucose and poly(ethylene glycol)s in agarose gels).
    - Conversion to micrometers: 6 * 10^-10 m^2/s = 600 micrometers^2/s.

    '''

    defaults = {
        'molecules': ['glc'],
        'initial_state': {},
        'n_bins': [10, 10],
        'size': [10, 10],
        'depth': 3000.0,  # um
        'diffusion': 5e-1,
        'gradient': {},
        'agents': {},
    }

    def __init__(self, initial_parameters={}):

        # initial state
        self.molecule_ids = initial_parameters.get('molecules', self.defaults['molecules'])
        self.initial_state = initial_parameters.get('initial_state', self.defaults['initial_state'])

        # parameters
        self.n_bins = initial_parameters.get('n_bins', self.defaults['n_bins'])
        self.size = initial_parameters.get('size', self.defaults['size'])
        depth = initial_parameters.get('depth', self.defaults['depth'])

        # diffusion
        diffusion = initial_parameters.get('diffusion', self.defaults['diffusion'])
        bins_x = self.n_bins[0]
        bins_y = self.n_bins[1]
        length_x = self.size[0]
        length_y = self.size[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        self.diffusion = diffusion / dx2
        self.diffusion_dt = 0.01
        # self.diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * self.diffusion * (dx ** 2 + dy ** 2))

        # volume, to convert between counts and concentration
        total_volume = (depth * length_x * length_y) * 1e-15 # (L)
        self.bin_volume = total_volume / (length_x * length_y)

        # initialize gradient fields
        gradient = initial_parameters.get('gradient', self.defaults['gradient'])
        if gradient:
            gradient_fields = make_gradient(gradient, self.n_bins, self.size)
            self.initial_state.update(gradient_fields)

        # agents
        self.initial_agents = initial_parameters.get('agents', self.defaults['agents'])

        # make ports
        ports = {
             'fields': self.molecule_ids,
             'agents': ['agents']}

        parameters = {}
        parameters.update(initial_parameters)

        super(DiffusionField, self).__init__(ports, parameters)

    def default_settings(self):
        state = {
            'fields': self.initial_state,
            'agents': {'agents': self.initial_agents}
        }
        schema = {'agents': {'agents': {'updater': 'merge'}}}
        default_emitter_keys = {
            port_id: keys for port_id, keys in self.ports.items()}

        return {
            'state': state,
            'schema': schema,
            'emitter_keys': default_emitter_keys,
        }

    def next_update(self, timestep, states):
        fields = states['fields'].copy()
        agents = states['agents']['agents']

        # uptake/secretion from agents
        delta_exchanges = self.apply_exchanges(agents)
        for field_id, delta in delta_exchanges.items():
            fields[field_id] += delta

        # diffuse field
        delta_fields = self.diffuse(fields, timestep)

        # get each agent's local environment
        local_environments = self.get_local_environments(agents, fields)
        agent_update = {
            agent_id: {'local_environment': local_env}
                for agent_id, local_env in local_environments.items()}

        return {
            'fields': delta_fields,
            'agents': {'agents': agent_update}}

    def count_to_concentration(self, count):
        return count / (self.bin_volume * AVOGADRO)

    def get_bin_site(self, location):
        bin = np.array([
            location[0] * self.n_bins[0] / self.size[0],
            location[1] * self.n_bins[1] / self.size[1]])
        bin_site = tuple(np.floor(bin).astype(int))
        return bin_site

    def get_single_local_environments(self, specs, fields):
        bin_site = self.get_bin_site(specs['location'])
        local_environment = {}
        for mol_id, field in fields.items():
            local_environment[mol_id] = field[bin_site]
        return local_environment

    def get_local_environments(self, agents, fields):
        local_environments = {}
        for agent_id, specs in agents.items():
            local_environments[agent_id] = self.get_single_local_environments(specs, fields)
        return local_environments

    def apply_single_exchange(self, delta_fields, specs):
        exchange = specs.get('exchange', {})
        bin_site = self.get_bin_site(specs['location'])

        for mol_id, count in exchange.items():
            concentration = self.count_to_concentration(count)
            delta_fields[mol_id][bin_site[0], bin_site[1]] += concentration

    def apply_exchanges(self, agents):
        # initialize delta_fields with zero array
        delta_fields = {
            mol_id: np.zeros((self.n_bins[0], self.n_bins[1]), dtype=np.float64)
            for mol_id in self.molecule_ids}

        # apply exchanges to delta_fields
        for agent_id, specs in agents.items():
            self.apply_single_exchange(delta_fields, specs)

        return delta_fields

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

def get_random_field_config(n_bins=(10, 10)):
    return {
        'molecules': ['glc'],
        'initial_state': {
            'glc': np.random.rand(n_bins[0], n_bins[1])},
        'n_bins': n_bins,
        'size': n_bins}

def get_gaussian_config(n_bins=(10, 10)):
    return {
        'molecules': ['glc'],
        'n_bins': n_bins,
        'size': n_bins,
        'gradient': {
            'type': 'gaussian',
            'molecules': {
                'glc': {
                    'center': [0.5, 0.5],
                    'deviation': 1}}}}

def get_secretion_agent_config(molecules=['glc'], n_bins=[10, 10]):
    agent = {
        'location': [
                np.random.uniform(0, n_bins[0]),
                np.random.uniform(0, n_bins[1])],
        'exchange': {
            mol_id: 1e2 for mol_id in molecules}}
    agents = {'1': agent}

    config = {
        'molecules': molecules,
        'n_bins': n_bins,
        'size': n_bins,
        'agents': agents}

    return exchange_agent_config(config)

def exchange_agent_config(config):
    molecules = config['molecules']
    size = config['size']
    n_bins = config['n_bins']
    agents = config['agents']
    return {
        'molecules': molecules,
        'initial_state': {
            mol_id: np.ones((n_bins[0], n_bins[1]))
            for mol_id in molecules},
        'n_bins': n_bins,
        'size': size,
        'agents': agents}

def test_diffusion_field(config=get_gaussian_config(), time=10):
    diffusion = DiffusionField(config)
    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}
    return simulate_process(diffusion, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'diffusion_field')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_random_field_config()
    timeseries = test_diffusion_field(config, 10)
    plot_field_output(timeseries, config, out_dir, 'random_field')

    gaussian_config = get_gaussian_config()
    gaussian_timeseries = test_diffusion_field(gaussian_config, 10)
    plot_field_output(gaussian_timeseries, gaussian_config, out_dir, 'gaussian_field')

    secretion_config = get_secretion_agent_config()
    secretion_timeseries = test_diffusion_field(secretion_config, 10)
    plot_field_output(secretion_timeseries, secretion_config, out_dir, 'secretion')
