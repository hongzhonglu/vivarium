from __future__ import absolute_import, division, print_function

import os
import random

import numpy as np

from vivarium.compartment.process import Process
from vivarium.compartment.composition import process_in_compartment, simulate_with_environment



class CellularPotts(Process):
    def __init__(self, initial_parameters={}):

        grid_size = initial_parameters.get('grid_size', (10, 10))
        n_initial = initial_parameters.get('n_agents', 1)

        cpm_config = {
            'n_initial': n_initial,
            'grid_size': grid_size}
        self.cpm = CPM(cpm_config)


        ports = {agent_id: ['volume'] for agent_id in self.cpm.agent_ids}

        # parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(CellularPotts, self).__init__(ports, parameters)

    def default_settings(self):
        initial_state = {agent_id: {
            'volume': volume}
            for agent_id, volume in self.cpm.get_agent_volumes().items()}

        return {
            'state': initial_state}

    def next_update(self, timestep, states):


        import ipdb; ipdb.set_trace()

        return {}


class CPM(object):

    def __init__(self, config):
        n_initial = config.get('n_initial')
        self.agent_ids = [
            agent_id for agent_id in range(1, n_initial+1)]
        grid_size = config.get('grid_size')

        # make the grid
        self.grid = np.zeros(grid_size, dtype=np.int)
        for agent_id in self.agent_ids:

            # pick a site, fill if it is empty
            filled = False
            while not filled:
                x = random.randint(0,grid_size[0]-1)
                y = random.randint(0,grid_size[1]-1)
                if self.grid[x][y] == 0:
                    self.grid[x][y] = agent_id
                    filled = True

    def get_agent_volumes(self):
        volumes = {}
        for agent_id in self.agent_ids:
            n_sites = np.count_nonzero(self.grid == agent_id)
            volumes[agent_id] = n_sites
        return volumes

    def update(self):
        pass






# test functions
def get_cpm_config():
    config = {
        'n_agents': 1,
        'grid_size': (10, 10)
    }

    return config

def test_CPM(cpm_config = get_cpm_config(), time=10):

    # load process
    expression = CellularPotts(cpm_config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12}

    compartment = process_in_compartment(expression)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'cellular_potts')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cpm_config = get_cpm_config()
    saved_data = test_CPM(cpm_config, 20)

