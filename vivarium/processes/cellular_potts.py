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

        # configure CPM
        cpm_config = {
            'n_initial': n_initial,
            'grid_size': grid_size}
        self.cpm = CPM(cpm_config)

        # make ports
        ports = {
            agent_id: [
                'volume',
                'vol_targetume']
                for agent_id in self.cpm.agent_ids}

        # parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(CellularPotts, self).__init__(ports, parameters)

    def default_settings(self):
        initial_state = {agent_id: {
            'volume': volume,
            'vol_targetume': 10}  # TODO -- configure this
            for agent_id, volume in self.cpm.get_agent_volumes().items()}

        return {
            'state': initial_state}

    def next_update(self, timestep, states):

        self.cpm.update()
        import ipdb; ipdb.set_trace()

        return {}


class CPM(object):

    def __init__(self, config):

        # CPM parameters
        self.temperature = config.get('temperature', 1)
        self.lambda_volume = 40
        self.adhesion_matrix = [[60, 60], [60, 1]]

        # make the grid
        self.grid_size = config.get('grid_size')
        self.grid = np.zeros(self.grid_size, dtype=np.int)
        self.n_sites = self.grid.size

        # make agents, place in grid
        n_initial = config.get('n_initial')
        self.agent_ids = [
            agent_id for agent_id in range(1, n_initial+1)]

        for agent_id in self.agent_ids:
            # pick a site, fill if it is empty
            filled = False
            while not filled:
                x, y = self.random_site()
                if self.grid[x][y] == 0:
                    self.grid[x][y] = agent_id
                    filled = True

    def random_site(self):
        x = random.randint(0, self.grid_size[0] - 1)
        y = random.randint(0, self.grid_size[1] - 1)
        return (x, y)

    def random_neighbor(self, x, y):
        """
        return a random neighbor, without wrapping
        """
        if x == self.grid_size[0]:
            # exclude north
            choice = random.choice([1, 2, 3])
        elif x == 0:
            # exclude south
            choice = random.choice([0, 2, 3])
        elif y == self.grid_size[1]:
            # exclude east
            choice = random.choice([0, 1, 3])
        elif y == 0:
            # exclude west
            choice = random.choice([0, 1, 2])
        else:
            choice = random.choice([0, 1, 2, 3])

        if choice == 0:
            return (x, y+1) # north
        elif choice == 1:
            return (x, y-1) # south
        elif choice == 2:
            return (x+1, y) # east
        elif choice == 3:
            return (x-1, y) # west

    def get_agent_volumes(self):
        volumes = {}
        for agent_id in self.agent_ids:
            n_sites = np.count_nonzero(self.grid == agent_id)
            volumes[agent_id] = n_sites
        return volumes

    def effective_energy_delta(self, volume, adhesion):
        vol_target = volume.get('target')
        vol_before = volume.get('before')
        vol_after = volume.get('after')
        adhesion_before = adhesion.get('before')
        adhesion_after = adhesion.get('after')

        H_before = adhesion_before + self.lambda_volume * ((vol_before - vol_target)**2)
        H_after = adhesion_after + self.lambda_volume * ((vol_after - vol_target)**2)

        return H_after - H_before

    def boltzmann_acceptance_function(self, energy):
        accept = False
        if energy < 0:
            accept = True
        elif np.rand() < np.exp(-energy / self.temperature):
            accept = True

        return accept

    def mutate(self, grid):
        """
        choose a random site and a random neighboring site.
        If they have the same values, change the site to its neighbor's value.
        If a successful mutation is made, return True
        """
        x, y = self.random_site()
        n_x, n_y = self.random_neighbor(x, y)

        value = grid[x][y]
        n_value = grid[n_x][n_y]
        if value != n_value:
            grid[x][y] = n_value
            return True
        else:
            return False

    def update(self):
        """
        Metropolis Monte Carlo.
        Attempt as many updates as there are sites in the grid
        """
        grid_copy = self.grid.copy()

        for update in range(self.n_sites):
            if self.mutate(grid_copy):

                import ipdb; ipdb.set_trace()

                dH = self.effective_energy_delta()

                if self.boltzmann_acceptance_function(dH):
                    # update grid
                    self.grid = grid_copy






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

