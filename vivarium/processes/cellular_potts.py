from __future__ import absolute_import, division, print_function

import os
import random

import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment)



class CellularPotts(Process):
    """
    Cellular Potts model
    """

    defaults = {
        'grid_size': (10, 10),
        'n_agents': 1,
        'target_area': 10
    }

    def __init__(self, initial_parameters={}):

        grid_size = initial_parameters.get('grid_size', self.defaults['grid_size'])
        n_initial = initial_parameters.get('n_agents', self.defaults['n_agents'])
        self.init_target_area = initial_parameters.get('target_area', self.defaults['target_area'])

        # configure CPM
        cpm_config = {
            'n_initial': n_initial,
            'grid_size': grid_size,
            'target_area': self.init_target_area}
        self.cpm = CPM(cpm_config)

        # animate (for debugging)
        self.animate = initial_parameters.get('animate', False)
        if self.animate:
            plt.ion()
            self.animate_frame()

        # make ports
        ports = {
            agent_id: [
                'area',
                'area_target']
                for agent_id in self.cpm.agent_ids}

        # parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(CellularPotts, self).__init__(ports, parameters)

    def default_settings(self):
        initial_state = {agent_id: {
            'area': area,
            'area_target': self.init_target_area}  # TODO -- configure this
            for agent_id, area in self.cpm.get_agents_areas(self.cpm.grid).items()}

        return {
            'state': initial_state}

    def animate_frame(self):
        if self.animate:
            plt.imshow(self.cpm.grid)
            plt.axis('off')
            plt.pause(0.0001)

    def next_update(self, timestep, states):

        area_target = {
            agent_id: states[agent_id]['area_target']
            for agent_id in self.cpm.agent_ids}  # TODO -- get target areas from state
        self.cpm.update_target_areas(area_target)
        self.cpm.update()
        self.animate_frame()

        return {}



class CPM(object):

    def __init__(self, config):
        # CPM parameters
        self.temperature = config.get('temperature', 1e2)
        self.area_constant = config.get('area_constant', 40)
        self.adhesion_baseline = 60
        self.adhesion_matrix = np.array([])

        # make the grid
        self.grid_size = config.get('grid_size')
        self.grid = np.zeros(self.grid_size, dtype=np.int)
        self.n_sites = self.grid.size

        # make agents, place in grid
        n_initial = config.get('n_initial')
        target_area = config.get('target_area')
        self.agent_ids = [
            agent_id for agent_id in range(1, n_initial+1)]

        for agent_id in self.agent_ids:
            # pick a site, fill if it is empty
            filled = False
            while not filled:
                x, y = self.random_site()
                if self.grid[x][y] == 0:
                    self.grid[x][y] = agent_id
                    neighbors = self.neighbor_sites((x, y))

                    # set all neighbors to same agent_id
                    for neighbor in neighbors:
                        self.grid[neighbor[0],neighbor[1]] = agent_id
                    filled = True

        # target areas
        self.target_areas = {
            agent_id: target_area for agent_id in self.agent_ids}

        # make adhesion matrix for agent_ids and medium
        self.update_adhesion_matrix()

    def update_adhesion_matrix(self):
        self.adhesion_matrix = np.full((len(self.agent_ids)+1, len(self.agent_ids)+1), self.adhesion_baseline)
        for agent_id in self.agent_ids:
            self.adhesion_matrix[agent_id, agent_id] = 1  # TODO -- set self-adhesion

    def random_site(self):
        x = random.randint(0, self.grid_size[0] - 1)
        y = random.randint(0, self.grid_size[1] - 1)
        return (x, y)

    def neighbor_sites(self, site):
        """
        return the (x, y) of all neighboring sites
        """
        x, y = site
        neighbors = [
            (n_x, n_y)
            for n_x in range(x-1,x+2)
            for n_y in range(y-1,y+2)
            if n_x>=0 and n_x<self.grid_size[0] and n_y>=0 and n_y<self.grid_size[1] and (n_x, n_y) != site]
        return neighbors

    def random_neighbor(self, site):
        """
        return a random neighbor, without wrapping
        """
        neighbors = self.neighbor_sites(site)
        return random.choice(neighbors)

    def get_agent_area(self, agent_id, grid):
        return np.count_nonzero(grid == agent_id)

    def get_agents_areas(self, grid):
        areas = {}
        for agent_id in self.agent_ids:
            areas[agent_id] = self.get_agent_area(agent_id, grid)
        return areas

    def inverse_kronecker_delta(self, value1, value2):
        """
        Returns 0 if the values are the same, and 1 if they are different.
        Keeps neighboring sites with the same value from contributing to the effective energy
        """
        if value1 == value2:
            return 0
        else:
            return 1

    def neighbor_values(self, site, grid):
        neighbors = self.neighbor_sites(site)
        values = [grid[site] for site in neighbors]
        return values

    def get_interactions(self, site, grid):
        interactions = 0
        site_value = grid[site]
        neighbor_values = self.neighbor_values(site, grid)
        for value in neighbor_values:
            interactions += self.adhesion_matrix[site_value, value] * self.inverse_kronecker_delta(site_value, value)
        return interactions

    def effective_energy(self, grid):
        hamiltonian = 0.0
        for x in range(self.grid_size[0]):
            for y in range(self.grid_size[1]):
                site = (x, y)

                # interaction constraints
                hamiltonian += self.get_interactions(site, grid)

                # area constraints
                areas = self.get_agents_areas(grid)
                area_constraints = [
                    self.area_constant*(areas[agent_id] - self.target_areas[agent_id])**2
                    for agent_id in self.agent_ids]
                hamiltonian += sum(area_constraints)

        return hamiltonian

    def boltzmann_acceptance_function(self, energy):
        accept = False
        if energy < 0:
            accept = True
        elif random.random() < np.exp(-energy / self.temperature):
            accept = True
        return accept

    def mutate(self, grid):
        """
        choose a random site and a random neighboring site.
        If they have the same values, change the site to its neighbor's value.
        If a successful mutation is made, return True
        """
        x, y = self.random_site()
        n_x, n_y = self.random_neighbor((x, y))

        value = grid[x][y]
        n_value = grid[n_x][n_y]
        if value != n_value:
            grid[x][y] = n_value
            return True
        else:
            return False

    def update_target_areas(self, areas):
        self.target_areas.update(areas)

    def update(self):
        """
        Metropolis Monte Carlo.
        Attempt as many updates as there are sites in the grid
        """
        H_before = self.effective_energy(self.grid)
        for update in range(self.n_sites):
            grid_copy = self.grid.copy()
            if self.mutate(grid_copy):
                H_after = self.effective_energy(grid_copy)

                # change in effective energy
                dH = H_after - H_before
                if self.boltzmann_acceptance_function(dH):
                    self.grid = grid_copy
                    H_before = H_after



# test functions
def get_cpm_config():
    return {
        'n_agents': 4,
        'grid_size': (30, 30),
        'target_area': 25,
        'animate': True}

def get_cpm_minimum_config():
    return {
        'n_agents': 2,
        'grid_size': (20, 20),
        'target_area': 10,
        'animate': True}

def run_CPM(cpm_config = get_cpm_minimum_config(), time=5):
    # load process
    expression = CellularPotts(cpm_config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_area': 1e-12}

    compartment = process_in_compartment(expression)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'cellular_potts')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cpm_config = get_cpm_config()
    saved_data = run_CPM(cpm_config, 30)

