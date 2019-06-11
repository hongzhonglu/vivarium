from __future__ import absolute_import, division, print_function

import time

from agent.inner import CellSimulation


DEFAULT_COLOR = [color/255 for color in [255, 0 , 0]]

class Division(CellSimulation):
    ''' '''

    def __init__(self, state):
        self.initial_time = state.get('time', 0.0)
        self.local_time = 0.0
        self.timestep = 1.0

        self.environment_change = {}
        self.volume = 1.0
        self.growth = 0.05
        self.division_volume = 2 * self.volume

        self.division = []


    def update_state(self):
        # update state based on internal and external concentrations
        self.volume += self.growth * self.timestep

    def check_division(self):
        # update division state based on time since initialization
        if self.volume >= self.division_volume:
            self.division = [{'time': self.local_time}, {'time': self.local_time}]

        return self.division

    def time(self):
        return self.local_time

    def apply_outer_update(self, update):
        self.external_concentrations = update['concentrations']
        self.environment_change = {}
        for molecule in self.external_concentrations.iterkeys():
            self.environment_change[molecule] = 0

    def run_incremental(self, run_until):
        while self.time() < run_until:
            self.local_time += self.timestep
            self.update_state()
            self.check_division()

        time.sleep(0.2)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': DEFAULT_COLOR,
            }
