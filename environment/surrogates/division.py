from __future__ import absolute_import, division, print_function

import time

from agent.inner import CellSimulation


DEFAULT_COLOR = [color/255 for color in [255, 0 , 0]]

INITIAL_VOLUME = 1.0

class Division(CellSimulation):
    ''' '''

    def __init__(self, state):
        self.initial_time = state.get('time', 0.0)
        self.volume = state.get('volume', INITIAL_VOLUME)
        self.local_time = 0.0
        self.timestep = 1.0

        self.growth = 0.02
        self.division_volume = 2 * INITIAL_VOLUME

        self.environment_change = {}
        self.division = []

        print('initial volume: ' + str(self.volume))


    def update_state(self):
        # update state based on internal and external concentrations
        self.volume += self.growth * self.timestep

    def check_division(self):
        # update division state based on time since initialization
        if self.volume >= self.division_volume:
            daughter_config = {
                'time': self.local_time,
                'volume': self.volume/2}

            self.division = [daughter_config, daughter_config]

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

            if self.division:
                break

        time.sleep(0.2)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': DEFAULT_COLOR,
            }
