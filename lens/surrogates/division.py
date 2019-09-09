from __future__ import absolute_import, division, print_function

import time
import random
import uuid

from lens.actor.inner import Simulation


DEFAULT_COLOR = [color/255 for color in [255, 0 , 0]]

INITIAL_VOLUME = 1.0

class Division(Simulation):
    ''' '''

    def __init__(self, boot_config):
        self.initial_time = boot_config.get('time', 0.0)
        self.volume = boot_config.get('volume', INITIAL_VOLUME)
        self.local_time = self.initial_time
        self.timestep = 0.5

        self.growth = 0.05
        self.division_volume = random.uniform(1.9 * INITIAL_VOLUME, 2.1 * INITIAL_VOLUME)

        self.environment_change = {}
        self.division = []

    def update_state(self):
        # update state based on internal and external concentrations
        self.volume += self.growth * self.timestep

    def daughter_config(self):
        config1 = {
            'id': str(uuid.uuid4()),
            'time': self.time(),
            'volume': self.volume * 0.5}
        config2 = {
            'id': str(uuid.uuid4()),
            'time': self.time(),
            'volume': self.volume * 0.5}
        return [config1, config2]

    def check_division(self):
        if self.volume >= self.division_volume:
            self.division = self.daughter_config()
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

        time.sleep(0.3)  # pause for better coordination with Lens visualization. TODO: remove this

    def generate_inner_update(self):
        return {
            'volume': self.volume,
            'environment_change': self.environment_change,
            'division': self.division,
            'color': DEFAULT_COLOR,
            }
