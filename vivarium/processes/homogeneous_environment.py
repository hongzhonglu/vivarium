from __future__ import absolute_import, division, print_function

import copy

from vivarium.utils.units import units
from vivarium.utils.dict_utils import deep_merge
from vivarium.core.process import Process

from vivarium.processes.derive_globals import AVOGADRO


class HomogeneousEnvironment(Process):
    ''' A minimal, non-spatial environment with volume'''

    defaults = {
        'volume': 1e-12,
        'states': [],
        'environment_port': 'environment',
        'exchange_port': 'exchange',
    }

    def __init__(self, initial_parameters={}):

        volume = initial_parameters.get('volume', self.defaults['volume']) * units.L
        self.mmol_to_counts = (AVOGADRO.to('1/mmol') * volume).to('L/mmol').magnitude

        environment_states = initial_parameters.get('states', self.defaults['states'])
        self.environment_port = initial_parameters.get('environment_port', self.defaults['environment_port'])
        self.exchange_port = initial_parameters.get('exchange_port', self.defaults['exchange_port'])

        ports = {
            self.environment_port: environment_states,
            self.exchange_port: environment_states}

        parameters = initial_parameters
        super(HomogeneousEnvironment, self).__init__(ports, parameters)

    def ports_schema(self):
        schema = {
            self.exchange_port: {
                mol_id: {
                    '_default': 0.0}
                for mol_id in self.ports['exchange']}}
        return schema

    def next_update(self, timestep, states):
        exchange = states[self.exchange_port]  # units: counts

        ## apply exchange to environment
        # get counts, convert to concentration change
        update = {
            self.exchange_port: {},
            self.environment_port: {}}

        for mol_id, delta_count in exchange.items():
            delta_concs = delta_count / self.mmol_to_counts
            if delta_concs != 0:
                update[self.environment_port][mol_id] = delta_concs

                # reset exchange
                update[self.exchange_port][mol_id] = {  # 0.0
                    '_value': 0.0,
                    '_updater': 'set'}

        return update
