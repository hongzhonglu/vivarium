from __future__ import absolute_import, division, print_function

import copy

from vivarium.utils.dict_utils import deep_merge
from vivarium.compartment.process import Process

from vivarium.processes.derive_globals import AVOGADRO


class Environment(Process):
    ''' A minimal, non-spatial environment with volume'''

    defaults = {
        'volume': 1e-12,
        'states': [],
        'environment_port': 'environment',
        'exchange_port': 'exchange',
    }

    def __init__(self, initial_parameters={}):

        self.nAvogadro = AVOGADRO
        volume = initial_parameters.get('volume', self.defaults['volume'])
        environment_states = initial_parameters.get('states', self.defaults['states'])
        self.environment_port = initial_parameters.get('environment_port', self.defaults['environment_port'])
        self.exchange_port = initial_parameters.get('exchange_port', self.defaults['exchange_port'])

        ports = {
            self.environment_port: environment_states,
            self.exchange_port: environment_states}

        parameters = initial_parameters
        super(Environment, self).__init__(ports, parameters)

    def ports_schema(self):
        return {}

    def next_update(self, timestep, states):

        import ipdb;
        ipdb.set_trace()

        update = {}
        # ## apply exchange to environment
        # # get counts, convert to change in concentration
        # if exchange:
        #     delta_counts = exchange.state_for(exchange_ids)
        #     mmol_to_counts = (self.nAvogadro.to('1/mmol') * env_volume).to('L/mmol').magnitude
        #     delta_concs = {mol_id: counts / mmol_to_counts for mol_id, counts in delta_counts.items()}
        #     environment.apply_update(delta_concs)
        #
        #     # reset exchange
        #     reset_exchange = {key: 0 for key in exchange_ids}
        #     exchange.assign_values(reset_exchange)



        return update