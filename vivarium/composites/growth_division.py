from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import (
    initialize_state)

from vivarium.compartment.tree import (
    process_derivers,
    Compartment,
)

from vivarium.compartment.composition import (
    get_derivers,
    simulate_with_environment,
    plot_simulation_output, load_compartment)

# processes
from vivarium.processes.growth import Growth
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.division import (
    Division,
    divide_condition
)
from vivarium.processes.meta_division import MetaDivision
from vivarium.processes.convenience_kinetics import (
    ConvenienceKinetics,
    get_glc_lct_config
)
from vivarium.utils.dict_utils import deep_merge


class GrowthDivision(Compartment):
    def __init__(self, config):
        self.config = config
        self.global_key = ('..', 'global')
        self.external_key = ('..',) + self.config.get('external_key', ('external',))
        self.cells_key = ('..', '..') + self.config.get('cells_key', ('cells',))

        self.transport_config = self.config.get('transport', get_glc_lct_config())
        self.transport_config['global_deriver_config'] = {
            'type': 'globals',
            'source_port': 'global',
            'derived_port': 'global',
            'global_port': self.global_key,
            'keys': []}

    def generate_processes(self, config):
        # declare the processes
        agent_id = config.get('agent_id')

        transport_config = deep_merge(
            config.get('transport', {}),
            self.transport_config)

        division_config = dict(
            config.get('division', {}),
            cell_id=agent_id,
            compartment=self)

        growth = Growth(config.get('growth', {}))
        transport = ConvenienceKinetics(transport_config)
        division = MetaDivision(division_config)
        expression = MinimalExpression(config.get('expression', {}))

        # place processes in layers,
        # let compartment handle derivers
        return {
            'transport': transport,
            'growth': growth,
            'expression': expression,
            'division': division}

    def generate_topology(self, config):
        # make the topology.
        # for each process, map process ports to store ids
        external_key = config.get('external_key', self.external_key)
        global_key = config.get('global_key', self.global_key)
        cells_key = config.get('cells_key', self.cells_key)

        return {
            'transport': {
                'internal': ('cell',),
                'external': external_key,
                'exchange': external_key,
                # 'fluxes': ['flux'], # just for testing
                'fluxes': None,
                'global': global_key},
            'growth': {
                'global': global_key},
            'division': {
                'global': global_key,
                'cells': cells_key},
            'expression': {
                'internal': ('cell',),
                'external': external_key,
                'concentrations': ('cell_concentrations',)}}


def growth_division(config):
    compartment = GrowthDivision(config)
    return compartment.generate({})


def compose_growth_division(config):
    agent = growth_division(config)
    processes = agent['processes']
    topology= agent['topology']

    # add derivers
    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = {}
    all_processes.update(processes)
    all_processes.update(derivers['deriver_processes'])
    topology.update(derivers['deriver_topology'])  # add derivers to the topology


    # initialize the states
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    options = {
        'name': 'growth_division_composite',
        # 'environment_port': BOUNDARY_STATE,
        # 'exchange_port': BOUNDARY_STATE,
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition}

    return {
        'processes': processes,
        'derivers': deriver_processes,
        'states': states,
        'options': options}


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'growth_division_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_growth_division)

    # settings for simulation and plot
    options = compartment.configuration
    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 100,
    }

    plot_settings = {
        'max_rows': 25,
        'skip_ports': ['prior_state']}

    timeseries = simulate_with_environment(compartment, settings)
    plot_simulation_output(timeseries, plot_settings, out_dir)
