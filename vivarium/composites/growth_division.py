from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import (
    initialize_state,
    BOUNDARY_STATE)
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
from vivarium.processes.convenience_kinetics import (
    ConvenienceKinetics,
    get_glc_lct_config
)



def growth_division(config):

    # declare the processes
    transport = ConvenienceKinetics(config.get('transport', get_glc_lct_config()))
    growth = Growth(config.get('growth', {}))
    division = Division(config.get('division', {}))
    expression = MinimalExpression(config.get('expression', {}))

    # place processes in layers
    processes = {
        'transport': transport,
        'growth': growth,
        'expression': expression,
        'division': division}

    external_key = '../' + config.get('external_key', 'external')
    cells_key = '../' + config.get('cells_key', 'cells')

    # make the topology.
    # for each process, map process ports to store ids
    topology = {
        'transport': {
            'internal': 'cell',
            'external': external_key,
            'exchange': external_key,
            'fluxes': None, # instead of 'null'
            'global': 'global'},
        'growth': {
            'global': 'global'},
        'division': {
            'global': 'global',
            'cells': cells_key},
        'expression': {
            'internal': 'cell',
            'external': external_key,
            'concentrations': 'cell_concentrations'}}

    return {
        'processes': processes,
        'topology': topology}


def compose_growth_division(config):
    agent = growth_division(config)
    processes = agent['processes']
    topology= agent['topology']

    # add derivers
    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])  # add derivers to the topology


    # initialize the states
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    options = {
        'name': 'growth_division_composite',
        'environment_port': BOUNDARY_STATE,
        'exchange_port': BOUNDARY_STATE,
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
