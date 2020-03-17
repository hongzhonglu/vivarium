from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import (
    initialize_state,
    load_compartment)
from vivarium.compartment.composition import (
    get_derivers,
    simulate_with_environment,
    plot_simulation_output)

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



def compose_growth_division(config):

    # declare the processes
    transport = ConvenienceKinetics(get_glc_lct_config())
    growth = Growth(config)
    division = Division(config)
    expression = MinimalExpression(config)

    # place processes in layers
    processes = [
        {'transport': transport,
         'growth': growth,
         'expression': expression},
        {'division': division}]

    # make the topology.
    # for each process, map process ports to store ids
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'exchange',
            'fluxes': 'null',
            'global': 'global'},
        'growth': {
            'global': 'global'},
        'division': {
            'global': 'global'},
        'expression': {
            'internal': 'cell',
            'external': 'environment',
            'concentrations': 'cell_concs'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])  # add deriver processes
    topology.update(derivers['deriver_topology'])  # add deriver topology

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'growth_division_composite',
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition}

    return {
        'processes': processes,
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
        'emit_timeseries': True,}

    plot_settings = {
        'max_rows': 25,
        'skip_ports': ['prior_state']}

    timeseries = simulate_with_environment(compartment, settings)
    plot_simulation_output(timeseries, plot_settings, out_dir)