from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state
from vivarium.actor.composition import get_derivers, get_schema

# processes
from vivarium.processes.growth import Growth
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.convenience_kinetics import ConvenienceKinetics, get_glc_lct_config



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
    # for each process, map process roles to compartment roles
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

    # get schema
    schema = get_schema(processes, topology)

    # initialize the states
    states = initialize_state(processes, topology, schema, config.get('initial_state', {}))

    options = {
        'name': 'growth_division_composite',
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'schema': schema,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}


if __name__ == '__main__':
    from vivarium.actor.process import load_compartment
    from vivarium.actor.composition import simulate_with_environment, convert_to_timeseries, plot_simulation_output

    out_dir = os.path.join('out', 'tests', 'growth_division_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_growth_division)

    # settings for simulation and plot
    options = compartment.configuration
    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 100}

    plot_settings = {
        'max_rows': 25,
        'skip_roles': ['prior_state']}

    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)