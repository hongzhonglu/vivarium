from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.growth import Growth
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.convenience_kinetics import ConvenienceKinetics, get_glc_lct_config
from vivarium.processes.derive_globals import DeriveGlobals



def compose_growth_division(config):

    # declare the processes
    transport_config = get_glc_lct_config()
    transport = ConvenienceKinetics(transport_config)
    growth = Growth(config)
    division = Division(config)
    expression = MinimalExpression(config)
    deriver = DeriveGlobals(config)

    # place processes in layers
    processes = [
        {'transport': transport,
         'growth': growth,
         'expression': expression},
        {'deriver': deriver,
         'division': division}]

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
            'external': 'environment'},
        'deriver': {
            'counts': 'cell_counts',
            'state': 'cell',
            'global': 'global'}}

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'growth_division_composite',
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
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
        'total_time': 20}

    plot_settings = {
        'max_rows': 25,
        'skip_roles': ['prior_state']}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)