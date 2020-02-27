from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state, load_compartment
from vivarium.compartment.composition import get_derivers, get_schema, simulate_with_environment, \
    convert_to_timeseries, plot_simulation_output

# processes
from vivarium.processes.division import Division, divide_condition
from vivarium.processes.metabolism import Metabolism, get_e_coli_core_config
from vivarium.processes.convenience_kinetics import ConvenienceKinetics
from vivarium.processes.ode_expression import ODE_expression



def compose_ode_expression(config):
    """
    A composite with kinetic transport, metabolism, and ode-based gene expression
    """

    ## Declare the processes.
    # Transport
    # load the kinetic parameters
    transport_config = config.get('transport', {})
    transport = ConvenienceKinetics(transport_config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids

    # Metabolism
    # get target fluxes from transport
    # load regulation function
    metabolism_config = config.get('metabolism', get_metabolism_config())
    metabolism_config.update({'constrained_reaction_ids': target_fluxes})
    metabolism = Metabolism(metabolism_config)

    # expression/degradation
    expression = ODE_expression(config.get('expression', {}))

    # Division
    # get initial volume from metabolism
    division_config = config.get('division', {})
    division_config.update({'initial_state': metabolism.initial_state})
    division = Division(division_config)

    # Place processes in layers
    processes = [
        {'transport': transport,
         'expression': expression},
        {'metabolism': metabolism},
        {'division': division}]

    # Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'null',
            'fluxes': 'flux',
            'global': 'global'},
        'metabolism': {
            'internal': 'cell',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux',
            'global': 'global'},
        'expression' : {
            'counts': 'cell_counts',
            'internal': 'cell',
            'external': 'environment'},
        'division': {
            'global': 'global'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])  # add deriver processes
    topology.update(derivers['deriver_topology'])  # add deriver topology

    # get schema
    schema = get_schema(processes, topology)

    # Initialize the states
    states = initialize_state(processes, topology, schema, config.get('initial_state', {}))

    options = {
        'name': config.get('name', 'master_composite'),
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'schema': schema,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition}

    return {
        'processes': processes,
        'states': states,
        'options': options}



# toy functions/ defaults
def get_metabolism_config():
    config = get_e_coli_core_config()

    # set flux bond tolerance for reactions in ode_expression's lacy_config
    metabolism_config = {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]}}

    config.update(metabolism_config)

    return config


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'ode_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_ode_expression, boot_config)

    # settings for simulation and plot
    options = compartment.configuration
    timeline = [(2520, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux'},
        'skip_roles': ['prior_state', 'null']}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    plot_simulation_output(timeseries, plot_settings, out_dir)
