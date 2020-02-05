from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state
from vivarium.environment.make_media import Media
from vivarium.utils.units import units

# processes
from vivarium.processes.deriver import Deriver
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.metabolism import Metabolism
from vivarium.processes.convenience_kinetics import ConvenienceKinetics
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.minimal_degradation import MinimalDegradation



# the composite function
def compose_master(config):
    """
    A composite with kinetic transport, metabolism, and gene expression
    """

    ## Declare the processes.
    # Transport
    # load the kinetic parameters
    transport_config = config.get('transport', default_transport_config())
    transport = ConvenienceKinetics(transport_config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids

    # Metabolism
    # get target fluxes from transport
    # load regulation function
    metabolism_config = config.get('metabolism', default_metabolism_config())
    metabolism_config.update({'constrained_reaction_ids': target_fluxes})
    metabolism = Metabolism(metabolism_config)

    # expression/degradation
    expression_config = config.get('expression', {})
    expression = MinimalExpression(expression_config)

    degradation_config = config.get('degradation', {})
    degradation = MinimalDegradation(degradation_config)

    # Division
    # get initial volume from metabolism
    division_config = config.get('division', {})
    division_config.update({'initial_state': metabolism.initial_state})
    division = Division(division_config)

    # Other processes
    deriver_config = config.get('deriver', {})
    deriver = Deriver(deriver_config)

    # Place processes in layers
    processes = [
        {'transport': transport,
         'expression': expression,
         'degradation': degradation},
        {'metabolism': metabolism},
        {'deriver': deriver,
         'division': division}
    ]

    # Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds'},
        'metabolism': {
            'internal': 'cell',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds'},
        'expression' : {
            'internal': 'cell_counts',  # updates counts, which the deriver converts to concentrations
            'external': 'environment'},
        'degradation': {
            'internal': 'cell_counts'},
        'division': {
            'internal': 'cell'},
        'deriver': {
            'counts': 'cell_counts',
            'state': 'cell',
            'prior_state': 'prior_state'},
    }

    # Initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': config.get('name', 'master_composite'),
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



# toy functions/ defaults
def default_metabolism_config():

    regulation = {
        'EX_lac__D_e': 'if not (external, glc__D_e) > 0.1'}

    metabolism_file = os.path.join('models', 'e_coli_core.json')

    # initial state
    # internal
    mass = 1339 * units.fg
    density = 1100 * units.g/units.L
    volume = mass.to('g') / density
    internal = {
            'mass': mass.magnitude,  # fg
            'volume': volume.to('fL').magnitude}

    # external
    # TODO -- generalize external to whatever BiGG model is loaded
    make_media = Media()
    external = make_media.get_saved_media('ecoli_core_GLC')
    initial_state = {
        'internal': internal,
        'external': external}

    return {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0]},
        'model_path': metabolism_file,
        'regulation': regulation,
        'initial_state': initial_state}

def default_transport_config():
    transport_reactions = {
        'EX_glc__D_e': {
            'stoichiometry': {
                ('internal', 'g6p_c'): 1.0,
                ('external', 'glc__D_e'): -1.0,
            },
            'is reversible': False,
            'catalyzed by': [('internal', 'PTSG')]}}

    transport_kinetics = {
        'EX_glc__D_e': {
            ('internal', 'PTSG'): {
                ('external', 'glc__D_e'): 1e-1,
                ('internal', 'pep_c'): None,
                'kcat_f': -3e5}}}

    transport_initial_state = {
        'internal': {
            'g6p_c': 0.0,
            'pep_c': 1.8e-1,
            'pyr_c': 0.0,
            'PTSG': 1.8e-6},
        'external': {
            'glc__D_e': 12.0},
        'fluxes': {
            'EX_glc__D_e': 0.0}}  # TODO -- is this needed?

    transport_roles = {
        'internal': ['g6p_c', 'pep_c', 'pyr_c', 'PTSG'],
        'external': ['glc__D_e']}

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, \
        simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'master_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_master, boot_config)

    # settings for simulation and plot
    options = compose_master({})['options']

    # define timeline
    timeline = [(1000, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'skip_roles': ['prior_state', 'null']
        }

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
