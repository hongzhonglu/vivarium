from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import initialize_state
from vivarium.compartment.composition import get_derivers, get_schema

# processes
from vivarium.processes.division import Division, divide_condition
from vivarium.processes.metabolism import Metabolism, get_e_coli_core_config
from vivarium.processes.convenience_kinetics import ConvenienceKinetics
from vivarium.processes.transcription import Transcription
from vivarium.processes.translation import Translation
from vivarium.processes.degradation import RnaDegradation



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

    # expression
    transcription_config = config.get('transcription', {})
    translation_config = config.get('translation', {})
    degradation_config = config.get('degradation', {})
    transcription = Transcription(transcription_config)
    translation = Translation(translation_config)
    degradation = RnaDegradation(degradation_config)

    # Division
    # get initial volume from metabolism
    division_config = config.get('division', {})
    division_config.update({'initial_state': metabolism.initial_state})
    division = Division(division_config)

    # Place processes in layers
    processes = [
        {'transport': transport,
         'transcription': transcription,
         'translation': translation,
         'degradation': degradation},
        {'metabolism': metabolism},
        {'division': division}]

    # Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds',
            'global': 'global'},
        'metabolism': {
            'internal': 'cell',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds',
            'global': 'global'},
        'transcription': {
            'chromosome': 'chromosome',
            'molecules': 'cell',
            'transcripts': 'transcripts'},
        'translation': {
            'ribosomes': 'ribosomes',
            'molecules': 'cell',
            'transcripts': 'transcripts',
            'proteins': 'proteins'},
        'degradation': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'cell',
            'global': 'global'},
        'division': {
            'global': 'global'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])  # add deriver processes
    topology.update(derivers['deriver_topology'])  # add deriver topology

    # get schema
    schema = get_schema(processes, topology)

    # initialize the states
    states = initialize_state(processes, topology, schema, config.get('initial_state', {}))

    options = {
        'name': config.get('name', 'master_composite'),
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'topology': topology,
        'schema': schema,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition}

    return {
        'processes': processes,
        'states': states,
        'options': options}



# toy functions/ defaults
def default_metabolism_config():
    config = get_e_coli_core_config()

    # set flux bond tolerance for reactions in ode_expression's lacy_config
    metabolism_config = {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]}}

    config.update(metabolism_config)

    return config

def default_transport_config():
    return {}



if __name__ == '__main__':
    from vivarium.compartment.process import load_compartment
    from vivarium.compartment.composition import simulate_with_environment, convert_to_timeseries, plot_simulation_output
    from vivarium.composites.gene_expression import plot_gene_expression_output

    out_dir = os.path.join('out', 'tests', 'master_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {}  # {'emitter': 'null'}
    compartment = load_compartment(compose_master, boot_config)

    # settings for simulation and plot
    options = compartment.configuration

    # define timeline
    timeline = [(2520, {})] # 2520 sec (42 min) is the expected doubling time in minimal media

    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'skip_ports': ['prior_state', 'null']}

    expression_plot_settings = {
        'name': 'gene_expression',
        'roles': {
            'transcripts': 'transcripts',
            'molecules': 'cell',
            'proteins': 'proteins'}}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    plot_gene_expression_output(timeseries, expression_plot_settings, out_dir)
    plot_simulation_output(timeseries, plot_settings, out_dir)
