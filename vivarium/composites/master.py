from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import (
    initialize_state)
from vivarium.compartment.composition import (
    get_derivers,
    simulate_with_environment,
    plot_simulation_output, load_compartment)
from vivarium.composites.gene_expression import plot_gene_expression_output

# processes
from vivarium.processes.division import Division, divide_condition
from vivarium.processes.metabolism import Metabolism, get_iAF1260b_config
from vivarium.processes.convenience_kinetics import ConvenienceKinetics, get_glc_lct_config
from vivarium.processes.transcription import Transcription
from vivarium.processes.translation import Translation
from vivarium.processes.degradation import RnaDegradation
from vivarium.processes.complexation import Complexation



# the composite function
def compose_master(config):
    """
    A composite with kinetic transport, metabolism, and gene expression
    """

    t_m_config = default_transport_metabolism_config()

    ## Declare the processes.
    # Transport
    # load the kinetic parameters
    transport_config = config.get('transport', t_m_config['transport'])
    transport = ConvenienceKinetics(transport_config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids

    # Metabolism
    # get target fluxes from transport
    # load regulation function
    metabolism_config = config.get('metabolism', t_m_config['metabolism'])
    metabolism_config.update({'constrained_reaction_ids': target_fluxes})
    metabolism = Metabolism(metabolism_config)

    # expression
    transcription_config = config.get('transcription', {})
    translation_config = config.get('translation', {})
    degradation_config = config.get('degradation', {})
    transcription = Transcription(transcription_config)
    translation = Translation(translation_config)
    degradation = RnaDegradation(degradation_config)
    complexation = Complexation(config.get('complexation', {}))

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
         'degradation': degradation,
         'complexation': complexation},
        {'metabolism': metabolism},
        {'division': division}]

    # Make the topology
    # for each process, map process ports to store ids
    topology = {
        'transport': {
            'internal': 'metabolites',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds',
            'global': 'global'},

        'metabolism': {
            'internal': 'metabolites',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds',
            'global': 'global'},

        'transcription': {
            'chromosome': 'chromosome',
            'molecules': 'metabolites',
            'proteins': 'proteins',
            'transcripts': 'transcripts',
            'factors': 'concentrations'},

        'translation': {
            'ribosomes': 'ribosomes',
            'molecules': 'metabolites',
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'concentrations': 'concentrations'},

        'degradation': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'metabolites',
            'global': 'global'},

        'complexation': {
            'monomers': 'proteins',
            'complexes': 'proteins'},

        'division': {
            'global': 'global'}}

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
        'name': config.get('name', 'master_composite'),
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition}

    return {
        'processes': processes,
        'derivers': deriver_processes,
        'states': states,
        'options': options}



# toy functions/ defaults
def default_transport_metabolism_config():
    transport_config = get_glc_lct_config()
    metabolism_config = get_iAF1260b_config()

    # set flux bond tolerance for reactions in ode_expression's lacy_config
    metabolism_config.update({
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lcts_e': [1.05, 1.0]}})

    return {
        'transport': transport_config,
        'metabolism': metabolism_config
    }




if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'master_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_master)

    # settings for simulation and plot
    options = compartment.configuration

    # define timeline
    timeline = [(2520, {})] # 2520 sec (42 min) is the expected doubling time in minimal media

    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline,
    }

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'skip_ports': ['prior_state', 'null']}

    expression_plot_settings = {
        'name': 'gene_expression',
        'ports': {
            'transcripts': 'transcripts',
            'molecules': 'metabolites',
            'proteins': 'proteins'}}

    # saved_state = simulate_compartment(compartment, settings)
    timeseries = simulate_with_environment(compartment, settings)
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    plot_gene_expression_output(timeseries, expression_plot_settings, out_dir)
    plot_simulation_output(timeseries, plot_settings, out_dir)
