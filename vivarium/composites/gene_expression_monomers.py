from __future__ import absolute_import, division, print_function

import os
import argparse

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import networkx as nx

from vivarium.compartment.process import initialize_state
from vivarium.compartment.composition import (
    get_derivers,
    load_compartment,
    simulate_compartment,
    plot_compartment_topology
)
from vivarium.utils.make_network import save_network

# chromosomes
from vivarium.data.chromosomes.lac_chromosome import LacChromosome

# processes
from vivarium.processes.transcription import Transcription, UNBOUND_RNAP_KEY
from vivarium.processes.translation import Translation, UNBOUND_RIBOSOME_KEY
from vivarium.processes.degradation import RnaDegradation
from vivarium.processes.complexation import Complexation
from vivarium.processes.division import Division, divide_condition
from vivarium.data.amino_acids import amino_acids
from vivarium.data.nucleotides import nucleotides


def compose_gene_expression_monomers(config):
    # declare the processes
    transcription = Transcription(config.get('transcription', {}))
    translation = Translation(config.get('translation', {}))
    degradation = RnaDegradation(config.get('degradation', {}))
    division = Division(config)

    # place processes in layers
    processes = {
        'transcription': transcription,
        'translation': translation,
        'degradation': degradation,
        'complexation': complexation,
        'division': division}

    # make the topology
    topology = {
        'transcription': {
            'chromosome': 'chromosome',
            'molecules': 'molecules',
            'proteins': 'proteins',
            'transcripts': 'transcripts',
            'factors': 'concentrations'},

        'translation': {
            'ribosomes': 'ribosomes',
            'molecules': 'molecules',
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'concentrations': 'concentrations',
            'global': 'global'},

        'degradation': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'molecules',
            'global': 'global'},

        'division': {
            'global': 'global'}}

    # add derivers
    deriver_config = {
        'mass': {
            'dark_mass': 1339,  # fg
            'ports': {'global': 'global'}}}
    derivers = get_derivers(processes, topology, deriver_config)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes.copy()
    all_processes.update(derivers['deriver_processes'])
    topology.update(derivers['deriver_topology'])  # add derivers to the topology

    # initialize the states
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    options = {
        'name': 'gene_expression_composite',
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


# configurations
def get_lac_operon_config(config):
    lac_data = LacChromosome(config)
    chromosome_config = lac_data.chromosome_config
    sequences = lac_data.chromosome.product_sequences()

    molecules = {}
    for nucleotide in nucleotides.values():
        molecules[nucleotide] = 5000000
    for amino_acid in amino_acids.values():
        molecules[amino_acid] = 1000000

    molecules['GLC-6-P'] = 10000  # TODO -- make sure this is read by transcription

    lac_config = {

        'transcription': {

            'sequence': chromosome_config['sequence'],
            'templates': chromosome_config['promoters'],
            'genes': chromosome_config['genes'],
            'transcription_factors': lac_data.transcription_factors,
        # TODO -- split keys between protein_TFs, metabolite_TFs
            'promoter_affinities': lac_data.promoter_affinities,
            'polymerase_occlusion': 30,
            'elongation_rate': 50},

        'translation': {

            'sequences': lac_data.protein_sequences,
            'templates': lac_data.transcript_templates,
            # 'concentration_keys': [],
            'transcript_affinities': lac_data.transcript_affinities,
            'elongation_rate': 22,
            'polymerase_occlusion': 50},

        'degradation': {

            'sequences': sequences,
            'catalysis_rates': {
                'endoRNAse': 0.01},
            'degradation_rates': {
                'transcripts': {
                    'endoRNAse': {
                        transcript: 1e-23
                        for transcript in chromosome_config['genes'].keys()}}}},

        'initial_state': {
            'molecules': molecules,
            'transcripts': {
                gene: 0
                for gene in chromosome_config['genes'].keys()},
            'proteins': {
                'lacY': 0,
                'lacZ': 0,
                'lacA': 0,
                'endoRNAse': 1,
                UNBOUND_RIBOSOME_KEY: 200,  # e. coli has ~ 20000 ribosomes
                UNBOUND_RNAP_KEY: 200}}}

    return lac_config


def generate_lac_compartment(config):
    lac_operon_config = get_lac_operon_config(config)
    return compose_gene_expression_monomers(lac_operon_config)


def run_lac_operon(config={}, out_dir='out'):
    lac_config = get_lac_operon_config(config)
    compartment = load_compartment(generate_lac_compartment, lac_config)

    # plot topology
    settings = {'show_ports': True}
    plot_compartment_topology(
        compartment,
        settings,
        out_dir,
        'lac_operon_topology')

    # run simulation
    settings = {
        'total_time': 100,  # 2400
        'verbose': True}
    timeseries = simulate_compartment(compartment, settings)

    plot_config = {
        'name': 'lac_operon_expression',
        'ports': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'molecules'}}

    plot_gene_expression_output(
        timeseries,
        plot_config,
        out_dir)




if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'gene_expression_monomers_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='gene_expression_monomers')
    parser.add_argument('--lac', '-l', action='store_true', default=False)
    args = parser.parse_args()

    if args.lac:
        run_lac_operon({}, out_dir)

    else:
        # load the compartment
        gene_expression_compartment = load_compartment(compose_gene_expression)

        # run simulation
        sim_settings = {
            'total_time': 100,
        }
        timeseries = simulate_compartment(gene_expression_compartment, sim_settings)

        plot_settings = {
            'name': 'gene_expression_monomers',
            'ports': {
                'transcripts': 'transcripts',
                'molecules': 'molecules',
                'proteins': 'proteins'}}

        plot_gene_expression_output(timeseries, plot_settings, out_dir)
