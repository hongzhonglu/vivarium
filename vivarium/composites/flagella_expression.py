import os

import matplotlib.pyplot as plt
import numpy as np

from vivarium.data.nucleotides import nucleotides
from vivarium.data.amino_acids import amino_acids
from vivarium.data.chromosomes.flagella_chromosome import FlagellaChromosome
from vivarium.states.chromosome import Chromosome, rna_bases, sequence_monomers
from vivarium.processes.transcription import UNBOUND_RNAP_KEY
from vivarium.processes.translation import UNBOUND_RIBOSOME_KEY
from vivarium.composites.gene_expression import compose_gene_expression, plot_gene_expression_output



def degradation_sequences(sequence, promoters):
    return {
        promoter.last_terminator().product[0]: rna_bases(sequence_monomers(
            sequence,
            promoter.position,
            promoter.last_terminator().position))
        for promoter_key, promoter in promoters.items()}

def generate_flagella_compartment(config):
    flagella = FlagellaChromosome()
    plasmid = Chromosome(flagella.config)
    sequences = plasmid.product_sequences()

    print(sequences)

    molecules = {}
    for nucleotide in nucleotides.values():
        molecules[nucleotide] = 1000000
    for amino_acid in amino_acids.values():
        molecules[amino_acid] = 100000

    flagella_expression_config = {

        'transcription': {

            'sequence': flagella.config['sequence'],
            'templates': flagella.config['promoters'],
            'genes': flagella.config['genes'],
            'transcription_factors': flagella.transcription_factors,
            'promoter_affinities': flagella.promoter_affinities,
            'polymerase_occlusion': 30,
            'elongation_rate': 50},

        'translation': {

            'sequences': flagella.operon_sequences,
            'templates': flagella.transcript_templates,
            'concentration_keys': ['CRP', 'flhDC'],
            'transcript_affinities': flagella.transcript_affinities,

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
                        for transcript in flagella.config['genes'].keys()}}}},

        'initial_state': {
            'molecules': molecules,
            'transcripts': {
                gene: 0
                for gene in flagella.config['genes'].keys()},
            'proteins': {
                'CpxR': 10,
                'CRP': 10,
                'Fnr': 10,
                'endoRNAse': 1,
                UNBOUND_RIBOSOME_KEY: 10,
                UNBOUND_RNAP_KEY: 10}}}

    return compose_gene_expression(flagella_expression_config)


def plot_flagella_just_in_time(timeseries, config, out_dir='out'):

    name = config.get('name', 'gene_expression')
    ports = config.get('ports', {})
    transcripts = timeseries[ports['transcripts']]
    time = timeseries['time']

    # make timeseries heatmap
    transcript_ids = list(transcripts.keys())
    transcript_ids.reverse()  # reverse to get proper labeling with imshow
    n_transcripts = len(transcript_ids)
    transcripts_ts_map = [
        [value/ max(max(series), 1)
            for value in series]
            for transcript, series in transcripts.items()]

    # make figure
    fig = plt.figure(figsize=(5, 0.6*len(transcript_ids)))
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(transcripts_ts_map,
        extent=[time[0], time[-1], 0, n_transcripts],
        interpolation='nearest',
        aspect='auto',
        cmap='cividis'
        )
    ax.locator_params(axis='y', nbins=n_transcripts)

    # set y ticks locations and labels
    y_tick_locs = np.asarray([loc+0.5 for loc in range(n_transcripts)])
    ax.set_yticks(y_tick_locs)
    ax.set_yticklabels(transcript_ids)
    ax.set_xlabel('time (s)')

    # colorbar
    cbar = fig.colorbar(im)
    cbar.set_label('relative flourescence', rotation=270, labelpad=20)

    # save figure
    fig_path = os.path.join(out_dir, name)
    plt.savefig(fig_path, bbox_inches='tight')



if __name__ == '__main__':
    from vivarium.compartment.process import load_compartment, simulate_compartment
    from vivarium.compartment.composition import convert_to_timeseries

    out_dir = os.path.join('out', 'tests', 'flagella_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the compartment
    flagella_expression_compartment = load_compartment(generate_flagella_compartment)

    # run simulation
    settings = {
        'total_time': 1000}
    saved_state = simulate_compartment(flagella_expression_compartment, settings)
    del saved_state[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_state)

    plot_config = {
        'name': 'flagella_expression',
        'ports': {
            'transcripts': 'transcripts',
            'molecules': 'molecules',
            'proteins': 'proteins'}}

    plot_gene_expression_output(
        timeseries,
        plot_config,
        out_dir)

    # just-in-time figure
    plot_config2 = plot_config.copy()
    plot_config2.update({'name': 'flagella_just_in_time'})

    plot_flagella_just_in_time(
        timeseries,
        plot_config2,
        out_dir)
