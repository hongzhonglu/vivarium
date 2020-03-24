import os
<<<<<<< HEAD

import matplotlib.pyplot as plt
import numpy as np
=======
>>>>>>> 9a604c059c91f9e7133e6157e046e7d877a6c15a

from vivarium.compartment.process import load_compartment, simulate_compartment
from vivarium.data.nucleotides import nucleotides
from vivarium.data.amino_acids import amino_acids
from vivarium.data.chromosomes.flagella_chromosome import FlagellaChromosome
from vivarium.states.chromosome import Chromosome, rna_bases, sequence_monomers
from vivarium.processes.transcription import UNBOUND_RNAP_KEY
from vivarium.processes.translation import UNBOUND_RIBOSOME_KEY
from vivarium.composites.gene_expression import compose_gene_expression, plot_gene_expression_output

<<<<<<< HEAD

=======
>>>>>>> 9a604c059c91f9e7133e6157e046e7d877a6c15a

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
            'concentration_keys': ['CRP', 'flhDC', 'fliA'],
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


def plot_timeseries_heatmaps(timeseries, config, out_dir='out'):
    ''' make a timeseries heatmap for each port specified in config['plot_ports'] '''

    name = config.get('name', 'timeseries')
    plot_ports = config.get('plot_ports', [])
    ports = config.get('ports', {})
    time = timeseries['time']

    # make timeseries heatmaps
    ts_heatmap = {}
    for port_id in plot_ports:
        port = timeseries[ports[port_id]]
        var_keys = list(port.keys())
        var_keys.reverse()  # reverse to get proper labeling with imshow
        var_series = [
            [value / max(max(series), 1)
                for value in series]
                for variable, series in port.items()]

        ts_heatmap[port_id] = {
            'keys': var_keys,
            'timeseries': var_series}

    # make figure for each port in plot_ports
    for port_id, heatmap in ts_heatmap.items():
        n_cols = 1
        n_vars = len(heatmap['keys'])

        fig = plt.figure(figsize=(4 * n_cols, 0.6 * n_vars))

        var_keys = heatmap['keys']
        var_series = heatmap['timeseries']
        n_vars = len(var_keys)
        ax = fig.add_subplot(111)

        im = ax.imshow(var_series,
            extent=[time[0], time[-1], 0, n_vars],
            interpolation='nearest',
            aspect='auto',
            cmap='cividis'
            )
        ax.locator_params(axis='y', nbins=n_vars)

        # set y ticks locations and labels
        y_tick_locs = np.asarray([loc+0.5 for loc in range(n_vars)])
        ax.set_yticks(y_tick_locs)
        ax.set_yticklabels(var_keys)
        ax.set_xlabel('time (s)')

        # colorbar
        cbar = fig.colorbar(im)
        cbar.set_label('relative flourescence', rotation=270, labelpad=20)

        # save figure
        figname = name + '_' + port_id
        fig_path = os.path.join(out_dir, figname)
        plt.savefig(fig_path, bbox_inches='tight')



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'flagella_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the compartment
    flagella_expression_compartment = load_compartment(generate_flagella_compartment)

    # run simulation
    settings = {
        'total_time': 2400,
    }
    timeseries = simulate_compartment(flagella_expression_compartment, settings)

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
    plot_config2.update({
        'name': 'flagella',
        'plot_ports': ['transcripts', 'proteins']})

    plot_timeseries_heatmaps(
        timeseries,
        plot_config2,
        out_dir)
