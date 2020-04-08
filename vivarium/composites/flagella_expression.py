import os

import matplotlib.pyplot as plt
import numpy as np

from vivarium.compartment.composition import load_compartment, simulate_compartment
from vivarium.data.nucleotides import nucleotides
from vivarium.data.amino_acids import amino_acids
from vivarium.data.chromosomes.flagella_chromosome import flagella_chromosome
from vivarium.states.chromosome import Chromosome, rna_bases, sequence_monomers
from vivarium.processes.transcription import UNBOUND_RNAP_KEY
from vivarium.processes.translation import UNBOUND_RIBOSOME_KEY
from vivarium.composites.gene_expression import compose_gene_expression, plot_gene_expression_output

def get_flg_expression_config():
    plasmid = Chromosome(flagella_chromosome.config)
    sequences = plasmid.product_sequences()

    molecules = {}
    for nucleotide in nucleotides.values():
        molecules[nucleotide] = 1000000
    for amino_acid in amino_acids.values():
        molecules[amino_acid] = 100000

    flagella_expression_config = {

        'transcription': {

            'sequence': flagella_chromosome.config['sequence'],
            'templates': flagella_chromosome.config['promoters'],
            'genes': flagella_chromosome.config['genes'],
            'transcription_factors': flagella_chromosome.transcription_factors,
            'promoter_affinities': flagella_chromosome.promoter_affinities,
            'polymerase_occlusion': 30,
            'elongation_rate': 50},

        'translation': {

            'sequences': flagella_chromosome.protein_sequences,
            'templates': flagella_chromosome.transcript_templates,
            'concentration_keys': ['CRP', 'flhDC', 'fliA'],
            'transcript_affinities': flagella_chromosome.transcript_affinities,

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
                        for transcript in flagella_chromosome.config['genes'].keys()}}}},

        'complexation': {
            'monomer_ids': flagella_chromosome.complexation_monomer_ids,
            'complex_ids': flagella_chromosome.complexation_complex_ids,
            'stoichiometry': flagella_chromosome.complexation_stoichiometry,
            'rates': flagella_chromosome.complexation_rates},

        'initial_state': {
            'molecules': molecules,
            'transcripts': {
                gene: 0
                for gene in flagella_chromosome.config['genes'].keys()},
            'proteins': {
                'CpxR': 10,
                'CRP': 10,
                'Fnr': 10,
                'endoRNAse': 1,
                UNBOUND_RIBOSOME_KEY: 10,
                UNBOUND_RNAP_KEY: 10}}}

    return flagella_expression_config


def generate_flagella_compartment(config):
    flagella_expression_config = get_flg_expression_config()

    return compose_gene_expression(flagella_expression_config)


def plot_timeseries_heatmaps(timeseries, config, out_dir='out'):
    ''' make a timeseries heatmap for each port specified in config['plot_ports'] '''

    name = config.get('name', 'timeseries')
    plot_ports = config.get('plot_ports', {})
    ports = config.get('ports', {})
    time = timeseries['time']

    def relative_to_max(series):
        relative = max(max(series), 1)
        return [
            value / relative
            for value in series]

    # make timeseries heatmaps
    ts_heatmap = {}
    for port_id, order in plot_ports.items():
        port = timeseries[ports[port_id]]
        var_keys = list(order)

        var_series = [
            relative_to_max(port[key])
            for key in var_keys]

        var_keys.reverse()  # reverse to get proper labeling with imshow

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
        'verbose': True}
    timeseries = simulate_compartment(flagella_expression_compartment, settings)

    plot_config = {
        'name': 'flagella_expression',
        'ports': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'molecules'}}

    plot_gene_expression_output(
        timeseries,
        plot_config,
        out_dir)

    # just-in-time figure
    plot_config2 = plot_config.copy()
    plot_config2.update({
        'name': 'flagella',
        'plot_ports': {
            'transcripts': list(flagella_chromosome.config['genes'].keys()),
            'proteins': flagella_chromosome.complexation_monomer_ids + flagella_chromosome.complexation_complex_ids,
            'molecules': list(nucleotides.values()) + list(amino_acids.values())}})

    plot_timeseries_heatmaps(
        timeseries,
        plot_config2,
        out_dir)
