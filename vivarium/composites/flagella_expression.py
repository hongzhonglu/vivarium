import os

from vivarium.compartment.composition import load_compartment, simulate_compartment
from vivarium.data.nucleotides import nucleotides
from vivarium.data.amino_acids import amino_acids
from vivarium.data.chromosomes.flagella_chromosome import FlagellaChromosome
from vivarium.states.chromosome import Chromosome, Promoter, rna_bases, sequence_monomers
from vivarium.processes.transcription import UNBOUND_RNAP_KEY
from vivarium.processes.translation import generate_template, UNBOUND_RIBOSOME_KEY
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
            'concentration_keys': ['CRP'],
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
                'CRP': 10,
                'endoRNAse': 1,
                UNBOUND_RIBOSOME_KEY: 10,
                UNBOUND_RNAP_KEY: 10}}}

    return compose_gene_expression(flagella_expression_config)

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'flagella_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the compartment
    flagella_expression_compartment = load_compartment(generate_flagella_compartment)

    # run simulation
    settings = {
        'total_time': 800,
        'emit_timeseries': True,
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
