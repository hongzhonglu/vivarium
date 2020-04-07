import os
from vivarium.utils.units import units

from vivarium.compartment.composition import load_compartment, simulate_compartment
from vivarium.data.proteins import GFP
from vivarium.data.chromosomes.gfp_chromosome import gfp_plasmid_config
from vivarium.states.chromosome import Chromosome, Promoter, rna_bases, sequence_monomers
from vivarium.utils.polymerize import generate_template
from vivarium.composites.gene_expression import compose_gene_expression, plot_gene_expression_output
from vivarium.processes.transcription import UNBOUND_RNAP_KEY
from vivarium.processes.translation import UNBOUND_RIBOSOME_KEY
from vivarium.environment.make_media import Media

from vivarium.processes.derive_globals import AVOGADRO

def degradation_sequences(sequence, promoters):
    return {
        promoter.last_terminator().product[0]: rna_bases(sequence_monomers(
            sequence,
            promoter.position,
            promoter.last_terminator().position))
        for promoter_key, promoter in promoters.items()}

def generate_gfp_compartment(config):
    media = Media()
    PURE = {
        key: value * units.mmol / units.L
        for key, value in media.get_saved_media('PURE_Fuji_2014').items()}

    # TODO: deal with volume
    volume = 1e-15 * units.L
    mmol_to_counts = AVOGADRO.to('1/mmol') * volume.to('L')
    
    print(mmol_to_counts)

    PURE_counts = {
        key: int(value * mmol_to_counts)
        for key, value in PURE.items()}

    print(PURE)
    print(PURE_counts)

    plasmid = Chromosome(gfp_plasmid_config)
    sequences = plasmid.product_sequences()

    print(sequences)

    proteins = {
        UNBOUND_RNAP_KEY: PURE_counts[UNBOUND_RNAP_KEY],
        UNBOUND_RIBOSOME_KEY: PURE_counts[UNBOUND_RIBOSOME_KEY],
        'GFP': 0,
        'endoRNAse': 1}

    gfp_config = {

        'transcription': {

            'sequence': gfp_plasmid_config['sequence'],
            'templates': gfp_plasmid_config['promoters'],
            'genes': gfp_plasmid_config['genes'],
            'promoter_affinities': {
                ('T7',): 0.001},

            'polymerase_occlusion': 30,
            'elongation_rate': 50},

        'translation': {

            'sequences': {
                'GFP_RNA': GFP.sequence},
            'templates': {
                'GFP_RNA': generate_template(
                    'GFP_RNA', len(GFP.sequence), ['GFP'])},
            'transcript_affinities': {
                'GFP_RNA': 0.001},

            'elongation_rate': 22,
            'polymerase_occlusion': 50},

        'degradation': {
            
            'sequences': sequences,
            'catalysis_rates': { # TODO: provide kcat for each RNA variety
                'endoRNAse': 0},
            'degradation_rates': {
                'transcripts': {
                    'endoRNAse': {
                        'GFP_RNA': 1e-23}}}},

        'initial_state': {
            'molecules': PURE_counts,
            'transcripts': {'GFP_RNA': 0},
            'proteins': proteins}}

    return compose_gene_expression(gfp_config)

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'gfp_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the compartment
    gfp_expression_compartment = load_compartment(generate_gfp_compartment)

    # run simulation
    settings = {
        'total_time': 100,
    }
    timeseries = simulate_compartment(gfp_expression_compartment, settings)

    plot_config = {
        'name': 'gfp_expression',
        'ports': {
            'transcripts': 'transcripts',
            'molecules': 'molecules',
            'proteins': 'proteins'}}

    plot_gene_expression_output(
        timeseries,
        plot_config,
        out_dir)
