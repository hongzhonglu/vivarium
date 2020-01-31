import os
from vivarium.utils.units import units

from scipy import constants

from vivarium.data.amino_acids import amino_acid_records
from vivarium.data.proteins import GFP
from vivarium.data.chromosome import gfp_plasmid_config
from vivarium.processes.translation import generate_template
from vivarium.composites.gene_expression import compose_gene_expression, plot_gene_expression_output
from vivarium.environment.make_media import Media

def generate_gfp_compartment(config):
    media = Media()
    PURE = {
        key: value * units.mmol / units.L
        for key, value in media.get_saved_media('PURE_Fuji_2014').items()}

    # TODO: deal with volume
    volume = 1e-15 * units.L
    avogadro = constants.N_A * 1 / units.mol
    mmol_to_count = avogadro.to('1/mmol') * volume.to('L')
    
    print(mmol_to_count)

    PURE_counts = {
        key: int(value * mmol_to_count)
        for key, value in PURE.items()}

    print(PURE)
    print(PURE_counts)

    gfp_config = {

        'transcription': {

            'sequence': gfp_plasmid_config['sequence'],
            'templates': gfp_plasmid_config['promoters'],
            'genes': gfp_plasmid_config['genes'],
            'promoter_affinities': {
                'T7': 1.0},

            'advancement_rate': 1.0,
            'elongation_rate': 50, # TODO: find real elongation rate
            'monomer_ids': ['A', 'T', 'G', 'C'],
            'transcript_ids': ['GFP_RNA'],
            'promoter_order': ['T7']},

        'translation': {

            'sequences': {
                'GFP_RNA': GFP.protein_sequence},
            'templates': {
                'GFP_RNA': generate_template(
                    'GFP_RNA', len(GFP.sequence), 'GFP')},
            'transcript_affinities': {
                'GFP_RNA': 1.0},

            'elongation_rate': 22,
            'advancement_rate': 1.0,
            'monomer_ids': [
                record['abbreviation'] for record in amino_acid_records],
            'protein_ids': ['GFP'],
            'transcript_order': ['GFP_RNA']}}

    return compose_gene_expression(gfp_config)

if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, simulate_compartment, convert_to_timeseries

    out_dir = os.path.join('out', 'tests', 'gfp_expression_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the compartment
    gfp_expression_compartment = load_compartment(generate_gfp_compartment)

    # run simulation
    settings = {
        'total_time': 1}
    saved_state = simulate_compartment(gfp_expression_compartment, settings)
    del saved_state[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_state)
    plot_gene_expression_output(timeseries, out_dir)
    
