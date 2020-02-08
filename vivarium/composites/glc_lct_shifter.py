from __future__ import absolute_import, division, print_function

import os

from vivarium.environment.make_media import Media
from vivarium.processes.convenience_kinetics import get_glc_lct_config
from vivarium.utils.units import units
from vivarium.composites.ode_expression import compose_ode_expression


# processes configurations
def get_metabolism_config():

    metabolism_file = os.path.join('models', 'e_coli_core.json')

    # initial state
    # internal
    mass = 1339 * units.fg
    density = 1100 * units.g / units.L
    volume = mass.to('g') / density
    internal = {
        'mass': mass.magnitude,  # fg
        'volume': volume.to('fL').magnitude}

    # external
    make_media = Media()
    external = make_media.get_saved_media('ecoli_core_GLC')  # TODO -- generalize external to whatever BiGG model is loaded
    initial_state = {
        'internal': internal,
        'external': external}

    return {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]},
        'model_path': metabolism_file,
        'initial_state': initial_state}

def get_expression_config():
    # define regulation
    regulators = [('external', 'glc__D_e')]
    regulation = {'lacy_RNA': 'if not (external, glc__D_e) > 0.1'}

    transcription_rates = {'lacy_RNA': 1e-20}
    translation_rates = {'LacY': 1e-18}
    protein_map = {'LacY': 'lacy_RNA'}
    degradation_rates = {
        'lacy_RNA': 1e-15,
        'LacY': 1e-30}
    initial_state = {}

    molecules = list(transcription_rates.keys()) + list(translation_rates.keys())

    return {
        'regulators': regulators,
        'regulation': regulation,
        'transcription_rates': transcription_rates,
        'translation_rates': translation_rates,
        'degradation_rates': degradation_rates,
        'protein_map': protein_map,
        'initial_state': initial_state,
        'counted_molecules': molecules}

# composite configuration
def compose_glc_lct_shifter(config):
    """
    Load a glucose/lactose diauxic shift configuration into the kinetic_FBA composite
    TODO (eran) -- fit glc/lct uptake rates to growth rates
    """
    shifter_config = {
        'name': 'glc_lct_shifter',
        'transport': get_glc_lct_config(),
        'metabolism': get_metabolism_config(),
        'expression': get_expression_config(),
        'deriver': get_expression_config()
    }
    config.update(shifter_config)

    return compose_ode_expression(config)



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, \
        simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'glc_lct_shifter')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_glc_lct_shifter, boot_config)

    # settings for simulation and plot
    options = compartment.configuration

    # define timeline
    timeline = [
        (0, {'environment': {
            'lac__D_e': 12.0}
        }),
        (200, {'environment': {
            'glc__D_e': 0.0}
        }),
        (800, {'environment': {
            'glc__D_e': 12.0}
        }),
        (1200, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-12,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'show_state': [
            ('environment', 'glc__D_e'),
            ('environment', 'lac__D_e'),
            ('reactions', 'GLCpts'),
            ('reactions', 'EX_glc__D_e'),
            ('reactions', 'EX_lac__D_e'),
            ('cell', 'g6p_c'),
            ('cell', 'PTSG'),
            ('cell', 'lac__D_c'),
            ('cell', 'lacy_RNA'),
            ('cell', 'LacY')],
        'skip_roles': ['prior_state', 'null']
    }

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
