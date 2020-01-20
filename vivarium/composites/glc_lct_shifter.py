from __future__ import absolute_import, division, print_function

import os

from vivarium.composites.ecoli_master import compose_ecoli_master



def get_transport_config():
    """
    Convenience kinetics configuration for simplified glucose transport.
    This abstracts the PTS/GalP system to a single uptake kinetic
    with glc__D_e_external as the only cofactor.
    """
    transport_reactions = {
        'EX_glc__D_e': {
            'stoichiometry': {
                ('internal', 'g6p_c'): 1.0,
                ('external', 'glc__D_e'): -1.0,
                # ('internal', 'pep_c'): -1.0,  # TODO -- PEP needs mechanism for homeostasis to avoid depletion
                # ('internal', 'pyr_c'): 1.0
            },
            'is reversible': False,
            'catalyzed by': [('internal', 'PTSG')]},
        'EX_lac__D_e': {
            'stoichiometry': {
                ('external', 'lac__D_e'): -1.0,
                ('external', 'h_e'): -1.0,
                ('internal', 'lac__D_c'): 1.0,
                ('internal', 'h_c'): 1.0,
            },
            'is reversible': False,
            'catalyzed by': [('internal', 'LacY')]},
        }

    transport_kinetics = {
        'EX_glc__D_e': {
            ('internal', 'PTSG'): {
                ('external', 'glc__D_e'): 1e-1,
                ('internal', 'pep_c'): None,
                'kcat_f': -3e5}},
        'EX_lac__D_e': {   # TODO -- GLC inhibition?
            ('internal', 'LacY'): {
                ('external', 'lac__D_e'): 1e-1,
                ('external', 'h_e'): None,
                'kcat_f': -3e5}},
        }

    transport_initial_state = {
        'internal': {
            'PTSG': 1.8e-6,
            'g6p_c': 0.0,
            'pep_c': 1.8e-1,
            'pyr_c': 0.0,
            'LacY': 0.0,
            'lac__D_c': 0.0,
            'h_c': 100.0},
        'external': {
            'glc__D_e': 12.0,
            'lac__D_e': 0.0,
            'h_e': 100.0},
        'fluxes': {
            'EX_glc__D_e': 0.0,
            'EX_lac__D_e': 0.0}}  # TODO -- is this needed?

    transport_roles = {
        'internal': ['g6p_c', 'pep_c', 'pyr_c', 'h_c', 'PTSG', 'LacY'],
        'external': ['glc__D_e', 'lac__D_e', 'h_e']}

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}

def get_metabolism_config():
    metabolism_file = os.path.join('models', 'e_coli_core.json')
    regulation_logic = {'EX_lac__D_e': 'IF not (glc__D_e_external)'}
    return {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]},
        'model_path': metabolism_file,
        'regulation_logic': regulation_logic}

def get_expression_config():
    expression_rates = {'LacY': 1e-1}
    return {
        'expression_rates': expression_rates}

def get_degradation_config():
    degradation_rates = {'LacY': 1e-2}
    return {
        'degradation_rates': degradation_rates}

def get_default_config():


    # TODO -- in reconciler, make sure degradation never lowers state below 0


    return {
        'name': 'glc_lct_shifter',
        'transport': get_transport_config(),
        'metabolism': get_metabolism_config(),
        'expression': get_expression_config(),
        'degradation': get_degradation_config()
    }


def GlcLctShifter(config):
    """
    Load a glucose/lactose diauxic shift configuration into the kinetic_FBA composite
    """
    shifter_config = get_default_config()
    config.update(shifter_config)

    return compose_ecoli_master(config)



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, \
        simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'glc_lct_shifter')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(GlcLctShifter)

    # settings for simulation and plot
    options = GlcLctShifter({})['options']

    # define timeline
    timeline = [
        (0, {'environment': {
            'lac__D_e': 12.0}
        }),
        (500, {'environment': {
            'glc__D_e': 0.0}
        }),
        (1000, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'show_state': [
            ('environment', 'glc__D_e'),
            ('environment', 'lac__D_e'),
            ('environment', 'h_e'),
            ('reactions', 'GLCpts'),
            ('reactions', 'EX_glc__D_e'),
            ('reactions', 'EX_lac__D_e'),
            ('cell', 'g6p_c'),
            ('cell', 'pep_c'),
            ('cell', 'pyr_c'),
            ('cell', 'PTSG'),
            ('cell', 'lac__D_c'),
            ('cell', 'h_c'),
            ('cell', 'LacY'),
        ]}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
