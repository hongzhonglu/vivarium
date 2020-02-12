from __future__ import absolute_import, division, print_function

import os

# composite
from vivarium.composites.ode_expression import compose_ode_expression

# process configurations
from vivarium.processes.metabolism import get_e_coli_core_config
from vivarium.processes.convenience_kinetics import get_glc_lct_config
from vivarium.processes.ode_expression import get_lacy_config


# processes configurations
def get_metabolism_config():
    config = get_e_coli_core_config()

    # set flux bond tolerance for reactions in ode_expression's lacy_config
    metabolism_config = {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]}}

    config.update(metabolism_config)

    return config

def get_expression_config():
    # glc lct config from ode_expression
    config = get_lacy_config()

    # define regulation
    regulators = [('external', 'glc__D_e')]
    regulation = {'lacy_RNA': 'if not (external, glc__D_e) > 0.1'}
    reg_config = {
        'regulators': regulators,
        'regulation': regulation}

    config.update(reg_config)

    return config

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
        'expression': get_expression_config()}

    config.update(shifter_config)

    return compose_ode_expression(config)



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment
    from vivarium.actor.composition import simulate_with_environment, convert_to_timeseries, plot_simulation_output

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
        (400, {'environment': {
            'glc__D_e': 0.0}
        }),
        (1200, {'environment': {
            'glc__D_e': 12.0}
        }),
        (1800, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-12,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux'},
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
        'skip_roles': ['prior_state', 'null']}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
