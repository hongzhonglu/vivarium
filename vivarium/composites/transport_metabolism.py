from __future__ import absolute_import, division, print_function

import os

from vivarium.environment.make_media import Media
from vivarium.utils.units import units
from vivarium.composites.master import compose_master


# processes configurations
def get_transport_config():
    """
    Convenience kinetics configuration for simplified glucose/lactose transport.
    This abstracts the PTS/GalP system to a single uptake kinetic
    with glc__D_e_external as the only cofactor.
    """
    transport_reactions = {
        'EX_glc__D_e': {
            'stoichiometry': {
                ('internal', 'g6p_c'): 1.0,
                ('external', 'glc__D_e'): -1.0,
                ('internal', 'pep_c'): -1.0,
                ('internal', 'pyr_c'): 1.0
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
                'kcat_f': -1e5}},
        }

    transport_initial_state = {
        'internal': {
            'PTSG': 1.8e-6,  # concentration (mmol/L)
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

    metabolism_file = os.path.join('models', 'iAF1260b.json')

    # initial state
    mass = 1339 * units.fg
    density = 1100 * units.g / units.L
    volume = mass.to('g') / density
    internal = {
        'mass': mass.magnitude,  # fg
        'volume': volume.to('fL').magnitude}
    initial_state = {'internal': internal}

    return {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lac__D_e': [1.05, 1.0]},
        'model_path': metabolism_file,
        'initial_state': initial_state}


# composite configuration
def compose_transport_metabolism(config):

    composite_config = {
        'name': 'transport_metabolism',
        # 'transport': get_transport_config(),
        'metabolism': get_metabolism_config()
    }
    config.update(composite_config)

    return compose_master(config)



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, \
        simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'transport_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_transport_metabolism, boot_config)

    # settings for simulation and plot
    options = compose_transport_metabolism({})['options']

    # define timeline
    timeline = [(100, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'skip_roles': ['prior_state', 'null']
    }

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
