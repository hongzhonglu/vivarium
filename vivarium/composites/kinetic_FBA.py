from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.derive_volume import DeriveVolume
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.BiGG_metabolism import BiGGMetabolism
from vivarium.processes.convenience_kinetics import ConvenienceKinetics


# the composite function
def compose_kinetic_FBA(config):
    """
    A composite with kinetic transport, metabolism, and regulation
    TODO (eran) -- fit glc/lct uptake rates to growth rates
    """

    ## Declare the processes.
    # Transport
    # load the kinetic parameters
    transport_config = config.get('transport', default_transport_config())
    transport = ConvenienceKinetics(transport_config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids

    # Metabolism
    # get target fluxes from transport
    # load regulation function
    metabolism_config = config.get('metabolism', default_metabolism_config())
    metabolism_config.update({'constrained_reaction_ids': target_fluxes})
    metabolism = BiGGMetabolism(metabolism_config)

    # Division
    # get initial volume from metabolism
    division_config = config.get('division', {})
    division_config.update({'initial_state': metabolism.initial_state})
    division = Division(division_config)

    # Other processes
    deriver_config = config.get('volume', {})
    deriver = DeriveVolume(deriver_config)

    # Place processes in layers
    processes = [
        {'transport': transport},
        {'metabolism': metabolism},
        {'deriver': deriver,
         'division': division}
    ]

    # Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds'},
        'metabolism': {
            'internal': 'cell',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds'},
        'division': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
    }

    # Initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'kinetic_FBA_composite',
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}



def default_metabolism_config():
    def regulation(state):
        regulation_logic = {
            'EX_lac__D_e': bool(not state[('external', 'glc__D_e')] > 0.1),
        }
        return regulation_logic

    metabolism_file = os.path.join('models', 'e_coli_core.json')

    return {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0]},
        'model_path': metabolism_file,
        'regulation': regulation}

def default_transport_config():
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
            'catalyzed by': [('internal', 'PTSG')]}}

    transport_kinetics = {
        'EX_glc__D_e': {
            ('internal', 'PTSG'): {
                ('external', 'glc__D_e'): 1e-1,
                ('internal', 'pep_c'): None,
                'kcat_f': -3e5}}}

    transport_initial_state = {
        'internal': {
            'g6p_c': 0.0,
            'pep_c': 1.8e-1,
            'pyr_c': 0.0,
            'PTSG': 1.8e-6},
        'external': {
            'glc__D_e': 12.0},
        'fluxes': {
            'EX_glc__D_e': 0.0}}  # TODO -- is this needed?

    transport_roles = {
        'internal': ['g6p_c', 'pep_c', 'pyr_c', 'PTSG'],
        'external': ['glc__D_e']}

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}



if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, \
        simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'kinetic_FBA_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_kinetic_FBA)

    # settings for simulation and plot
    options = compose_kinetic_FBA({})['options']

    # define timeline
    timeline = [
        (0, {'environment': {
            'lac__D_e': 12.0}
        }),
        (1000, {'environment': {
            'glc__D_e': 0.0}
        }),
        (3500, {})]

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
            ('reactions', 'GLCpts'),
            ('reactions', 'EX_glc__D_e'),
            ('reactions', 'EX_lac__D_e'),
            ('cell', 'g6p_c'),
            ('cell', 'pep_c'),
            ('cell', 'pyr_c'),
            ('cell', 'PTSG'),
        ]}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove the first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
