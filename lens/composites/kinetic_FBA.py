from __future__ import absolute_import, division, print_function

import copy
import os

from lens.actor.process import initialize_state

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.division import Division, divide_condition, divide_state
from lens.processes.BiGG_metabolism import BiGGMetabolism
from lens.processes.convenience_kinetics import ConvenienceKinetics


# BiGG model for metabolism
METABOLISM_FILE = os.path.join('models', 'e_coli_core.json')

# convenience kinetics configuration for transport
def get_transport_config():
    # stoichiometry needs to match metabolism -- perhaps this can be extracted?
    transport_reactions = {
        'GLCpts': {
            'stoichiometry': {
                'g6p_c': 1.0,
                'glc__D_e': -1.0,
                'pep_c': -1.0,
                'pyr_c': 1.0},
            'is reversible': False,
            'catalyzed by': ['PTSG']}}

    # very simplified PTS
    transport_kinetics = {
        'GLCpts': {
            'PTSG': {
                'glc__D_e': 3.0,
                'pep_c': 2.0,
                'kcat_f': 1.0}}}

    transport_initial_state = {
        'internal': {
            'g6p_c': 1.0,
            'pep_c': 1.0,
            'pyr_c': 1.0,
            'PTSG': 1.0},
        'external': {
            'glc__D_e': 12.0},
        'fluxes': {
            'GLCpts': 1.0}}  # TODO -- initial fluxes can be set within kinetics process

    transport_roles = {
        'internal': ['g6p_c', 'pep_c', 'pyr_c', 'PTSG'],
        'external': ['glc__D_e'],
        }
    
    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}


# the composite function
def compose_kinetic_FBA(config):

    ## declare the processes
    # transport
    transport_config = copy.deepcopy(config)
    transport_config.update(get_transport_config())

    # transport_config.update({'target_fluxes': TARGET_FLUXES})
    transport = ConvenienceKinetics(transport_config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids

    # metabolism
    metabolism_config = copy.deepcopy(config)
    metabolism_config.update({
        'model_path': METABOLISM_FILE,
        'constrained_reactions': target_fluxes})
    metabolism = BiGGMetabolism(metabolism_config)

    # other processes
    deriver = DeriveVolume(config)
    division = Division(config)

    # place processes in layers.
    processes = [
        {'transport': transport},
        {'metabolism': metabolism},
        {'deriver': deriver,
        'division': division}
    ]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            # 'exchange': 'null',  # metabolism's exchange is used. TODO -- add exchange to convenience_kinetics
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

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
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


if __name__ == '__main__':
    from lens.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'kinetic_FBA_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # TODO -- load print emitter
    compartment = load_compartment(compose_kinetic_FBA)

    # settings for simulation and plot
    options = compose_kinetic_FBA({})['options']
    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 20}

    plot_settings = {
        'overlay': {'reactions': 'flux_bounds'},
        'max_rows': 25}

    # saved_state = simulate_compartment(compartment, settings)
    saved_state = simulate_with_environment(compartment, settings)
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # TODO -- make a flux network with metabolism

