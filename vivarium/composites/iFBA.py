from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.derive_volume import DeriveVolume
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.BiGG_metabolism import BiGGMetabolism
from vivarium.processes.Kremling2007_transport import Transport
from vivarium.processes.CovertPalsson2002_regulation import Regulation

# target flux reaction names come from BiGG models
TARGET_FLUXES = ['GLCpts', 'PPS', 'PYK']  #, 'glc__D_e'] # TODO -- add exchange constraints


def compose_iFBA(config):
    '''
    TODO -- transport and metabolism have different names of environmental molecules.
    '''

    ## declare the processes
    # transport
    transport_config = copy.deepcopy(config)
    transport_config.update({'target_fluxes': TARGET_FLUXES})
    transport = Transport(transport_config)
    target_fluxes = transport.target_fluxes  # TODO -- just use TARGET_FLUXES?

    # metabolism
    metabolism_config = copy.deepcopy(config)
    metabolism_config.update({'constrained_reaction_ids': target_fluxes})
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
    from vivarium.actor.process import load_compartment, simulate_compartment, plot_simulation_output, simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'iFBA_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_iFBA)

    # get options
    options = compose_iFBA({})['options']
    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-9,  # L
        'timestep': 1,
        'total_time': 10}

    # saved_state = simulate_compartment(compartment, settings)
    saved_state = simulate_with_environment(compartment, settings)
    plot_simulation_output(saved_state, out_dir)
