from __future__ import absolute_import, division, print_function

import copy

from lens.actor.process import initialize_state

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.division import Division, divide_condition, divide_state
from lens.processes.BiGG_metabolism import BiGGMetabolism
from lens.processes.Kremling2007_transport import Transport
from lens.processes.CovertPalsson2002_regulation import Regulation

# target flux reaction names come from BiGG models
TARGET_FLUXES = ['GLCpts', 'PPS', 'PYK']  #, 'glc__D_e'] # TODO -- add exchange constraints


def compose_iFBA(config):

    ## declare the processes
    # transport
    transport_config = copy.deepcopy(config)
    transport_config.update({'target_fluxes': TARGET_FLUXES})
    transport = Transport(transport_config)
    target_fluxes = transport.target_fluxes  # TODO -- just use TARGET_FLUXES?

    # metabolism
    metabolism_config = copy.deepcopy(config)
    metabolism_config.update({'constrained_reactions': target_fluxes})
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


def test_iFBA():
    import numpy as np
    from lens.actor.process import Compartment
    from lens.environment.lattice_compartment import LatticeCompartment

    boot_config = {}
    composite_config = compose_iFBA(boot_config)
    processes = composite_config['processes']
    states = composite_config['states']
    options = composite_config['options']

    # make compartment
    compartment = LatticeCompartment(processes, states, options)

    print(compartment.current_parameters())
    print(compartment.current_state())

    # test compartment
    compartment = Compartment(processes, states, options)

    print('compartment current_parameters: {}'.format(compartment.current_parameters()))
    print('compartment current_state: {}'.format(compartment.current_state()))

    # make lattice_compartment
    lattice_compartment = LatticeCompartment(processes, states, options)

    print(lattice_compartment.current_parameters())
    print(lattice_compartment.current_state())

    # evaluate compartment
    timestep = 1
    for steps in np.arange(10):
        lattice_compartment.update(timestep)
        print('lattice_compartment current_state: {}'.format(lattice_compartment.current_state()))


if __name__ == '__main__':
    saved_state = test_iFBA()
