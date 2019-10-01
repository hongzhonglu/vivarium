from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_initial_states
from lens.actor.emitter import get_emitter, configure_emitter
from lens.environment.lattice_compartment import LatticeCompartment

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.growth import Growth


def initialize_growth_division(config):
    config.update({
        'emitter': {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
        }
    })

    # declare the processes
    growth = Growth(config)
    deriver = DeriveVolume(config)
    # division = Division(config)
    processes = {
        'growth': growth,
        'deriver': deriver,
        # 'division': division
    }

    # configure the states to the roles for each process
    topology = {
        'growth': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        # 'division': {
        #     'internal': 'cell'},
        }

    # initialize the states
    # TODO (Eran) -- this can be refactored
    initial_state = merge_initial_states(processes)
    states = {
        # 'environment': State(initial_state['external']),
        'cell': State(initial_state['internal'])}

    # configure emitter
    emitter = configure_emitter(config, processes, topology)

    options = {
        'topology': topology,
        'emitter': emitter,
        'environment': 'environment',
        'compartment': 'cell',
        'exchange_key': config['exchange_key'],
        'environment_ids': initial_state['environment_ids'],
        'environment_deltas': initial_state['environment_deltas']}

    # create the compartment
    return LatticeCompartment(processes, states, options)