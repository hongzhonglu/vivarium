from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states
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
    default_state = merge_default_states(processes)
    initial_state = config.get('initial_state', {})
    initial_time = config.get('initial_time', 0.0)
    default_state['internal'].update(initial_state['cell'])

    states = {
        # 'environment': State(initial_state['external']),
        'cell': State(initial_state=default_state['internal'])}

    # configure emitter
    emitter = configure_emitter(config, processes, topology)

    options = {
        'topology': topology,
        'emitter': emitter,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'exchange_key': config['exchange_key'],
        'environment_ids': initial_state['environment_ids'],
        'environment_deltas': initial_state['environment_deltas']}

    # create the compartment
    return LatticeCompartment(processes, states, options)
