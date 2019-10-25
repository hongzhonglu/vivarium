from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.growth import Growth


def compose_growth_division(config):

    # declare the processes
    growth = Growth(config)
    deriver = DeriveVolume(config)
    processes = {
        'growth': growth,
        'deriver': deriver,
    }

    # configure the states to the roles for each process
    topology = {
        'growth': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # initialize the states
    default_state = merge_default_states(processes)
    initial_state = config.get('initial_state', {})
    initial_time = config.get('initial_time', 0.0)
    default_state['internal'].update(initial_state['cell'])

    states = {
        # 'environment': State(initial_state['external']),
        'cell': State(initial_state=default_state['internal'])}

    options = {
        'topology': topology,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': initial_state['environment_ids']}

    return {
        'processes': processes,
        'states': states,
        'options': options}
