from __future__ import absolute_import, division, print_function

from lens.actor.process import merge_default_updaters, initialize_state
from lens.utils.dict_utils import merge_dicts

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity


def compose_vladimirov_chemotaxis(config):
    exchange_key = config.get('exchange_key')

    # declare the processes
    receptor = ReceptorCluster(config)
    motor = MotorActivity(config)

    # place processes in layers
    processes = [
        {'receptor': receptor},
        {'motor': motor}]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cell'},
        'motor': {
            'external': 'environment',
            'internal': 'cell'},
        }

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    # get environment_ids -- TODO remove this!
    default_updaters = merge_default_updaters(processes)
    environment_ids = []
    for role, state_ids in default_updaters.items():
        for state_id, updater in state_ids.items():
            if updater is 'accumulate':
                environment_ids.append(state_id)

    options = {
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
