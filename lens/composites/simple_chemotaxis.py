from __future__ import absolute_import, division, print_function

from lens.actor.process import initialize_state

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity


def compose_simple_chemotaxis(config):

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

    options = {
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment_role': 'environment',
        # 'exchange_role': 'exchange',
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
