from __future__ import absolute_import, division, print_function

from vivarium.actor.process import initialize_state, get_compartment_timestep

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity


def compose_simple_chemotaxis(config):

    # declare the processes
    receptor = ReceptorCluster(config)
    motor = MotorActivity(config)

    # place processes in layers
    processes_layers = [
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
    states = initialize_state(processes_layers, topology, config.get('initial_state', {}))

    # get the time step
    time_step = get_compartment_timestep(processes_layers)

    options = {
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'time_step': time_step,
        'environment_role': 'environment',
        # 'exchange_role': 'exchange',
    }

    return {
        'processes': processes_layers,
        'states': states,
        'options': options}
