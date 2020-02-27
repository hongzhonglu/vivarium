from __future__ import absolute_import, division, print_function

from vivarium.actor.process import initialize_state, get_compartment_timestep
from vivarium.compartment.composition import get_schema

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity


def compose_simple_chemotaxis(config):

    receptor_parameters = {'ligand': 'MeAsp'}
    receptor_parameters.update(config)

    # declare the processes
    receptor = ReceptorCluster(receptor_parameters)
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
            'internal': 'cell'}}

    # get schema
    schema = get_schema(processes, topology)

    # initialize the states
    states = initialize_state(processes, topology, schema, config.get('initial_state', {}))

    # get the time step
    time_step = get_compartment_timestep(processes)

    options = {
        'name': 'simple_chemotaxis_composite',
        'topology': topology,
        'schema': schema,
        'initial_time': config.get('initial_time', 0.0),
        'time_step': time_step,
        'environment_role': 'environment',
        # 'exchange_role': 'exchange',
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
