from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.process import (
    initialize_state,
    flatten_process_layers
)
from vivarium.compartment.composition import (
    simulate_compartment,
    load_compartment,
    get_derivers
)

from vivarium.processes.multibody_physics import plot_snapshots
from vivarium.processes.diffusion_field import plot_field_output

# compartments
from vivarium.composites.lattice_environment import get_environment
from vivarium.composites.growth_division import growth_division


BOUNDARY_ID = 'boundary'

# TODO -- this can be made into a general function
def get_agents(n_agents):
    processes = []
    topologies = {}
    agent_ids = []
    for agent in range(n_agents):

        agent_id = str(uuid.uuid1())
        agent_ids.extend(agent_id)

        # make the agent
        agent = growth_division(config.get('agents', {}))  # TODO -- pass in compartment

        # processes -- each process id is associated with its agent id with a tuple
        a_processes = flatten_process_layers(agent['processes'])
        a_processes = {
            (agent_id, process_id): process
            for process_id, process in a_processes.items()}

        # topology
        a_topology = {}
        for process_id, topology in agent['topology'].items():
            ports = {}
            for port_id, store_id in topology.items():
                if store_id == BOUNDARY_ID:
                    ports[port_id] = store_id
                else:
                    ports[port_id] = (agent_id, store_id)

            a_topology[(agent_id, process_id)] = ports

        # save processes and topology
        processes.append(a_processes)
        topologies[agent_id] = a_topology

    return {
        'ids': agent_ids,
        'processes': processes,
        'topologies': topologies}


# TODO -- this can move to a separate experiments directory
def lattice_experiment(config):
    # configure the experiment
    n_agents = config.get('n_agents')

    # get the environment
    environment = get_environment(config.get('environment', {}))
    environment_processes = flatten_process_layers(environment['processes'])
    environment_topology = environment['topology']
    inner_key = 'agents'  # TODO -- get this from config of each env process

    # get agent processes and topologies
    agents = get_agents(n_agents)
    agent_processes = agents['processes']
    agent_topologies = agents['topologies']
    agent_ids = agents['ids']

    ## make processes and topology for experiment
    processes = []
    topology = {}

    # add environment
    processes.append(environment_processes)
    topology.update(environment_topology)

    # add agents
    processes.extend(agent_processes)

    # combine agent and environment topologies
    for agent_id, agent_topology in agent_topologies.items():
        topology.update(agent_topology)

    ## add derivers
    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])  # add derivers to the topology

    # initialize the states
    # TODO -- pull out each agent_boundary, make a special initialize_state that can connect these up
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    options = {
        'name': config.get('name', 'lattice_environment'),
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0)}

    return {
        'processes': processes,
        'derivers': deriver_processes,
        'states': states,
        'options': options}



# toy functions/ defaults
def get_lattice_config():

    environment_config = {
        'molecules': ['glc'],
        'bounds': [10, 10],
        'size': [10, 10],
    }

    agent_config = {}

    return {
        'n_agents': 5,
        'environment': environment_config,
        'agents': agent_config
    }

def test_lattice_experiment(config=get_lattice_config(), time=10):
    lattice_environment = load_compartment(lattice_experiment, config)
    settings = {'total_time': time}
    return simulate_compartment(lattice_environment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    timeseries = test_lattice_experiment(config, 10)
    plot_field_output(timeseries, config, out_dir, 'lattice_field')
    plot_snapshots(timeseries, config, out_dir, 'lattice_bodies')
    