from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.tree import generate_state, Experiment
from vivarium.compartment.composition import (
    simulate_compartment,
    load_compartment,
    get_derivers
)

from vivarium.processes.multibody_physics import plot_snapshots
from vivarium.processes.diffusion_field import plot_field_output

# compartments
from vivarium.composites.lattice_environment import (
    make_lattice_environment,
)
from vivarium.composites.growth_division import growth_division

from vivarium.composites.growth import growth_compartment



def make_agents(count, compartment, config):
    processes = {}
    topology = {}

    for agent in range(count):
        # agent_id = str(uuid.uuid1())
        agent_id = str(agent)

        # make the agent
        agent_config = config.copy()
        agent_config.update({'agent_id': agent_id})
        agent = compartment(dict(
            agent_config,
            agent_id=agent_id))

        # save processes and topology
        processes[agent_id] = {
            'cell': agent['processes']}
        topology[agent_id] = {
            'cell': agent['topology']}

    return {
        'processes': processes,
        'topology': topology}


def lattice_experiment(config):
    # configure the experiment
    count = config.get('count')

    # get the environment
    environment = make_lattice_environment(config.get('environment', {}))
    processes = environment['processes']
    topology = environment['topology']

    agents = make_agents(count, growth_compartment, {
        'cells_key': ['..', 'agents']})
    processes['agents'] = agents['processes']
    topology['agents'] = agents['topology']

    experiment = Experiment({
        'processes': processes,
        'topology': topology,
        'initial_state': config.get('initial_state', {})})

    print('processes ------------------------')
    print(experiment.processes)

    print('topology ------------------------')
    print(experiment.topology)

    print('before ------------------------')
    print(experiment.state.get_value())

    experiment.update(10.0)

    print('after ----------------------------------')
    print(experiment.state.get_value())

    import ipdb; ipdb.set_trace()


# toy functions/ defaults
def get_lattice_config():

    environment_config = {
        'molecules': ['glc'],
        'bounds': [10, 10],
        'size': [10, 10],
    }

    agent_config = {}

    return {
        'count': 3,
        'environment': environment_config,
        'agents': agent_config
    }



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    lattice_experiment(config)