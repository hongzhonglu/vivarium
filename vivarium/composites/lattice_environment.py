from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import initialize_state
from vivarium.compartment.composition import (
    simulate_compartment,
    load_compartment
)

# processes
from vivarium.processes.multibody_physics import (
    Multibody,
    random_body_config,
    plot_snapshots,
    DEFAULT_BOUNDS
)
from vivarium.processes.diffusion_field import (
    DiffusionField,
    exchange_agent_config,
)



def get_lattice_environment(config):
    """
    Lattice environment

    A two-dimensional lattice environmental model

    """

    # declare the processes.
    multibody = Multibody(config.get('multibody', {}))
    diffusion = DiffusionField(config.get('diffusion_field', {}))

    # place processes in layers
    processes = [
        {'multibody': multibody,
        'diffusion': diffusion}]

    # topology
    topology = {
        'multibody': {
            'agents': 'boundary',
        },
        'diffusion': {
            'agents': 'boundary',
            'fields': 'fields'}}

    return {
        'processes': processes,
        'topology': topology}



def compose_lattice_environment(config):
    """"""
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    size = config.get('size', DEFAULT_BOUNDS)
    molecules = config.get('molecules', ['glc'])

    # configure the agents
    agents_config = config.get('agents', {})
    agent_ids = list(agents_config.keys())
    n_agents = len(agent_ids)
    multibody_config = {
        'n_agents': n_agents,
        'bounds': bounds}
    multibody_config.update(random_body_config(multibody_config))

    # config for the diffusion proces
    agents = {
        agent_id: {
        'location': boundary['location'],
        'exchange': {
            mol_id: 1e2 for mol_id in molecules}}  # TODO -- don't hardcode exchange
            for agent_id, boundary in multibody_config['agents'].items()}
    exchange_config = {
        'molecules': molecules,
        'n_bins': bounds,
        'size': size,
        'agents': agents
    }
    diffusion_config = exchange_agent_config(exchange_config)

    # environment gets both process configs
    environment_config = {
        'multibody': multibody_config,
        'diffusion_field': diffusion_config
    }

    # get the environment compartment
    environment_compartment = get_lattice_environment(environment_config)
    processes = environment_compartment['processes']
    topology = environment_compartment['topology']

    # add derivers
    deriver_processes = []

    # initialize the states
    states = initialize_state(
        processes,
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
def get_lattice_config(config={'n_agents':2}):
    body_config = random_body_config(config)
    body_config.update({
        'molecules': ['glc'],
        'bounds': DEFAULT_BOUNDS,
        'size': DEFAULT_BOUNDS})  # no boundary store

    return body_config

def test_lattice_environment(config=get_lattice_config(), time=10):
    initial_agents_state = config['agents']
    compartment = load_compartment(compose_lattice_environment, config)

    # initialize agent boundary state
    compartment.send_updates({'boundary': [{'agents': initial_agents_state}]})

    settings = {
        'return_raw_data': True,
        'total_time': time}
    return simulate_compartment(compartment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_environment_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    data = test_lattice_environment(config, 10)

    # make snapshot plot
    agents = {time: time_data['boundary']['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_snapshots(agents, fields, config, out_dir, 'snapshots')
