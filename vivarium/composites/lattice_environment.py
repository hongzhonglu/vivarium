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
    get_n_dummy_agents,
    random_body_config,
    plot_snapshots,
)
from vivarium.processes.diffusion_field import (
    DiffusionField,
    exchange_agent_config,
    plot_field_output,
)



def lattice_environment(config):
    bounds = config.get('bounds', [10, 10])
    size = config.get('size', [10, 10])
    molecules = config.get('molecules', ['glc'])

    # get the agents
    agents_config = config.get('agents', {})
    agent_ids = list(agents_config.keys())

    ## Declare the processes.
    # multibody physics
    multibody_config = random_body_config(agents_config, bounds)
    multibody = Multibody(multibody_config)

    # diffusion field
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
        'agents': agents}
    diffusion_config = exchange_agent_config(exchange_config)
    diffusion = DiffusionField(diffusion_config)

    # Place processes in layers
    processes = [
        {'multibody': multibody,
        'diffusion': diffusion}]

    # topology
    topology = {
        'multibody': {
            'agents': 'agents',
        },
        'diffusion': {
            'agents': 'agents',
            'fields': 'fields'}}

    return {
        'processes': processes,
        'topology': topology}




def lattice_environment_experiment(config):
    environment = lattice_environment(config)
    processes = environment['processes']
    topology = environment['topology']

    # add derivers
    deriver_processes = []

    # initialize the states
    # TODO -- pull out each agent_boundary, make a special initialize_state that can connect these up
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
def get_lattice_config():
    return {
        'molecules': ['glc'],
        'bounds': [10, 10],
        'size': [10, 10],
        'agents': get_n_dummy_agents(6)}  # no boundary store

def test_lattice_environment(config=get_lattice_config(), time=10):
    lattice_environment = load_compartment(lattice_environment_experiment, config)
    settings = {'total_time': time}
    return simulate_compartment(lattice_environment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_environment_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    timeseries = test_lattice_environment(config, 10)
    plot_field_output(timeseries, config, out_dir, 'lattice_field')
    plot_snapshots(timeseries, config, out_dir, 'lattice_bodies')
    