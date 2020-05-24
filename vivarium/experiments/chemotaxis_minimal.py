from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.core.tree import (
    generate_state,
    Experiment)
from vivarium.core.composition import (
    make_agents,
    simulate_experiment
)

# compartments
from vivarium.compartments.lattice import (
    Lattice,
    get_lattice_config
)
from vivarium.compartments.chemotaxis_minimal import (
    ChemotaxisMinimal,
    get_chemotaxis_config
)

# processes
from vivarium.processes.multibody_physics import plot_snapshots



def make_chemotaxis_experiment(config={}):

    # get the environment
    env_config = config.get('environment', {})
    environment = Lattice(env_config)
    network = environment.generate({})
    processes = network['processes']
    topology = network['topology']

    # get the agents
    n_agents = config.get('n_agents', 1)
    chemotaxis = ChemotaxisMinimal({
        'external_key': ('..', 'external')})
    agents = make_agents(range(n_agents), chemotaxis, {})
    processes['agents'] = agents['processes']
    topology['agents'] = agents['topology']

    emitter = {'type': 'timeseries'}
    return Experiment({
        'processes': processes,
        'topology': topology,
        'emitter': emitter,
        'initial_state': config.get('initial_state', {})})


def run_chemotaxis_experiment(out_dir):
    time = 60
    config = {
        'environment': get_lattice_config(),
        'chemotaxis': get_chemotaxis_config({})}

    experiment = make_chemotaxis_experiment(config)

    settings = {
        'timestep': 1,
        'total_time': time,
        'return_raw_data': True}
    data = simulate_experiment(experiment, settings)

    # make snapshot plot
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_snapshots(agents, fields, config, out_dir, 'snapshots')

if __name__ == '__main__':
    out_dir = os.path.join('out', 'experiments', 'minimal_chemotaxis')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_chemotaxis_experiment(out_dir)
