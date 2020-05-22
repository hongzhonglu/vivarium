from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.tree import (
    generate_state,
    Experiment)
from vivarium.compartment.composition import (
    make_agents,
    simulate_experiment
)

# compartments
from vivarium.composites.lattice_environment import Lattice
from vivarium.composites.chemotaxis_minimal import ChemotaxisMinimal



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

def run_chemotaxis_experiment(config={}, time=10):
    experiment = make_chemotaxis_experiment(config)

    settings = {
        'timestep': 1,
        'total_time': time,
        'return_raw_data': True}
    return simulate_experiment(experiment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'experiments', 'minimal_chemotaxis')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    timeseries = run_chemotaxis_experiment()
