from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.process import load_compartment

# compartments
from vivarium.composites.lattice_environment import compose_lattice_environment
from vivarium.composites.growth_division import compose_growth_division


def compose_lattice_experiment(config):

    # TODO -- lattice_config needs to get the correct agent_id linking it to growth_division
    # TODO -- lattice_config needs a list of subcompartment
    # TODO -- we need a initialize_embedded_compartment() composition function
    # growth_division_config = config.get('agents', {})
    lattice_config = config.get('lattice', {})
    agents_config = config.get('agents', {})


    # initialize agent compartments and get their boundaries
    n_agents = agents_config.get('n_agents', 0)
    agents = {}
    for agent in range(n_agents):
        compartment = load_compartment(compose_growth_division)
        boundary_store = compartment.states['environment']
        agents[str(uuid.uuid1())] = {
            'boundary': boundary_store}

    # load agent boundaries in lattice
    lattice_config.update({'agents': agents})
    lattice_compartment = load_compartment(compose_lattice_environment, lattice_config)


    import ipdb; ipdb.set_trace()


def get_ecoli_core_glc_config():
    # todo -- implement timeline
    lattice_config = {

        'size': (15, 15),
        'n_bins': (10, 10),
        'diffusion': 1e-3,
        'depth': 2e-2,
        'jitter_force': 1e-1}

    agent_config = {
        'n_agents': 1
    }

    return {
        'lattice': lattice_config,
        'agents': agent_config}

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_ecoli_core_glc_config()
    experiment = compose_lattice_experiment(config)
