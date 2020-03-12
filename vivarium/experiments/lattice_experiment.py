from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import load_compartment

# compartments
from vivarium.composites.lattice_environment import compose_lattice_environment
from vivarium.composites.growth_division import compose_growth_division


def compose_lattice_experiment(config):

    # TODO -- lattice_config needs to get the correct agent_id linking it to growth_division
    # TODO -- lattice_config needs a list of subcompartment
    # TODO -- we need a initialize_embedded_compartment() composition function
    growth_division_config = {}
    lattice_config = {}

    # make the compartments
    n_agents = config.get('n_agents', 0)

    agents = {}
    for agent in range(n_agents):
        # TODO -- uuid?
        agents[agent] = load_compartment(compose_growth_division, growth_division_config)

    # TODO -- load agent ids in lattice
    lattice_compartment = load_compartment(compose_lattice_environment, lattice_config)


    import ipdb;
    ipdb.set_trace()


def get_ecoli_core_glc_config():
    timeline_str = '0 ecoli_core_GLC 1.0 L + lac__D_e 1.0 mmol 0.1 L, 21600 end'
    lattice_config = {
        'name': 'ecoli_core',
        'timeline_str': timeline_str,  # todo -- implement timeline
        'size': (15, 15),
        'n_bins': (10, 10),
        'diffusion': 1e-3,
        'depth': 2e-2,
        'jitter_force': 1e-1,
        'n_agents': 1,
        # 'run_for': 5.0,
        # 'emit_fields': [
        #     'co2_e',
        #     'o2_e',
        #     'glc__D_e',
        #     'lac__D_e']
    }

    return lattice_config

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_ecoli_core_glc_config()
    experiment = compose_lattice_experiment(config)
