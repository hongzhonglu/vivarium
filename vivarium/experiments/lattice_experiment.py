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
    growth_division_compartment = load_compartment(compose_growth_division, growth_division_config)
    lattice_compartment = load_compartment(compose_lattice_environment, lattice_config)


    import ipdb;
    ipdb.set_trace()


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = {}
    experiment = compose_lattice_experiment(config)
