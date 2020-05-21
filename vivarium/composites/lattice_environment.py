from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.tree import Compartment
from vivarium.compartment.composition import compartment_in_experiment

# processes
from vivarium.processes.multibody_physics import (
    Multibody,
    random_body_config,
    plot_snapshots,
)
from vivarium.processes.diffusion_field import (
    DiffusionField,
    exchange_agent_config,
)



class Lattice(Compartment):
    """
    Lattice:  A two-dimensional lattice environmental model with multibody physics and diffusing molecular fields.
    """

    defaults = {
        'multibody': {
            'bounds': [10, 10],
            'size': [10, 10],
            'agents': {}
        },
        'diffusion': {
            'molecules': ['glc'],
            'n_bins': [10, 10],
            'size': [10, 10],
            'depth': 3000.0,  # um
            'diffusion': 5e-1,
        }
    }

    def __init__(self, config):
        self.multibody_config = config.get('multibody', self.defaults['multibody'])
        self.diffusion_config = config.get('diffusion', self.defaults['diffusion'])

    def generate_processes(self, config):
        multibody = Multibody(config.get(
            'multibody',
            self.multibody_config))
        diffusion = DiffusionField(config.get(
            'diffusion_field',
            self.diffusion_config))

        return {
            'multibody': multibody,
            'diffusion': diffusion}

    def generate_topology(self, config):
        return {
            'multibody': {
                'agents': ('agents',)},
            'diffusion': {
                'agents': ('agents',),
                'fields': ('fields',)}}


def get_lattice_config(bounds=[25,25]):

    # multibody confid
    mbp_config = {
        # 'animate': True,
        'jitter_force': 1e0,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'n_agents': 2}
    mbp_config.update(random_body_config(body_config))

    # diffusion config
    dff_mol = 'glc'
    dff_config = {
        'molecules': [dff_mol],
        'n_bins': bounds,
        'size': bounds,
        'diffusion': 1e-10,
        'gradient': {
            'type': 'gaussian',
            'molecules': {
                dff_mol:{
                    'center': [0.5, 0.5],
                    'deviation': 5}}}}

    return {
        'bounds': bounds,
        'multibody': mbp_config,
        'diffusion': dff_config}

def test_lattice(config=get_lattice_config(), end_time=10):

    # configure the compartment
    compartment = Lattice(config)

    # configure experiment
    experiment_settings = {
        'compartment': config}
    experiment = compartment_in_experiment(
        compartment,
        experiment_settings)

    # run experiment
    timestep = 1
    time = 0
    while time < end_time:
        experiment.update(timestep)
        time += timestep
    return experiment.emitter.get_data()



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_compartment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    data = test_lattice(config, 40)

    # make snapshot plot
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    plot_snapshots(agents, fields, config, out_dir, 'snapshots')
