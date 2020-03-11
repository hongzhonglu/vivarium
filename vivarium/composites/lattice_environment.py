from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import initialize_state
from vivarium.compartment.composition import (
    simulate_compartment,
    convert_to_timeseries
)

# processes
from vivarium.processes.multibody_physics import Multibody, n_body_config
from vivarium.processes.diffusion_field import DiffusionField, exchange_body_config


def compose_lattice_environment(config):
    """"""

    # Declare the processes.

    # multibody physics
    n_bodies = 10
    bounds = [10, 10]
    multibody_config = n_body_config(n_bodies, bounds)
    multibody = Multibody(multibody_config)

    # diffusion field
    molecules = ['glc']
    # TODO -- get bodies from multibody
    bodies = {
        body_id: {
        'location': [bounds[0]/4, bounds[1]/4],
        'exchange': {
            mol_id: 1e2 for mol_id in molecules}}
            for body_id in range(n_bodies)}

    import ipdb; ipdb.set_trace()

    config = {
        'molecules': molecules,
        'n_bins': bounds,
        'size': bounds,
        'bodies': bodies}
    diffusion_config = exchange_body_config(config)
    diffusion = DiffusionField(diffusion_config)

    # Place processes in layers
    processes = [
        {'multibody': multibody},
        {'diffusion': diffusion}]

    import ipdb; ipdb.set_trace()
    # TODO -- link up locations to same ports

    # topology
    topology = {
        'diffusion': {},
        'multibody': {}}

    # Initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': config.get('name', 'lattice_environment'),
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0)}

    return {
        'processes': processes,
        'states': states,
        'options': options}



# toy functions/ defaults
def get_lattice_config():

    return {}

def test_lattice_environment(config=get_lattice_config(), time=10):
    lattice_environment = compose_lattice_environment(config)
    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-2}
    return simulate_compartment(lattice_environment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_environment_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    saved_data = test_lattice_environment(config, 20)
    timeseries = convert_to_timeseries(saved_data)

    # # settings for simulation and plot
    # options = compartment.configuration
    # timeline = [(2520, {})]
    #
    # settings = {
    #     'environment_port': options['environment_port'],
    #     'exchange_port': options['exchange_port'],
    #     'environment_volume': 1e-13,  # L
    #     'timeline': timeline}
    #
    # plot_settings = {
    #     'max_rows': 20,
    #     'remove_zeros': True,
    #     'overlay': {
    #         'reactions': 'flux'},
    #     'skip_ports': ['prior_state', 'null']}
    #
    # # saved_state = simulate_compartment(compartment, settings)
    # saved_data = simulate_with_environment(compartment, settings)
    # del saved_data[0]
    # timeseries = convert_to_timeseries(saved_data)
    # volume_ts = timeseries['global']['volume']
    # print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    # plot_simulation_output(timeseries, plot_settings, out_dir)
