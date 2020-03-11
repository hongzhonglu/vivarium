from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import (
    initialize_state,
    load_compartment,
)
from vivarium.compartment.composition import (
    simulate_compartment,
    convert_to_timeseries
)

# processes
from vivarium.processes.multibody_physics import (
    Multibody,
    n_body_config,
    plot_snapshots,
)
from vivarium.processes.diffusion_field import (
    DiffusionField,
    exchange_body_config,
    plot_field_output,
)


def compose_lattice_environment(config):
    """"""
    n_bodies = config.get('n_bodies', 10)
    bounds = config.get('bounds', [10, 10])
    size = config.get('size', [10, 10])
    molecules = config.get('molecules', ['glc'])

    ## Declare the processes.
    # multibody physics
    multibody_config = n_body_config(n_bodies, bounds)
    multibody = Multibody(multibody_config)

    # diffusion field
    bodies = {
        body_id: {
        'location': [bounds[0]/4, bounds[1]/4],
        'exchange': {
            mol_id: 1e2 for mol_id in molecules}}
            for body_id in range(n_bodies)}

    exchange_config = {
        'molecules': molecules,
        'n_bins': bounds,
        'size': size,
        'bodies': bodies}
    diffusion_config = exchange_body_config(exchange_config)
    diffusion = DiffusionField(diffusion_config)

    # Place processes in layers
    processes = [
        {'multibody': multibody,
        'diffusion': diffusion}]

    # topology
    topology = {
        'multibody': {
            'bodies': 'bodies',
        },
        'diffusion': {
            'bodies': 'bodies',
            'fields': 'fields'}}

    # initialize the states
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
    return {
        'molecules': ['glc'],
        'n_bodies': 10,
        'bounds': [10, 10],
        'size': [10, 10]}

def test_lattice_environment(config=get_lattice_config(), time=10):
    boot_config = {'emitter': 'null'}
    lattice_environment = load_compartment(compose_lattice_environment, boot_config)
    settings = {'total_time': time}
    return simulate_compartment(lattice_environment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_environment_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    saved_data = test_lattice_environment(config, 20)
    timeseries = convert_to_timeseries(saved_data)
    plot_field_output(timeseries, config, out_dir, 'lattice_field')
    plot_snapshots(timeseries, config, out_dir, 'lattice_bodies')
    