from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.compartment.composition import get_derivers
from vivarium.compartment.process import (
    initialize_state,
    load_compartment)
from vivarium.compartment.composition import (
    simulate_with_environment,
    convert_to_timeseries,
    plot_simulation_output)

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity

LIGAND_ID = 'MeAsp'


def compose_simple_chemotaxis(config):

    initial_ligand = config.get('concentrations', {}).get(LIGAND_ID, 0.1)
    receptor_parameters = {
        'ligand': LIGAND_ID,
        'initial_ligand': initial_ligand}
    receptor_parameters.update(config)

    # declare the processes
    receptor = ReceptorCluster(receptor_parameters)
    motor = MotorActivity(config)

    # place processes in layers
    processes = [
        {'receptor': receptor},
        {'motor': motor}]

    # make the topology.
    # for each process, map process ports to store ids
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cell'},
        'motor': {
            'external': 'environment',
            'internal': 'cell'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])  # add deriver processes
    topology.update(derivers['deriver_topology'])  # add deriver topology

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'simple_chemotaxis_composite',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment_port': 'environment',
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}


# testing function
def get_exponential_random_timeline(config):
    # exponential space with random direction changes
    time = config.get('time', 100)
    base = config.get('base', 1+1e-4)  # mM/um
    speed = config.get('speed', 14)     # um/s
    conc_0 = config.get('initial_conc', 0)  # mM

    conc = conc_0
    timeline = [(0, {'environment': {LIGAND_ID: conc}})]
    for t in range(time):
        conc += base**(random.choice((-1, 1)) * speed) - 1
        if conc<0:
            conc = 0
        timeline.append((t, {'environment': {LIGAND_ID: conc}}))

    return timeline

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'simple_chemotaxis_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_simple_chemotaxis, boot_config)

    # settings for simulation and plot
    options = compartment.configuration
    # timeline = [(10, {})]

    exponential_random_config = {
        'time': 60,
        'base': 1+4e-4,
        'speed': 14}
    timeline = get_exponential_random_timeline(exponential_random_config)

    settings = {
        'environment_port': options['environment_port'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux'},
        'skip_ports': ['prior_state', 'null', 'global']}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))
    plot_simulation_output(timeseries, plot_settings, out_dir)
