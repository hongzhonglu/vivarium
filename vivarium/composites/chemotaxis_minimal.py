from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.compartment.composition import get_derivers, load_compartment
from vivarium.compartment.process import (
    initialize_state)
from vivarium.compartment.composition import (
    simulate_with_environment,
    plot_simulation_output)

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity

LIGAND_ID = 'MeAsp'
ENVIRONMENT_PORT = 'environment'



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
            'external': ENVIRONMENT_PORT,
            'internal': 'cell'},
        'motor': {
            'external': ENVIRONMENT_PORT,
            'internal': 'cell'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])

    # initialize the states
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    options = {
        'name': 'simple_chemotaxis_composite',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment_port': 'environment'}

    return {
        'processes': processes,
        'derivers': deriver_processes,
        'states': states,
        'options': options}


# testing function
def get_exponential_random_timeline(config):
    # exponential space with random direction changes
    time_total = config.get('time', 100)
    time_step = config.get('time_step', 1)
    base = config.get('base', 1+1e-4)  # mM/um
    speed = config.get('speed', 14)     # um/s
    forward_prob = config.get('forward_probability', 0.5)
    reverse_prob = 1 - forward_prob
    assert 0 <= forward_prob <= 1
    conc_0 = config.get('initial_conc', 0)  # mM

    conc = conc_0
    t = 0
    timeline = [(t, {ENVIRONMENT_PORT: {LIGAND_ID: conc}})]
    while t<=time_total:
        t += time_step
        direction = random.choices(
            population=[1, -1],
            weights=[forward_prob, reverse_prob])
        conc += base**(direction[0] * speed) - 1
        if conc<0:
            conc = 0
        timeline.append((t, {ENVIRONMENT_PORT: {LIGAND_ID: conc}}))

    return timeline

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'simple_chemotaxis_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # plot settings for the simulations
    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux'},
        'skip_ports': ['prior_state', 'null', 'global']}

    # exponential random timeline simulation
    time_step = 0.1
    exponential_random_config = {
        'time': 60,
        'time_step': time_step,
        'initial_conc': 0.01,
        'base': 1+6e-4,
        'speed': 14,
        'forward_probability': 0.52}
    timeline = get_exponential_random_timeline(exponential_random_config)
    boot_config = {
        'initial_ligand': timeline[0][1][ENVIRONMENT_PORT][LIGAND_ID],  # set initial_ligand from timeline
        'time_step': time_step}
    compartment = load_compartment(compose_simple_chemotaxis, boot_config)

    settings = {
        'environment_port': ENVIRONMENT_PORT,
        'environment_volume': 1e-13,  # L
        'timeline': timeline,
    }

    timeseries = simulate_with_environment(compartment, settings)
    plot_simulation_output(
        timeseries,
        plot_settings,
        out_dir,
        'exponential_timeline')
