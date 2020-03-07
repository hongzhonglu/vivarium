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
    time_total = config.get('time', 100)
    timestep = config.get('timestep', 1)
    base = config.get('base', 1+1e-4)  # mM/um
    speed = config.get('speed', 14)     # um/s
    conc_0 = config.get('initial_conc', 0)  # mM

    conc = conc_0
    t = 0
    timeline = [(t, {ENVIRONMENT_PORT: {LIGAND_ID: conc}})]
    while t<=time_total:
        t += timestep
        conc += base**(random.choice((-1, 1)) * speed) - 1
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
    exponential_random_config = {
        'time': 60,
        'timestep': 0.1,
        'base': 1+4e-4,
        'speed': 14}
    timeline = get_exponential_random_timeline(exponential_random_config)
    boot_config = {
        'initial_ligand': timeline[0][1][ENVIRONMENT_PORT][LIGAND_ID],  # set initial_ligand from timeline
        'emitter': 'null'}
    compartment = load_compartment(compose_simple_chemotaxis, boot_config)

    settings = {
        'environment_port': ENVIRONMENT_PORT,
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir, 'exponential_timeline')


    # # null timeline simulation
    # null_timeline = [
    #     (0, {ENVIRONMENT_PORT: {LIGAND_ID: 0.0}}),
    #     (20, {})]
    # null_boot_config = {
    #     'initial_ligand': null_timeline[0][1][ENVIRONMENT_PORT][LIGAND_ID],  # set initial_ligand from timeline
    #     'emitter': 'null'}
    # compartment2 = load_compartment(compose_simple_chemotaxis, null_boot_config)
    #
    # null_settings = {
    #     'environment_port': ENVIRONMENT_PORT,
    #     'environment_volume': 1e-13,  # L
    #     'timeline': null_timeline}
    #
    # null_saved_data = simulate_with_environment(compartment2, null_settings)
    # del null_saved_data[0]
    # null_timeseries = convert_to_timeseries(null_saved_data)
    # plot_simulation_output(null_timeseries, plot_settings, out_dir, 'null_timeline')
