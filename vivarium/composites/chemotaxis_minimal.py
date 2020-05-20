from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.compartment.tree import Compartment
from vivarium.compartment.composition import (
    compartment_in_experiment,
    simulate_with_environment,
    plot_simulation_output)

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity



class ChemotaxisMinimal(Compartment):

    defaults = {
        'ligand_id': 'MeAsp',
        'initial_ligand': 0.1,
        'external_key': ('external',)
    }

    def __init__(self, config):
        self.config = config
        self.ligand_id = config.get(
            'ligand_id',
            self.defaults['ligand_id'])
        self.initial_ligand = config.get(
            'initial_ligand',
            self.defaults['initial_ligand'])
        self.external_key = ('..',) + self.config.get(
            'external_key',
            self.defaults['external_key'])

    def generate_processes(self, config):
        receptor_parameters = {
            'ligand': self.ligand_id,
            'initial_ligand': self.initial_ligand}
        receptor_parameters.update(config)

        # declare the processes
        receptor = ReceptorCluster(receptor_parameters)
        motor = MotorActivity(config)

        return {
            'receptor': receptor,
            'motor': motor}

    def generate_topology(self, config):
        external_key = config.get('external_key', self.external_key)
        return {
            'receptor': {
                'external': external_key,
                'internal': 'cell'},
            'motor': {
                'external': external_key,
                'internal': 'cell'}}


# testing functions
def get_exponential_random_timeline(config):
    # exponential space with random direction changes
    ligand_id = config.get('ligand_id', 'MeAsp')
    environment_port = config.get('environment_port', 'environment')
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
    timeline = [(t, {environment_port: {ligand_id: conc}})]
    while t<=time_total:
        t += time_step
        direction = random.choices(
            population=[1, -1],
            weights=[forward_prob, reverse_prob])
        conc += base**(direction[0] * speed) - 1
        if conc<0:
            conc = 0
        timeline.append((t, {environment_port: {ligand_id: conc}}))

    return timeline

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'chemotaxis_minimal')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make and exponential random timeline
    ligand_id = 'MeAsp'
    environment_port = 'environment'
    exponential_random_config = {
        'ligand_id': ligand_id,
        'environment_port': environment_port,
        'time': 60,
        'time_step': 0.1,
        'initial_conc': 0.01,
        'base': 1+6e-4,
        'speed': 14,
        'forward_probability': 0.52}
    timeline = get_exponential_random_timeline(exponential_random_config)

    # configure the compartment and experiment
    compartment_config = {
        'ligand_id': ligand_id,
        'initial_ligand': timeline[0][1][environment_port][ligand_id]}  # set initial_ligand from timeline
    compartment = ChemotaxisMinimal(compartment_config)

    experiment_settings = {}
    experiment = compartment_in_experiment(compartment, experiment_settings)

    settings = {
        'environment_port': environment_port,
        'environment_volume': 1e-13,  # L
        'timeline': timeline,
    }

    import ipdb; ipdb.set_trace()

    timeseries = simulate_with_environment(compartment, settings)

    # # plot settings for the simulations
    # plot_settings = {
    #     'max_rows': 20,
    #     'remove_zeros': True,
    #     'overlay': {
    #         'reactions': 'flux'},
    #     'skip_ports': ['prior_state', 'null', 'global']}
    # plot_simulation_output(
    #     timeseries,
    #     plot_settings,
    #     out_dir,
    #     'exponential_timeline')
