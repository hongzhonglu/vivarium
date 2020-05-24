from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.core.tree import Compartment
from vivarium.core.composition import (
    compartment_in_experiment,
    simulate_with_environment,
    plot_simulation_output)

# processes
from vivarium.processes.Endres2006_chemoreceptor import (
    ReceptorCluster,
    get_exponential_random_timeline
)
from vivarium.processes.Vladimirov2008_motor import MotorActivity



class ChemotaxisMinimal(Compartment):

    defaults = {
        'ligand_id': 'MeAsp',
        'initial_ligand': 0.1,
        'external_key': ('..', 'external',)
    }

    def __init__(self, config):
        self.config = config
        self.ligand_id = config.get(
            'ligand_id',
            self.defaults['ligand_id'])
        self.initial_ligand = config.get(
            'initial_ligand',
            self.defaults['initial_ligand'])
        self.external_key = self.config.get(
            'external_key',
            self.defaults['external_key'])

    def generate_processes(self, config):
        receptor_parameters = {
            'ligand': self.ligand_id,
            'initial_ligand': self.initial_ligand}

        # declare the processes
        receptor = ReceptorCluster(receptor_parameters)
        motor = MotorActivity({})

        return {
            'receptor': receptor,
            'motor': motor}

    def generate_topology(self, config):
        return {
            'receptor': {
                'external': self.external_key,
                'internal': ('cell',)},
            'motor': {
                'external': self.external_key,
                'internal': ('cell',)}}




def get_chemotaxis_config(config={}):
    ligand_id = config.get('ligand_id', 'MeAsp')
    initial_ligand = config.get('initial_ligand', 5.0)
    external_key = config.get('external_key', 'external')
    # configure the compartment
    return {
        'external_key': (external_key,),
        'ligand_id': ligand_id,
        'initial_ligand': initial_ligand}  # set initial_ligand from timeline


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'chemotaxis_minimal')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ligand_id = 'MeAsp'
    environment_port = 'external'

    exponential_random_config = {
        'ligand': ligand_id,
        'environment_port': environment_port,
        'time': 60,
        'timestep': 1,
        'base': 1+4e-4,
        'speed': 14}
    timeline = get_exponential_random_timeline(exponential_random_config)
    end_time = timeline[-1][0]


    config = {
        'ligand_id': ligand_id,
        'initial_ligand': timeline[0][1][(environment_port, ligand_id)],
        'external_key': environment_port}
    compartment = ChemotaxisMinimal(get_chemotaxis_config(config))

    # configure experiment
    experiment_settings = {
        'timeline': timeline,
        'timeline_port_mapping': {
            environment_port: (environment_port,)}}
    experiment = compartment_in_experiment(compartment, experiment_settings)

    # run experiment
    timestep = 1
    time = 0
    while time < end_time:
        experiment.update(timestep)
        time += timestep
    timeseries = experiment.emitter.get_timeseries()

    # plot settings for the simulations
    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux'},
        'skip_ports': ['prior_state', 'null', 'global']}
    plot_simulation_output(
        timeseries,
        plot_settings,
        out_dir,
        'exponential_timeline')
