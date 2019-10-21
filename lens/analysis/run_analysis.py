from __future__ import absolute_import, division, print_function

import os
import argparse
from pymongo import MongoClient

from lens.analysis.analyze_compartment import Compartment
from lens.analysis.analyze_lattice import LatticeTrace
from lens.analysis.chemotaxis_trace import ChemotaxisTrace
from lens.analysis.snapshots import Snapshots


# classes for run_analysis to cycle through
analysis_classes = {
    'chemotaxis': ChemotaxisTrace,
    'compartment': Compartment,
    'lattice': LatticeTrace,
    'snapshots': Snapshots,
}

# mongoDB local url
url='localhost:27017'

def get_experiment(data):
    experiment_config = {}
    for row in data:
        if row.get('type') == 'lattice':
            experiment_config['edge_length'] = row['edge_length']
            experiment_config['patches_per_edge'] = row['patches_per_edge']
            experiment_config['cell_radius'] = row['cell_radius']
        elif row.get('topology'):
            experiment_config['topology'] = row['topology']

    return experiment_config

def get_sims_from_exp(client, experiment_id):
    # given a database client and experiment id, return a list of simulation ids
    simulation_ids = set()
    query_dict = {'experiment_id': experiment_id}  # TODO -- narrow query further. {'time': 1.0} not sufficient
    data = client.find(query_dict)
    for row in data:
        simulation_id = row.get('simulation_id')
        if simulation_id:
            simulation_ids.add(simulation_id)

    return list(simulation_ids)

class Analyze(object):
    '''
    example use:
        python lens/analysis/run_analysis.py -e lattice
    '''

    def __init__(self):

        parser = argparse.ArgumentParser(description='analyze rate laws')
        parser = self.add_arguments(parser)
        args = parser.parse_args()

        self.client = MongoClient(url)

        self.path = args.path
        self.experiment_id = args.experiment
        self.analyses = args.analyses

    def run_analysis(self):

        config_client = self.client.simulations.configuration
        history_client = self.client.simulations.history


        simulation_ids = get_sims_from_exp(self.client.simulations.history, self.experiment_id)
        output_dir = os.path.join(self.path, self.experiment_id)

        # make plot output directories
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        for sim_id in simulation_ids:
            sim_out_dir = os.path.join(output_dir, sim_id)
            if not os.path.isdir(sim_out_dir):
                os.makedirs(sim_out_dir)

        # get the experiment configuration
        query = {'experiment_id': self.experiment_id}
        config_data = config_client.find(query)
        experiment_config = get_experiment(config_data)
        active_processes = experiment_config['topology'].keys()


        # make list of analysis objects to try
        if self.analyses:
            run_analyses = [analysis_classes[analysis_id] for analysis_id in self.analyses]
        else:
            run_analyses = analysis_classes.values()

        # run analyses that declare the active processes
        for analysis_class in run_analyses:
            analysis = analysis_class()
            required_processes = analysis.requirements()

            if all(processes in active_processes for processes in required_processes):
                data = analysis.get_data(history_client)

                if analysis.analysis_type is 'compartment':
                    for sim_id in simulation_ids:
                        sim_out_dir = os.path.join(output_dir, sim_id)
                        analysis.analyze(experiment_config, data, sim_out_dir)

                elif analysis.analysis_type is 'lattice':
                    analysis.analyze(experiment_config, data, output_dir)

                print('completed analysis: {}'.format(analysis_class.__name__))


    def add_arguments(self, parser):

        parser.add_argument(
            '-p', '--path',
            type=str,
            default='out/',
            help='the experiment output path')

        parser.add_argument(
            '-e', '--experiment',
            type=str,
            default='',
            help='the experiment id')

        parser.add_argument(
            '-a', '--analyses',
            nargs = '+',
            type=str,
            default='',
            help='names of analyses to run')

        return parser

if __name__ == '__main__':
    analyze = Analyze()
    analyze.run_analysis()
