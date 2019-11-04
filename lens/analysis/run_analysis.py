from __future__ import absolute_import, division, print_function

import os
import argparse
from pymongo import MongoClient

from lens.analysis.analyze_compartment import Compartment
from lens.analysis.multigen_compartment import MultigenCompartment
from lens.analysis.analyze_lattice import LatticeTrace
from lens.analysis.chemotaxis_trace import ChemotaxisTrace
from lens.analysis.snapshots import Snapshots


# classes for run_analysis to cycle through
analysis_classes = {
    # 'chemotaxis': ChemotaxisTrace,
    'compartment': Compartment,
    'multigen': MultigenCompartment,
    'lattice': LatticeTrace,
    'snapshots': Snapshots,
}

# mongoDB local url
url='localhost:27017'


class AnalysisError(Exception):
    pass

def get_phylogeny(data):
    # given the phylogeny client and experiment id, return a dict of {parent_id: [daughter1_id, daughter2_id]}
    phylogeny_data = {}
    for row in data:
        simulation_id = row.get('simulation_id')
        daughters = row.get('daughters')
        phylogeny_data[simulation_id] = daughters

    return phylogeny_data

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
        self.tags = args.tags

    def run_analysis(self):
        # get the tables
        config_client = self.client.simulations.configuration
        phylogeny_client = self.client.simulations.phylogeny
        history_client = self.client.simulations.history

        # get simulations ids
        simulation_ids = get_sims_from_exp(history_client, self.experiment_id)

        # make plot output directories, for experiment and for each simulation
        output_dir = os.path.join(self.path, self.experiment_id)
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
        if not experiment_config:
            raise AnalysisError('database has no experiment id: {}'.format(self.experiment_id))
        active_processes = experiment_config['topology'].keys()

        # get the phylogenetic tree
        query = {'experiment_id': self.experiment_id}
        phylogeny_data = phylogeny_client.find(query)
        experiment_config['phylogeny'] = get_phylogeny(phylogeny_data)

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

                # run the compartment analysis for each simulation in simulation_ids
                if analysis.analysis_type is 'compartment':
                    # A compartment analysis is run on a single compartment.
                    # It expects to run queries on the compartment tables in the DB.
                    # Output is saved to the compartment's directory.

                    for sim_id in simulation_ids:
                        compartment_query = query.copy()
                        compartment_query.update({'simulation_id': sim_id})
                        data = analysis.get_data(history_client, compartment_query)

                        sim_out_dir = os.path.join(output_dir, sim_id)
                        analysis.analyze(experiment_config, data, sim_out_dir)

                elif analysis.analysis_type is 'environment':
                    # A environment analysis is run on the environment.
                    # It expects to run queries on the environment (lattice) tables in the DB.
                    # Output is saved to the experiment's base directory.

                    environment_data = analysis.get_data(history_client, query.copy())
                    analysis.analyze(experiment_config, environment_data, output_dir)

                elif analysis.analysis_type is 'both':
                    # A both analysis is run on the environment AND compartments.
                    # It expects to run queries on the environment (lattice) tables in the DB,
                    # but if an option is passed with simulation id, it will query those as well.
                    # Output is saved to the experiment's base directory.

                    if self.tags:
                        tag_data = {}
                        for sim_id in simulation_ids:
                            compartment_query = query.copy()
                            compartment_query.update({'simulation_id': sim_id})
                            options = {'tags': self.tags}
                            tag_data[sim_id] = analysis.get_data(history_client, compartment_query, options)

                    environment_data = analysis.get_data(history_client, query.copy())

                    data = {
                        'tags': tag_data,
                        'environment': environment_data}

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

        parser.add_argument(
            '-t', '--tags',
            nargs='+',
            type=str,
            default='',
            help='names of molecules to tag')

        return parser

if __name__ == '__main__':
    analyze = Analyze()
    analyze.run_analysis()
