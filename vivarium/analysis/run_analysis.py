from __future__ import absolute_import, division, print_function

import os
import argparse
from pymongo import MongoClient

from vivarium.analysis.analyze_compartment import Compartment
from vivarium.analysis.multigen_compartment import MultigenCompartment
from vivarium.analysis.location_trace import LatticeTrace
from vivarium.analysis.motor import Motor
from vivarium.analysis.snapshots import Snapshots
from vivarium.analysis.topology import Topology


# classes for run_analysis to cycle through
analysis_classes = {
    'motor': Motor,
    'compartment': Compartment,
    'multigen': MultigenCompartment,
    'location': LatticeTrace,
    'snapshots': Snapshots,
    'topology': Topology,
}

class AnalysisError(Exception):
    pass

def get_phylogeny(client, experiment_id):
    # given the data from the phylogeny table,
    # return a dict with {parent_id: [daughter1_id, daughter2_id]}
    query = {'experiment_id': experiment_id}
    data = client.find(query)
    phylogeny_data = {}
    for row in data:
        simulation_id = row.get('simulation_id')
        daughters = row.get('daughters')
        phylogeny_data[simulation_id] = daughters

    return phylogeny_data

def get_experiment(client, experiment_id):
    query = {'experiment_id': experiment_id}
    data = client.find(query)
    experiment_config = {'agents': {}}
    for row in data:
        if row.get('type') == 'lattice':
            experiment_config['edge_length_x'] = row['edge_length_x']
            experiment_config['edge_length_y'] = row['edge_length_y']
            experiment_config['patches_per_edge_x'] = row['patches_per_edge_x']
            experiment_config['patches_per_edge_y'] = row['patches_per_edge_y']
            experiment_config['cell_radius'] = row['cell_radius']
            experiment_config['description'] = row['description']

        elif row.get('type') == 'compartment':
            sim_id = row.get('simulation_id')
            experiment_config['agents'][sim_id] = {}
            experiment_config['agents'][sim_id]['topology'] = row.get('topology')
            experiment_config['agents'][sim_id]['name'] = row.get('name')

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

def compartment_query(query, analysis, history_client, sim_id, options={}):
    compartment_query = query.copy()
    compartment_query.update({'simulation_id': sim_id})
    data = analysis.get_data(history_client, compartment_query, options)
    return data

def environment_query(query, analysis, history_client, options={}):
    data = analysis.get_data(history_client, query.copy(), options)
    return data


class Analyze(object):
    '''
    example use:
        python vivarium/analysis/run_analysis.py -e lattice
    '''
    client = None

    def __init__(self):

        parser = argparse.ArgumentParser(description='analyze rate laws')
        parser = self.add_arguments(parser)
        args = parser.parse_args()

        # create singleton instance of mongo client
        if Analyze.client is None:
            Analyze.client = MongoClient(args.mongo_host)

        self.path = args.path
        self.experiment_id = args.experiment
        self.analyses = args.analyses
        self.tags = args.tags

    def run_analysis(self):
        # get the tables
        config_client = self.client.simulations.configuration
        phylogeny_client = self.client.simulations.phylogeny
        history_client = self.client.simulations.history

        # baseline query for experiment_id
        query = {'experiment_id': self.experiment_id}

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
        experiment_config = get_experiment(config_client, self.experiment_id)
        if not experiment_config:
            raise AnalysisError('database has no experiment id: {}'.format(self.experiment_id))

        # get all processes active in the agents.
        active_processes = set()
        for agent, specs in experiment_config['agents'].items():
            processes = list(specs['topology'].keys())
            active_processes.update(processes)

        # get the phylogenetic tree in experiment config
        experiment_config['phylogeny'] = get_phylogeny(phylogeny_client, self.experiment_id)

        # get list of analysis objects to run from self.analyses, or run all analyses
        if self.analyses:
            run_analyses = [analysis_classes[analysis_id] for analysis_id in self.analyses]
        else:
            run_analyses = analysis_classes.values()

        # run analyses that declare the active processes
        for analysis_class in run_analyses:
            analysis = analysis_class()
            required_processes = analysis.requirements()

            if all(processes in active_processes for processes in required_processes):

                if analysis.analysis_type is 'experiment':
                    data = analysis.get_data(history_client, query.copy())
                    analysis.analyze(experiment_config, data, output_dir)

                elif analysis.analysis_type is 'compartment':
                    # Run a compartment analysis for each simulation in simulation_ids
                    # A 'compartment' analysis is run on a single compartment.
                    # It expects to run queries on the compartment tables in the DB.
                    # Output is saved to the compartment's directory.

                    for sim_id in simulation_ids:
                        data = compartment_query(query, analysis, history_client, sim_id)
                        sim_out_dir = os.path.join(output_dir, sim_id)
                        analysis.analyze(experiment_config, data, sim_out_dir)

                elif analysis.analysis_type is 'environment':
                    # A 'environment' analysis is run on the environment.
                    # It expects to run queries on the environment (lattice) tables in the DB.
                    # Output is saved to the experiment's base directory.

                    data = environment_query(query, analysis, history_client)
                    analysis.analyze(experiment_config, data, output_dir)

                elif analysis.analysis_type is 'env_with_compartment':
                    # An 'env_with_compartment' analysis is run on the environment AND compartments,
                    # from the point of view of the environment.
                    # It expects to run queries on the environment (lattice) tables in the DB,
                    # but if an option is passed with simulation id, it will query those as well.
                    # A single output is saved to the experiment's base directory.

                    compartment_data = {}
                    for sim_id in simulation_ids:
                        options = {
                            'type': 'compartment',
                            'tags': self.tags}
                        compartment_data[sim_id] = compartment_query(query, analysis, history_client, sim_id, options)

                    options = {'type': 'environment'}
                    environment_data = environment_query(query, analysis, history_client, options)

                    data = {
                        'compartments': compartment_data,
                        'environment': environment_data}

                    analysis.analyze(experiment_config, data, output_dir)

                elif analysis.analysis_type is 'compartment_with_env':
                    # An 'compartment_with_env' analysis is run on the environment AND compartments,
                    # from the point of view of a compartment.
                    # It expects to run queries on the environment (lattice) tables in the DB,
                    # but if an option is passed with simulation id, it will query those as well.
                    # Output is saved to each compartment's directory.

                    for sim_id in simulation_ids:
                        sim_processes = list(experiment_config['agents'][sim_id]['topology'].keys())
                        required = analysis.requirements()

                        # if THIS sim has the required processes
                        if all(processes in sim_processes for processes in required):
                            env_query = query.copy()
                            env_query.update({'agent_id': sim_id})
                            options = {'type': 'environment'}
                            environment_data = environment_query(env_query, analysis, history_client, options)
                            data = {'environment': environment_data}

                            options = {
                                'type': 'compartment',
                                'tags': self.tags}
                            compartment_data = compartment_query(query, analysis, history_client, sim_id, options)
                            data['compartment'] = compartment_data

                            sim_out_dir = os.path.join(output_dir, sim_id)

                            analysis.analyze(experiment_config, data, sim_out_dir)

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
            '-m', '--mongo-host',
            type=str,
            default='localhost:27017')

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
