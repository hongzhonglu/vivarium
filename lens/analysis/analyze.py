from __future__ import absolute_import, division, print_function

import os
import argparse
from pymongo import MongoClient

import matplotlib

from lens.analysis.analyze_compartment import compartment_analysis
from lens.analysis.analyze_lattice import location_trace

# mongoDB local url
url='localhost:27017'

# helper functions
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

def db_to_dict_compartment(data):
    skip_keys = [
        'type',
        'simulation_id',
        'experiment_id',
        'time',
        '_id']

    # organize data into a dict
    data_dict = {}
    time_vec = []
    for row in data:
        time = row.get('time')
        time_vec.append(time)
        data_keys = [key for key in row.keys() if key not in skip_keys]
        for key in data_keys:
            value = row[key]
            if key not in data_dict.keys():
                data_dict[key] = {}
            for mol_id, v in value.iteritems():
                if mol_id in data_dict[key].keys():
                    data_dict[key][mol_id].append(v)
                else:
                    data_dict[key][mol_id] = [v]
    data_dict['time'] = time_vec
    return data_dict

def db_to_dict_lattice(data):
    skip_keys = [
        'type',
        'agent_id',
        'simulation_id',
        'experiment_id',
        'time',
        '_id']

    # organize data into a dict
    data_dict = {}
    time_vec = []
    for row in data:
        # get time
        time = row.get('time')
        time_vec.append(time)

        # get agent
        agent_id = row['agent_id']
        if agent_id not in data_dict.keys():
            data_dict[agent_id] = {}

        # get data
        data_keys = [key for key in row.keys() if key not in skip_keys]
        for key in data_keys:
            value = row[key]
            if key in data_dict[agent_id].keys():
                data_dict[agent_id][key].append(value)
            else:
                data_dict[agent_id][key] = [value]

    data_dict['time'] = time_vec
    return data_dict


class Analyze(object):
    '''
    example use:
        python lens/analysis/analyze.py -e lattice -s 2
    '''

    def __init__(self):

        parser = argparse.ArgumentParser(description='analyze rate laws')
        parser = self.add_arguments(parser)
        args = parser.parse_args()

        self.client = MongoClient(url)

        path = args.path
        experiment_id = args.experiment
        simulation_id = args.simulation
        analysis_type = args.type

        # query database
        query = {'type': analysis_type}
        if experiment_id:
            query.update({'experiment_id': experiment_id})
        if simulation_id:
            query.update({'simulation_id': simulation_id})

        data = self.client.simulations.output.find(query)
        data.sort('time')

        if analysis_type == 'compartment':
            # structure data into dictionary for single analysis
            data_dict = db_to_dict_compartment(data)

            if not simulation_id:
                # simulation_id not specified, analyze all simulations in this experiment
                simulation_ids = get_sims_from_exp(self.client.simulations.output, experiment_id)
                for simulation_id in simulation_ids:
                    compartment_analysis(data_dict, experiment_id, simulation_id, path)
            else:
                compartment_analysis(data_dict, experiment_id, simulation_id, path)

        elif analysis_type == 'lattice':
            # structure data into dictionary for lattice analysis
            data_dict = db_to_dict_lattice(data)
            location_trace(data_dict, experiment_id, path)  # TODO -- run this with a separate option?


    def add_arguments(self, parser):

        parser.add_argument(
            '-p', '--path',
            type=str,
            default='out/',
            help='the experiment output path')

        parser.add_argument(
            '-t', '--type',
            type=str,
            default='compartment',
            help='the type of analysis')

        parser.add_argument(
            '-e', '--experiment',
            type=str,
            default='',
            help='the experiment id')

        parser.add_argument(
            '-s', '--simulation',
            type=str,
            default='',
            help='the simulation/agent id. leave empty to analyze all simulations')

        return parser

if __name__ == '__main__':
    command = Analyze()
