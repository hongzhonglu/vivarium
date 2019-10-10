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

def get_compartment(data):
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

def get_lattice(data):
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

def get_experiment(data):
    experiment_config = {}
    for row in data:
        if row.get('type') == 'lattice':
            experiment_config['edge_length'] = row['edge_length']
        elif row.get('topology'):
            experiment_config['topology'] = row['topology']

    return experiment_config

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

        # get the experiment configuration
        query = {'experiment_id': experiment_id}
        config_data = self.client.simulations.configuration.find(query)
        exp_data = get_experiment(config_data)

        # get the time series data for compartments
        query.update({'type': 'compartment'})
        history_data = self.client.simulations.history.find(query)
        history_data.sort('time')
        compartment_hist = get_compartment(history_data)   # ['internal', 'external', 'time']

        # analyze all simulations in this experiment
        simulation_ids = get_sims_from_exp(self.client.simulations.history, experiment_id)
        for simulation_id in simulation_ids:
            compartment_analysis(compartment_hist, experiment_id, simulation_id, path)

        # run lattice analysis
        query.update({'type': 'lattice'})
        history_data = self.client.simulations.history.find(query)
        history_data.sort('time')
        lattice_hist = get_lattice(history_data)
        location_trace(exp_data, lattice_hist, experiment_id, path)  # TODO -- run this with a separate option?



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


        return parser

if __name__ == '__main__':
    command = Analyze()
