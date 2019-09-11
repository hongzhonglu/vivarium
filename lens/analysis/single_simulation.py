from __future__ import absolute_import, division, print_function

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pymongo import MongoClient


url='localhost:27017'

class AnalyzeState(object):
    '''
    example use:
        python lens/analysis/single_simulation.py -e lattice -s 2
    '''

    def __init__(self):

        parser = argparse.ArgumentParser(description='analyze rate laws')
        parser = self.add_arguments(parser)
        args = parser.parse_args()

        self.client = MongoClient(url)

        experiment_id = args.experiment
        simulation_id = args.simulation
        self.run_analysis(experiment_id, simulation_id)

    def run_analysis(self, experiment_id, simulation_id):

        # all_output = [row for row in self.client.simulations.output.find()]
        query = {'experiment_id': experiment_id,
                 'simulation_id': simulation_id}
        data = self.client.simulations.output.find(query) # TODO - order by time.
        data.sort('time')

        # organize data
        data_structured = {
            'internal': {},
            'external': {}}
        time_vec = []
        for row in data:
            internal = row['internal']
            external = row['external']
            time = row['time']
            time_vec.append(time)

            for mol_id, value in internal.iteritems():
                if mol_id in data_structured['internal'].keys():
                    data_structured['internal'][mol_id].append(value)
                else:
                    data_structured['internal'][mol_id] = [value]
            for mol_id, value in external.iteritems():
                if mol_id in data_structured['external'].keys():
                    data_structured['external'][mol_id].append(value)
                else:
                    data_structured['external'][mol_id] = [value]

        fig = plt.figure(figsize=(8,25))
        plot_idx = 1
        for mol_id, series in data_structured['external'].iteritems():

            ax = fig.add_subplot(30, 1, plot_idx)
            ax.plot(time_vec, series)
            ax.title.set_text(mol_id)
            plot_idx += 1

        plt.tight_layout()
        fig_path = 'out/' + experiment_id + '/' + simulation_id + '/single'
        plt.savefig(fig_path)
        plt.clf()


    def add_arguments(self, parser):

        parser.add_argument(
            '-e', '--experiment',
            type=str,
            default='',
            help='the experiment id')

        parser.add_argument(
            '-s', '--simulation',
            type=str,
            default='',
            help='the simulation id')

        return parser

if __name__ == '__main__':
    command = AnalyzeState()
