from __future__ import absolute_import, division, print_function

import os
import argparse
from pymongo import MongoClient

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from lens.utils.query import db_to_dict, get_sims_from_exp

# mongoDB local url
url='localhost:27017'


class AnalyzeState(object):
    '''
    example use:
        python lens/analysis/analyze_experiment.py -e lattice -s 2
    '''

    def __init__(self):

        parser = argparse.ArgumentParser(description='analyze rate laws')
        parser = self.add_arguments(parser)
        args = parser.parse_args()

        self.client = MongoClient(url)

        path = args.path
        experiment_id = args.experiment
        simulation_id = args.simulation

        if not simulation_id:
            simulation_ids = get_sims_from_exp(self.client.simulations.output, experiment_id)
            for simulation_id in simulation_ids:
                self.run_analysis(experiment_id, simulation_id, path)
        else:
            self.run_analysis(experiment_id, simulation_id, path)

    def run_analysis(self, experiment_id, simulation_id, output_dir='out'):

        query = {'experiment_id': experiment_id,
                 'simulation_id': simulation_id}
        data = self.client.simulations.output.find(query)
        data.sort('time')
        data_dict = db_to_dict(data)
        data_keys = [key for key in data_dict.keys() if key is not 'time']

        time_vec = [t / 3600 for t in data_dict['time']]  # convert to hours
        n_data = [len(data_dict[key].keys()) for key in data_keys]
        n_rows = sum(n_data)

        fig = plt.figure(figsize=(8,n_rows*1.5))
        plot_idx = 1

        for key in data_keys:

            for mol_id, series in data_dict[key].iteritems():
                ax = fig.add_subplot(n_rows, 1, plot_idx)
                ax.plot(time_vec, series)
                ax.title.set_text(str(key) + ': ' + mol_id)
                # ax.ticklabel_format(useOffset=False)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                ax.set_xlabel('time (hrs)')
                plot_idx += 1

        # make figure output directory and save figure
        fig_path = os.path.join(output_dir, experiment_id, simulation_id)
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)
        plt.tight_layout()
        plt.savefig(fig_path + '/analyze_compartment')
        plt.clf()


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
            '-s', '--simulation',
            type=str,
            default='',
            help='the simulation/agent id. leave empty to analyze all simulations')

        return parser

if __name__ == '__main__':
    command = AnalyzeState()
