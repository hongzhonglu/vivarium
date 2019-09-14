'''
Analysis plot

run by passing path to the output directory:
> python lens/analysis/location_trace.py -p 'out/out_dir'
'''
from __future__ import absolute_import

import os
import argparse

import matplotlib
matplotlib.use('TkAgg')  # solves "RuntimeError: Python is not installed as a framework."
import matplotlib.pyplot as plt

import lens.utils.filepath as fp

# TODO -- this analysis script needs to be rehabbed without TableWriter
class TraceLocationPlot(object):

    def __init__(self):
        parser = argparse.ArgumentParser(description='analysis plot')
        parser = self.add_arguments(parser)
        args = parser.parse_args()
        path = args.path

        self.plot(path)

    def add_arguments(self, parser):
        parser.add_argument(
            '-p', '--path',
            type=str,
            help='run analysis on the experiment specified by path')

        return parser

    def plot(self, path):

        # TODO -- replace with database readout
        # lattice_reader = TableReader(os.path.join(path))
        # edge_length = lattice_reader.readAttribute("edge_length")
        #
        edge_length = 0
        agent_locations = {}
        #
        # # get all agent_ids in this directory and read their locations if available.
        # agent_ids = os.listdir(path)
        # for agent_id in agent_ids:
        #     agent_path = os.path.join(path, agent_id, 'location')
        #     # if agent's location has been saved to a table, read it and add to locations
        #     if os.path.isfile(agent_path):
        #         agent_reader = TableReader(os.path.join(path, agent_id))
        #         start_time = agent_reader.readAttribute("start_time")
        #         location = agent_reader.readColumn("location")
        #         agent_locations[agent_id] = location

        plt.figure(figsize=(8, 8))
        for agent_id, locations in agent_locations.iteritems():

            x_coord = locations[:, 0]
            y_coord = locations[:, 1]
            angle = locations[:, 2]

            for i in range(len(locations)):
                plt.plot(x_coord[i:i + 2], y_coord[i:i + 2], 'b-')

        plt.xlim((0, edge_length))
        plt.ylim((0, edge_length))

        # save figure
        output_dir = fp.makedirs(path, 'plot_out')
        plt.savefig(os.path.join(output_dir, 'location_trace'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'location_trace.pdf'), bbox_inches='tight')

if __name__ == '__main__':
    command = TraceLocationPlot()