from __future__ import absolute_import, division, print_function

import json

from vivarium.analysis.analysis import Analysis

class Topology(Analysis):
    def __init__(self):
        super(Topology, self).__init__(analysis_type='experiment')

    def get_data(self, client, query, options={}):
        return

    def analyze(self, experiment_config, history_data, output_dir):
        topology = experiment_config.get('topology')

        # save topology
        json_out = json.dumps(topology)
        f = open(output_dir + '/topology.json', 'w')
        f.write(json_out)
        f.close()