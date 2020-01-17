from __future__ import absolute_import, division, print_function

import json

from vivarium.analysis.analysis import Analysis

class Topology(Analysis):
    def __init__(self):
        super(Topology, self).__init__(analysis_type='experiment')

    def get_data(self, client, query, options={}):
        return

    def analyze(self, experiment_config, history_data, output_dir):
        agents = experiment_config.get('agents')

        topologies = {}
        for agent_id, specs in agents.items():
            name = specs.get('name')
            topology = specs.get('topology')

            topologies[agent_id] = {
                'name': name,
                'topology': topology}

        # save topology
        json_out = json.dumps(topologies)
        f = open(output_dir + '/topology.json', 'w')
        f.write(json_out)
        f.close()