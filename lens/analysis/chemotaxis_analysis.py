from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis, get_lattice, get_compartment

class Chemotaxis(Analysis):
    def __init__(self):
        super(Chemotaxis, self).__init__(analysis_type='both')

    def requirements(self):
        return ['motor']

    def get_data(self, client, query, options={}):

        type = options.get('type')
        if type is 'compartment':
            return
        elif type is 'environment':
            # get data about agent locations ('type': 'lattice')
            query_lattice = {'type': 'lattice'}
            query_lattice.update(query)
            data_lattice = client.find(query_lattice)
            data_lattice.sort('time')

            # get data about concentrations ('type': 'lattice-field')
            query_field = {'type': 'lattice-field'}
            query_field.update(query)
            data_field = client.find(query_field)
            data_field.sort('time')

            # organize data into a dict by time
            time_dict = {}
            for row in data_lattice:
                time = row.get('time')
                if time not in time_dict:
                    time_dict[time] = {}
                    time_dict[time]['agents'] = {}

                agent_id = row['agent_id']
                location = row['location']
                volume = row['volume']

                time_dict[time]['agents'][agent_id] = {
                    'location': location,
                    'volume': volume,
                }

            # add fields
            for row in data_field:
                time = row.get('time')
                if time not in time_dict:
                    break

                fields = row.get('fields', [])
                time_dict[time].update({'fields': fields})

            return time_dict

    def analyze(self, experiment_config, data, output_dir):

        env_data = data['environment']
        compartments_data = data['compartments']

        time_vec = env_data.keys()
        agent_ids = compartments_data.keys()
        phylogeny = experiment_config['phylogeny']
        edge_length_x = experiment_config['edge_length_x']
        edge_length_y = experiment_config['edge_length_y']



        import ipdb;
        ipdb.set_trace()

        #
        #
        # plt.figure(figsize=(8, 8))
        #
        #
        #
        #
        # plt.savefig(output_dir + '/chemotaxis_analysis')
        # plt.clf()