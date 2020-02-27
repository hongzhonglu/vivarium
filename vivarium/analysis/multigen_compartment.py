from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
import numpy as np

from vivarium.analysis.analysis import Analysis, get_compartment

class MultigenCompartment(Analysis):
    def __init__(self):
        super(MultigenCompartment, self).__init__(analysis_type='env_with_compartment')

    def get_data(self, client, query, options={}):
        tags = options.get('tags')
        data = {}
        if tags:
            query.update({'type': 'compartment'})
            history_data = client.find(query)
            history_data.sort('time')
            compartment_history = get_compartment(history_data)

            times = compartment_history['time']

            # get history of all tags
            ports = [port for port in list(compartment_history.keys()) if port not in ['time']]
            tags_history = {}
            for tag in tags:
                for port in ports:
                    if tag in compartment_history[port]:
                        tags_history.update({tag: compartment_history[port][tag]})

            data['time'] = times
            data['tags'] = tags_history

        return data

    def analyze(self, experiment_config, data, output_dir):

        phylogeny = experiment_config['phylogeny']
        compartments = data['compartments']

        # get tag ids and range
        tag_ids = []
        for c_id, c_data in compartments.items():
            try:
                for tag_id, counts in c_data['tags'].items():
                    if tag_id not in tag_ids:
                        tag_ids.append(tag_id)
            except:
                print('No tags for multigen analysis. Specify tags with -t tag_id1 tag_id2')
                return

        if phylogeny:
            # find initial agents in phylogeny
            ancestors = list(phylogeny.keys())
            descendents = list(set([daughter
                for daughters in phylogeny.values() for daughter in daughters]))
            agent_ids = np.setdiff1d(ancestors,descendents)
        else:
            # if no phylogeny, all compartments must be initial agents
            agent_ids = np.array(list(compartments.keys()))


        # make figure
        n_rows = len(tag_ids)
        fig = plt.figure(figsize=(8, n_rows * 3))
        grid = plt.GridSpec(n_rows + 1, 1, wspace=0.4, hspace=1.5)

        # one row for each tag
        for plot_idx, tag_id in enumerate(tag_ids):
            ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)
            plt.title(tag_id, y=1.08)
            ax.set_xlabel('time (hrs)')

            for agent_id in agent_ids:
                simulation_ids = [agent_id]

                # plot all descendants
                while simulation_ids:
                    sim_id = simulation_ids[0]
                    simulation_ids.remove(sim_id)
                    daughter_ids = phylogeny.get(sim_id)
                    if daughter_ids:
                        simulation_ids.extend(daughter_ids)  # add agent's daughters
                    if compartments.get(sim_id):
                        plot_single(ax, compartments[sim_id], tag_id)


        plt.savefig(output_dir + '/multigen')
        plt.close(fig)

def plot_single(ax, data, tag_id):
    series = data['tags'][tag_id]
    time_data = data['time']
    time_vec = [t / 3600 for t in time_data]  # convert to hours
    ax.plot(time_vec, series)
