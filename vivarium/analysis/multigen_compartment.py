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
            tags_history = {tag: compartment_history['cell'][tag] for tag in tags}

            data['time'] = times
            data['tags'] = tags_history

        return data

    def analyze(self, experiment_config, data, output_dir):

        phylogeny = experiment_config['phylogeny']
        compartments = data['compartments']

        if not compartments:
            print('no tags for multigen_compartment analysis. specify tag with "-t tag_id"')
            return

        if phylogeny:
            # find initial agents in phylogeny
            ancestors = list(phylogeny.keys())
            descendents = list(set([daughter
                for daughters in phylogeny.values() for daughter in daughters]))
            initial_agents = np.setdiff1d(ancestors,descendents)
        else:
            # if no phylogeny, all compartments must be initial agents
            initial_agents = np.array(list(compartments.keys()))

        n_rows = len(initial_agents)  # 20 zero_state ids per additional subplot
        fig = plt.figure(figsize=(8, n_rows * 3))
        grid = plt.GridSpec(n_rows + 1, 1, wspace=0.4, hspace=1.5)
        for plot_idx, agent_id in enumerate(initial_agents):
            ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)
            simulation_ids = [agent_id]

            # plot all descendants of ancestor
            while simulation_ids:
                sim_id = simulation_ids[0]
                simulation_ids.remove(sim_id)
                daughter_ids = phylogeny.get(sim_id)
                if daughter_ids:
                    simulation_ids.extend(daughter_ids)
                if compartments.get(sim_id):
                    plot_single(ax, compartments[sim_id])

            ax.set_xlabel('time (hrs)')

        plt.savefig(output_dir + '/multigen')
        plt.close(fig)

def plot_single(ax, data):
    tag_data = data['tags']
    time_data = data['time']
    time_vec = [t / 3600 for t in time_data]  # convert to hours
    for mol_id, series in tag_data.items():
        ax.plot(time_vec, series)
