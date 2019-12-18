from __future__ import absolute_import, division, print_function

import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from lens.analysis.analysis import Analysis, get_compartment


def set_axes(ax, show_xaxis=False):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)

class Compartment(Analysis):
    def __init__(self):
        super(Compartment, self).__init__(analysis_type='compartment')

    def get_data(self, client, query, options={}):
        sim_id = query['simulation_id']
        query.update({'type': 'compartment'})
        history_data = client.find(query)
        history_data.sort('time')
        compartment_history = get_compartment(history_data)

        # put sim_id back into data
        compartment_history['sim_id'] = sim_id

        return compartment_history

    def analyze(self, experiment_config, data_dict, output_dir):

        skip_keys = ['time', 'sim_id']
        sim_id = data_dict['sim_id']

        # remove series with all zeros
        zero_state = []
        for key1 in data_dict.iterkeys():
            if key1 not in skip_keys:
                for key2, series in data_dict[key1].items():
                    if all(v == 0 for v in series):
                        zero_state.append((key1, key2))
        for (key1, key2) in zero_state:
            del data_dict[key1][key2]

        data_keys = [key for key in data_dict.keys() if key not in skip_keys]
        time_vec = [t / 3600 for t in data_dict['time']]  # convert to hours

        # make figure, with grid for subplots
        n_rows_base = 15
        n_data = [len(data_dict[key].keys()) for key in data_keys]
        n_zeros = len(zero_state)
        n_cols = int(math.ceil(sum(n_data)/float(n_rows_base)))
        n_rows = n_rows_base + int(math.ceil(n_zeros/20.0))  # zero_state ids per additional subplot

        fig = plt.figure(figsize=(n_cols*6, n_rows*2.5))
        fig.suptitle('{}'.format(sim_id), fontsize=12)
        plt.rcParams.update({'font.size': 12})
        grid = plt.GridSpec(n_rows+1, n_cols, wspace=0.3, hspace=1.2)

        # plot data
        row_idx = 0
        col_idx = 0
        for key in data_keys:
            for mol_id, series in sorted(data_dict[key].items()):
                ax = fig.add_subplot(grid[row_idx, col_idx])
                ax.plot(time_vec, series)
                ax.title.set_text(str(key) + ': ' + mol_id)

                row_idx += 1
                if row_idx > n_rows_base:
                    set_axes(ax, True)
                    ax.set_xlabel('time (hrs)')
                    row_idx = 0
                    col_idx += 1
                else:
                    set_axes(ax)

        # additional data as text
        if zero_state:
            zeros = ['{}[{}]'.format(state, role) for (role, state) in zero_state]
            ax = fig.add_subplot(grid[n_rows_base:, :])
            ax.text(0.02, 0.1, 'states with all zeros: {}'.format(zeros), wrap=True)
            ax.axis('off')

        plt.savefig(output_dir + '/compartment')
        plt.close(fig)
