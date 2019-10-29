from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from lens.analysis.analysis import Analysis, get_compartment

class Compartment(Analysis):
    def __init__(self):
        super(Compartment, self).__init__(analysis_type='compartment')

    def get_data(self, client, query):
        query.update({'type': 'compartment'})
        history_data = client.find(query)
        history_data.sort('time')
        compartment_history = get_compartment(history_data)

        return compartment_history

    def analyze(self, experiment_config, data_dict, output_dir):
        # remove series with all zeros
        zero_state = []
        for key1 in data_dict.iterkeys():
            if key1 is not 'time':
                for key2, series in data_dict[key1].iteritems():
                    if all(v == 0 for v in series):
                        zero_state.append((key1, key2))
        for (key1, key2) in zero_state:
            del data_dict[key1][key2]

        data_keys = [key for key in data_dict.keys() if key is not 'time']
        time_vec = [t / 3600 for t in data_dict['time']]  # convert to hours

        # make figure, with grid for subplots
        n_data = [len(data_dict[key].keys()) for key in data_keys]
        n_zeros = len(zero_state)
        n_rows = sum(n_data) + int(n_zeros/20)  # 20 zero_state ids per additional subplot

        fig = plt.figure(figsize=(8, n_rows * 2.5))
        grid = plt.GridSpec(n_rows+1, 1, wspace=0.4, hspace=1.5)

        # plot data
        plot_idx = 0
        for key in data_keys:
            for mol_id, series in sorted(data_dict[key].iteritems()):
                ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)
                ax.plot(time_vec, series)
                ax.title.set_text(str(key) + ': ' + mol_id)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                ax.set_xlabel('time (hrs)')
                plot_idx += 1

        # additional data as text
        if zero_state:
            zeros = ['{}[{}]'.format(state, role) for (role, state) in zero_state]
            ax = fig.add_subplot(grid[plot_idx:, 0])
            ax.text(0.02, 0.1, 'states with all zeros: {}'.format(zeros), wrap=True)
            ax.axis('off')

        plt.savefig(output_dir + '/compartment')
        plt.close(fig)
