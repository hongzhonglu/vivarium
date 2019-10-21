from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis

class Snapshots(Analysis):
    def __init__(self):
        super(Snapshots, self).__init__(analysis_type='lattice')

    def get_data(self, client):
        query = {'type': 'lattice'}
        data = client.find(query)
        data.sort('time')

        # organize data into a dict by time
        time_dict = {}
        for row in data:
            time = row.get('time')
            if time not in time_dict:
                time_dict[time] = {}

            agent_id = row['agent_id']
            location = row['location']
            volume = row['volume']

            time_dict[time][agent_id] = {
                'location': location,
                'volume': volume,
            }

        return time_dict

    def analyze(self, experiment_config, time_dict, output_dir):

        time_vec = time_dict.keys()
        edge_length = experiment_config['edge_length']
        patches_per_edge = experiment_config['patches_per_edge']  # TODO -- save patches_per_edge
        cell_radius = experiment_config['cell_radius']  # TODO -- save cell_radius

        # define number of snapshots to be plotted
        n_snapshots = 10

        # get the time steps that will be used
        plot_step = int(len(time_vec)/n_snapshots)
        snapshot_times = time_vec[::plot_step]

        fig = plt.figure(figsize=(20*n_snapshots, 20))
        for index, time in enumerate(snapshot_times):
            ax = fig.add_subplot(1, n_snapshots, index+1)
            ax.set_xlim([0, edge_length])
            ax.set_ylim([0, edge_length])

            snapshot_data = time_dict[time]
            for agent_id, agent_data in snapshot_data.iteritems():
                location = agent_data['location']
                volume = agent_data['volume']
                y = location[0]
                x = location[1]
                theta = location[2]
                length = self.volume_to_length(volume, cell_radius)

                # plot cell as a 2D line
                dx = length * np.sin(theta)
                dy = length * np.cos(theta)
                width = linewidth_from_data_units(cell_radius*2, ax)

                ax.plot([x-dx/2, x+dx/2], [y-dy/2, y+dy/2], linewidth=width, color='slateblue', solid_capstyle='round')

        # plt.subplots_adjust(wspace=0.7, hspace=0.1)
        plt.savefig(output_dir + '/snapshots', bbox_inches='tight')
        plt.clf()


    def volume_to_length(self, volume, cell_radius):
        '''
        get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
        a=length of cylinder without rounded caps, l=total length:
        V = (4/3)*PI*r^3 + PI*r^2*a
        l = a + 2*r
        '''
        pi = np.pi

        cylinder_length = (volume - (4 / 3) * pi * cell_radius ** 3) / (pi * cell_radius ** 2)
        total_length = cylinder_length + 2 * cell_radius

        return total_length

def linewidth_from_data_units(linewidth, axis, reference='y'):
    """
    Convert a linewidth in data units to linewidth in points.

    Parameters
    ----------
    linewidth: float
        Linewidth in data units of the respective reference-axis
    axis: matplotlib axis
        The axis which is used to extract the relevant transformation
        data (data limits and size must not change afterwards)
    reference: string
        The axis that is taken as a reference for the data width.
        Possible values: 'x' and 'y'. Defaults to 'y'.

    Returns
    -------
    linewidth: float
        Linewidth in points
    """
    fig = axis.get_figure()
    if reference == 'x':
        length = fig.bbox_inches.width * axis.get_position().width
        value_range = np.diff(axis.get_xlim())
    elif reference == 'y':
        length = fig.bbox_inches.height * axis.get_position().height
        value_range = np.diff(axis.get_ylim())
    # Convert length to points
    length *= 72
    # Scale linewidth to value range
    return linewidth * (length / value_range)
