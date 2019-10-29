from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

from lens.analysis.analysis import Analysis

DEFAULT_COLOR = [color/255 for color in [102, 178 , 255]]


class Snapshots(Analysis):
    def __init__(self):
        super(Snapshots, self).__init__(analysis_type='lattice')

    def get_data(self, client, query):

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
            fields = row.get('fields')
            time_dict[time].update({'fields': fields})

        return time_dict

    def analyze(self, experiment_config, time_data, output_dir):

        time_vec = time_data.keys()
        edge_length = experiment_config['edge_length']
        patches_per_edge = experiment_config['patches_per_edge']  # TODO -- save patches_per_edge
        cell_radius = experiment_config['cell_radius']  # TODO -- save cell_radius
        lattice_scaling = patches_per_edge / edge_length

        # define number of snapshots to be plotted
        n_snapshots = 6

        # get the time steps that will be used
        plot_step = int(len(time_vec)/(n_snapshots-1))
        snapshot_times = time_vec[::plot_step]

        # get number of fields
        field_ids = time_data[time_vec[0]]['fields'].keys()
        n_fields = max(len(field_ids),1)

        fig = plt.figure(figsize=(20*n_snapshots, 20*n_fields))
        plt.rcParams.update({'font.size': 36})
        for index, time in enumerate(snapshot_times):
            field_data = time_data[time]['fields']
            agent_data = time_data[time]['agents']

            if field_ids:
                # plot fields
                for field_id in field_ids:
                    ax = fig.add_subplot(1, n_snapshots, index + 1, adjustable='box')
                    ax.title.set_text('time = {}'.format(time))
                    ax.set_xlim([0, patches_per_edge])
                    ax.set_ylim([0, patches_per_edge])
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])

                    plt.imshow(field_data[field_id],
                               extent=[0,patches_per_edge,0,patches_per_edge],
                               interpolation='nearest',
                               cmap='YlGn')
                    self.plot_agents(ax, agent_data, lattice_scaling, cell_radius, DEFAULT_COLOR)
            else:
                ax = fig.add_subplot(1, n_snapshots, index + 1, adjustable='box')
                ax.title.set_text('time = {}'.format(time))
                ax.set_xlim([0, patches_per_edge])
                ax.set_ylim([0, patches_per_edge])
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                self.plot_agents(ax, agent_data, lattice_scaling, cell_radius, DEFAULT_COLOR)

        # plt.subplots_adjust(wspace=0.7, hspace=0.1)
        plt.savefig(output_dir + '/snapshots', bbox_inches='tight')
        plt.close(fig)

    def plot_agents(self, ax, agent_data, lattice_scaling, cell_radius, agent_color):
        for agent_id, agent_data in agent_data.iteritems():
            location = agent_data['location']
            volume = agent_data['volume']
            y = location[0] * lattice_scaling
            x = location[1] * lattice_scaling
            theta = location[2]
            length = self.volume_to_length(volume, cell_radius) * lattice_scaling

            # plot cell as a 2D line
            dx = length * np.sin(theta)
            dy = length * np.cos(theta)
            width = linewidth_from_data_units(cell_radius * 2, ax) * lattice_scaling
            body_width = width * 0.95

            # plot outline
            ax.plot([y - dy / 2, y + dy / 2], [x - dx / 2, x + dx / 2],
                    linewidth=width,
                    color='k',
                    solid_capstyle='round')

            # plot body
            ax.plot([y - dy / 2, y + dy / 2], [x - dx / 2, x + dx / 2],
                    linewidth=body_width,
                    color=agent_color,
                    solid_capstyle='round')


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