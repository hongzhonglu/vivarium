from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from lens.analysis.analysis import Analysis, get_compartment
from lens.actor.process import dict_merge

# DEFAULT_COLOR = [color/255 for color in [102, 178, 255]]
DEFAULT_COLOR = [210/360, 100.0/100.0, 0.0/100.0]  # HSV

class Snapshots(Analysis):
    def __init__(self):
        super(Snapshots, self).__init__(analysis_type='both')

    def get_data(self, client, query, options={}):

        tags = options.get('tags')
        if tags:
            query.update({'type': 'compartment'})
            history_data = client.find(query)
            history_data.sort('time')
            compartment_history = get_compartment(history_data)

            times = compartment_history['time']
            tags_history = {tag: compartment_history['cell'][tag] for tag in tags}

            # arrange data by time, for easy integration with environment data
            time_dict = {time: {} for time in times}
            for tag, series in tags_history.iteritems():
                tag_hist = {time: {'tags': {tag: state}} for time, state in zip(times,series)}
                time_dict = dict_merge(dict(time_dict), tag_hist)

            return time_dict

        else:
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
                    # time_dict[time] = {}

                fields = row.get('fields', [])
                time_dict[time].update({'fields': fields})

            return time_dict

    def analyze(self, experiment_config, data, output_dir):

        time_data = data['environment']
        tags_data = data['tags']

        time_vec = time_data.keys()
        edge_length = experiment_config['edge_length']
        patches_per_edge = experiment_config['patches_per_edge']
        cell_radius = experiment_config['cell_radius']
        lattice_scaling = patches_per_edge / edge_length

        # define number of snapshots to be plotted
        n_snapshots = 8

        # get the time steps that will be used
        plot_step = int(len(time_vec)/(n_snapshots-1))
        snapshot_times = time_vec[::plot_step]

        # get number of fields
        field_ids = time_data[time_vec[0]].get('fields',{}).keys()
        n_fields = max(len(field_ids),1)

        fig = plt.figure(figsize=(20*n_snapshots, 20*n_fields))
        plt.rcParams.update({'font.size': 36})
        for index, time in enumerate(snapshot_times):
            field_data = time_data[time].get('fields')
            agent_data = time_data[time]['agents']

            agent_tags = {agent_id: tags_data[agent_id].get(time, {}) for agent_id in agent_data.keys()}
            agent_data = dict_merge(dict(agent_data), agent_tags)

            if field_ids:
                # plot fields
                for field_id in field_ids:
                    ax = fig.add_subplot(1, n_snapshots, index + 1, adjustable='box')
                    ax.title.set_text('time = {:.2f} hr'.format(time/60/60))
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



    def plot_agents(self, ax, agent_data, lattice_scaling, cell_radius, baseline_color):

        for agent_id, data in agent_data.iteritems():
            location = data['location']
            volume = data['volume']
            tags = data.get('tags')
            intensity = 0.0
            if tags:
                intensity = tags[tags.keys()[0]]  # only use first tag TODO -- multiple tags?

            y = location[0] * lattice_scaling
            x = location[1] * lattice_scaling
            theta = location[2]
            length = self.volume_to_length(volume, cell_radius) * lattice_scaling

            # plot cell as a 2D line
            width = linewidth_from_data_units(cell_radius * 2, ax) * lattice_scaling
            body_width = width * 0.95
            # TODO get body_length

            dx = length * np.sin(theta)
            dy = length * np.cos(theta)

            # adjust color to tag intensity
            agent_color = flourescent_color(baseline_color, intensity)
            rgb = hsv_to_rgb(agent_color)

            # plot outline
            ax.plot([y - dy / 2, y + dy / 2], [x - dx / 2, x + dx / 2],
                    linewidth=width,
                    color='k',
                    solid_capstyle='butt')

            # plot body
            ax.plot([y - dy / 2, y + dy / 2], [x - dx / 2, x + dx / 2],
                    linewidth=body_width,
                    color=rgb,
                    solid_capstyle='butt')


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

def flourescent_color(baseline_hsv, intensity):
    new_hsv = baseline_hsv[:]

    # saturation
    new_hsv[1] = max((1-intensity/100), 0.3)

    # value
    new_hsv[2] = min(intensity/20, 1)

    return new_hsv

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
    return linewidth * (length / value_range[0])