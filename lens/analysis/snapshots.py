from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb

from lens.analysis.analysis import Analysis, get_compartment
from lens.actor.process import deep_merge

# DEFAULT_COLOR = [color/255 for color in [102, 178, 255]]
DEFAULT_COLOR = [220/360, 100.0/100.0, 50.0/100.0]  # HSV
FLOURESCENT_COLOR = [120/360, 100.0/100.0, 100.0/100.0]  # HSV

# TODO (Eran) -- min/max should be an argument
MIN_PROTEIN = 75  # any value less than this shows no flourescence
MAX_PROTEIN = 110  # if a tagged protein has a value over this, it is fully saturated

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
                time_dict = deep_merge(dict(time_dict), tag_hist)

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

                fields = row.get('fields', [])
                time_dict[time].update({'fields': fields})

            return time_dict

    def analyze(self, experiment_config, data, output_dir):

        n_snapshots = 6  # number of snapshots

        phylogeny = experiment_config['phylogeny']
        time_data = data['environment']
        tags_data = data['compartments']

        time_vec = time_data.keys()
        edge_length_x = experiment_config['edge_length_x']
        edge_length_y = experiment_config['edge_length_y']
        cell_radius = experiment_config['cell_radius']

        # time steps that will be used
        plot_steps = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
        snapshot_times = [time_vec[i] for i in plot_steps]

        # number of fields
        field_ids = time_data[time_vec[0]].get('fields',{}).keys()
        n_fields = max(len(field_ids),1)

        # agent colors based on phylogeny
        agent_colors = {agent_id: [] for agent_id in phylogeny.keys()}
        ancestors = phylogeny.keys()
        descendents = list(set([daughter for daughters in phylogeny.values() for daughter in daughters]))
        initial_agents = np.setdiff1d(ancestors,descendents)
        for agent_id in initial_agents:
            agent_colors.update(color_phylogeny(agent_id, phylogeny, DEFAULT_COLOR))

        # make figure
        fig = plt.figure(figsize=(20*n_snapshots, 10*n_fields))
        grid = plt.GridSpec(n_fields, n_snapshots, wspace=0.2, hspace=0.01)
        plt.rcParams.update({'font.size': 36})
        for index, time in enumerate(snapshot_times, 0):
            field_data = time_data[time].get('fields')
            agent_data = time_data[time]['agents']

            if tags_data:
                agent_tags = {}
                for agent_id in agent_data.keys():
                    tdata = tags_data[agent_id]
                    tags = tdata.get(time) or tdata[min(tdata.keys(), key=lambda k: abs(k-time))]  # get closest time key
                    agent_tags[agent_id] = tags
                agent_data = deep_merge(dict(agent_data), agent_tags)

            if field_ids:
                # plot fields
                for f_index, field_id in enumerate(field_ids, 0):  # TODO --multiple fields

                    ax = fig.add_subplot(grid[f_index, index])  # grid is (row, column)
                    ax.title.set_text('time: {:.4f} hr | field: {}'.format(float(time)/60./60., field_id))
                    ax.set_xlim([0, edge_length_x])
                    ax.set_ylim([0, edge_length_y])
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])

                    # rotate field and plot
                    field = np.rot90(np.array(field_data[field_id])).tolist()
                    plt.imshow(field,
                               origin='lower',
                               extent=[0, edge_length_x, 0, edge_length_y],
                               interpolation='nearest',
                               cmap='YlGn')
                    self.plot_agents(ax, agent_data, cell_radius, agent_colors)

            else:
                ax = fig.add_subplot(1, n_snapshots, index + 1, adjustable='box')
                ax.title.set_text('time = {}'.format(time))
                ax.set_xlim([0, edge_length_x])
                ax.set_ylim([0, edge_length_y])
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                self.plot_agents(ax, agent_data, cell_radius, agent_colors)

        # plt.subplots_adjust(wspace=0.7, hspace=0.1)
        figname = '/snapshots'
        if tags_data:
            figname = '/snapshots_tagged'
        plt.savefig(output_dir + figname, bbox_inches='tight')
        plt.close(fig)


    def plot_agents(self, ax, agent_data, cell_radius, agent_colors):

        for agent_id, data in agent_data.iteritems():

            # location, orientation, length
            volume = data['volume']
            x = data['location'][0]
            y = data['location'][1]
            theta = data['location'][2] / np.pi * 180 + 90 # rotate 90 degrees to match field
            length = volume_to_length(volume, cell_radius)
            width = cell_radius * 2

            # colors and flourescent tags
            agent_color = agent_colors.get(agent_id, DEFAULT_COLOR)
            tags = data.get('tags')
            if tags:
                intensity = max((tags[tags.keys()[0]] - MIN_PROTEIN), 0)
                intensity = min(intensity / (MAX_PROTEIN - MIN_PROTEIN), 1)  # only use first tag TODO -- multiple tags?
                agent_color = flourescent_color(DEFAULT_COLOR, intensity)
            rgb = hsv_to_rgb(agent_color)

            # Create a rectangle
            rect = patches.Rectangle((x, y), width, length, theta, linewidth=1, edgecolor='k', facecolor=rgb)
            ax.add_patch(rect)


def color_phylogeny(ancestor_id, phylogeny, baseline_hsv, phylogeny_colors={}):
    # get colors for all descendants of the ancestor through recursive calls to each generation
    phylogeny_colors.update({ancestor_id: baseline_hsv})
    daughter_ids = phylogeny.get(ancestor_id)
    if daughter_ids:
        for daughter_id in daughter_ids:
            daughter_color = mutate_color(baseline_hsv)
            color_phylogeny(daughter_id, phylogeny, daughter_color)
    return phylogeny_colors

def mutate_color(baseline_hsv):
    new_hsv = baseline_hsv[:]
    new_h = new_hsv[0] + np.random.normal(0, 0.08)  # (mu, sigma)
    new_hsv[0] = new_h % 1  # hue
    return new_hsv

def flourescent_color(baseline_hsv, intensity):
    # move color towards bright green when intensity = 1
    new_hsv = baseline_hsv[:]
    distance = [a - b for a, b in zip(FLOURESCENT_COLOR, new_hsv)]
    new_hsv = [a + intensity*b for a, b in zip(new_hsv, distance)]
    return new_hsv

def volume_to_length(volume, cell_radius):
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
    return linewidth * (length / value_range[0])