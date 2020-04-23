from __future__ import absolute_import, division, print_function

import random
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb
from mpl_toolkits.axes_grid1 import make_axes_locatable

from vivarium.analysis.analysis import Analysis, get_compartment
from vivarium.utils.dict_utils import deep_merge


DEFAULT_COLOR = [220/360, 100.0/100.0, 70.0/100.0]  # HSV

# colors for phylogeny initial agents
PHYLOGENY_HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]

# colors for flourescent tags
BASELINE_TAG_COLOR = [220/360, 100.0/100.0, 30.0/100.0]  # HSV
FLOURESCENT_COLOR = [120/360, 100.0/100.0, 100.0/100.0]  # HSV

N_SNAPSHOTS = 6  # number of snapshots

class Snapshots(Analysis):
    def __init__(self):
        super(Snapshots, self).__init__(analysis_type='env_with_compartment')

    def get_data(self, client, query, options={}):

        type = options.get('type')
        tags = options.get('tags')
        if type is 'compartment' and tags:
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

            # arrange data by time, for easy integration with environment data
            time_dict = {time: {} for time in times}
            for tag, series in tags_history.items():
                tag_hist = {time: {'tags': {tag: state}} for time, state in zip(times,series)}
                time_dict = deep_merge(dict(time_dict), tag_hist)

            return time_dict

        elif type is 'environment':
            # get data about agent locations ('type': 'lattice-agent')
            query_lattice = {'type': 'lattice-agent'}
            query_lattice.update(query)
            data_lattice = client.find(query_lattice)
            data_lattice.sort('time')

            # get data about concentrations ('type': 'lattice-fields')
            query_field = {'type': 'lattice-fields'}
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

        phylogeny = experiment_config['phylogeny']
        time_data = data['environment']
        compartments = data['compartments']
        edge_length_x = experiment_config['edge_length_x']
        edge_length_y = experiment_config['edge_length_y']
        cell_radius = experiment_config['cell_radius'] # TODO --get cell_radius from compartments

        time_vec = list(time_data.keys())

        # time steps that will be used
        plot_steps = np.round(np.linspace(0, len(time_vec) - 1, N_SNAPSHOTS)).astype(int)
        snapshot_times = [time_vec[i] for i in plot_steps]

        # number of fields
        field_ids = list(time_data[time_vec[0]].get('fields',{}).keys())
        n_fields = max(len(field_ids),1)

        ## get tag ids and range
        # TODO -- make into a function in vivarium.analysis.analysis
        tag_range = {}
        for c_id, c_data in compartments.items():
            # if this compartment has tags, get their ids and range
            if c_data:
                # get volume from time_data
                vol_data = {}
                for t, t_data in time_data.items():
                    # if c_id present, save its volume?
                    if c_id in t_data['agents']:
                        vol_data[t] = t_data['agents'][c_id]['volume']

                # go through each time point to get tag counts
                for t, t_data in c_data.items():
                    volume = vol_data.get(t) or \
                        vol_data[min(vol_data.keys(), key=lambda k: abs(k - t))]  # get closest time key
                    c_tags = t_data['tags']
                    for tag_id, count in c_tags.items():
                        conc = count / volume
                        if tag_id in tag_range:
                            tag_range[tag_id] = [
                                min(tag_range[tag_id][0], conc),
                                max(tag_range[tag_id][1], conc)]
                        else:
                            tag_range[tag_id] = [conc, conc]

        # get fields' range
        field_range = {}
        for field_id in field_ids:
            for t, t_data in time_data.items():
                field = t_data['fields'][field_id]
                concs = [conc for row in field for conc in row]
                min_conc = min(concs)
                max_conc = max(concs)
                if field_id in field_range:
                    field_range[field_id] = [
                        min(field_range[field_id][0], min_conc),
                        max(field_range[field_id][1], max_conc)]
                else:
                    field_range[field_id] = [min_conc, max_conc]

        # initial agent ids
        if phylogeny:
            # find initial agents in phylogeny
            ancestors = list(phylogeny.keys())
            descendents = list(set([daughter
                for daughters in phylogeny.values() for daughter in daughters]))
            initial_agents = np.setdiff1d(ancestors,descendents)
        else:
            # if no phylogeny, all compartments must be initial agents
            initial_agents = np.array(list(compartments.keys()))

        # agent colors based on phylogeny
        agent_colors = {agent_id: [] for agent_id in phylogeny.keys()}
        for agent_id in initial_agents:
            hue = random.choice(PHYLOGENY_HUES)  # select random initial hue
            initial_color = [hue] + DEFAULT_SV
            agent_colors.update(color_phylogeny(agent_id, phylogeny, initial_color))

        ## make the figure
        # fields and tag data are plotted in separate rows
        n_rows = len(tag_range) + n_fields
        n_cols = N_SNAPSHOTS + 2  # one column for text, one for the colorbar
        fig = plt.figure(figsize=(12*n_cols, 12*n_rows))
        grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
        plt.rcParams.update({'font.size': 36})

        # text in first column
        row_idx = 0
        for field_id in field_ids:
            fig.add_subplot(grid[row_idx, 0])
            plt.text(0.02, 0.9, '{}'.format(field_id), fontsize=36)
            plt.axis('off')
            row_idx+=1
        for tag_id in list(tag_range.keys()):
            fig.add_subplot(grid[row_idx, 0])
            plt.text(0.02, 0.9, '{}'.format(tag_id), fontsize=36)
            plt.axis('off')
            row_idx+=1

        # plot snapshot data in each subsequent column
        for col_idx, time in enumerate(snapshot_times, 1):
            row_idx = 0
            field_data = time_data[time].get('fields')
            agent_data = time_data[time]['agents']

            if not field_ids and not tag_range.keys():
                ax = init_axes(
                    fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time)
                plot_agents(ax, agent_data, cell_radius, agent_colors)

            else:
                # plot fields
                for field_id in field_ids:
                    ax = init_axes(
                        fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time)

                    # transpose field to align with agent
                    field = np.transpose(np.array(field_data[field_id])).tolist()
                    if field_range[field_id]:
                        vmin, vmax = field_range[field_id]
                        im = plt.imshow(field,
                                        origin='lower',
                                        extent=[0, edge_length_x, 0, edge_length_y],
                                        vmin=vmin,
                                        vmax=vmax,
                                        cmap='BuPu')
                    else:
                        im = plt.imshow(field,
                                        origin='lower',
                                        extent=[0, edge_length_x, 0, edge_length_y],
                                        interpolation='nearest',
                                        cmap='BuPu')
                    plot_agents(ax, agent_data, cell_radius, agent_colors)

                    # colorbar in new column after final snapshot
                    if col_idx == N_SNAPSHOTS:
                        cbar_col = col_idx+1
                        ax = fig.add_subplot(grid[row_idx, cbar_col])
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes("left", size="5%", pad=0.0)
                        fig.colorbar(im, cax=cax, format='%.6f')
                        ax.axis('off')

                    row_idx += 1

            # plot tags
            for tag_id in list(tag_range.keys()):
                ax = init_axes(
                    fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time)
                ax.set_facecolor('palegoldenrod')  # set background color

                # update agent colors based on tag_level
                agent_tag_colors = {}
                for agent_id in agent_data.keys():
                    agent_color = BASELINE_TAG_COLOR
                    if agent_id in compartments:
                        tdata = compartments[agent_id]
                        all_tags = tdata.get(time) \
                            or tdata[min(tdata.keys(), key=lambda k: abs(k-time))]  # get closest time key

                        # get current tag concentration, and determine color
                        counts = all_tags['tags'][tag_id]
                        volume = agent_data[agent_id]['volume']
                        level = counts / volume
                        min_tag, max_tag = tag_range[tag_id]
                        if min_tag != max_tag:
                            intensity = max((level - min_tag), 0)
                            intensity = min(intensity / (max_tag - min_tag), 1)
                            agent_color = flourescent_color(BASELINE_TAG_COLOR, intensity)

                    agent_tag_colors[agent_id] = agent_color

                plot_agents(ax, agent_data, cell_radius, agent_tag_colors)

                row_idx += 1

        if tag_range:
             figname = '/snap_out_tagged'
        else:
            figname = '/snap_out'

        plt.subplots_adjust(wspace=0.1, hspace=1.0)
        plt.savefig(output_dir + figname)
        plt.close(fig)


def plot_agents(ax, agent_data, cell_radius, agent_colors):

    for agent_id, data in agent_data.items():

        # location, orientation, length
        volume = data['volume']
        x_center = data['location'][0]
        y_center = data['location'][1]
        theta = data['location'][2] / np.pi * 180 + 90 # rotate 90 degrees to match field
        length = volume_to_length(volume, cell_radius)
        width = cell_radius * 2

        # get bottom left position
        x_offset = (width / 2)
        y_offset = (length / 2)
        theta_rad = math.radians(theta)
        dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
        dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

        x = x_center - dx
        y = y_center - dy

        # get color, convert to rgb
        agent_color = agent_colors.get(agent_id, DEFAULT_COLOR)
        rgb = hsv_to_rgb(agent_color)

        # Create a rectangle
        rect = patches.Rectangle((x, y), width, length, angle=theta, linewidth=2, edgecolor='w', facecolor=rgb)

        ax.add_patch(rect)

def init_axes(fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.4f} hr'.format(float(time) / 60. / 60.)
        plt.title(plot_title, y=1.08)
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    return ax

def color_phylogeny(ancestor_id, phylogeny, baseline_hsv, phylogeny_colors={}):
    """
    get colors for all descendants of the ancestor
    through recursive calls to each generation
    """
    phylogeny_colors.update({ancestor_id: baseline_hsv})
    daughter_ids = phylogeny.get(ancestor_id)
    if daughter_ids:
        for daughter_id in daughter_ids:
            daughter_color = mutate_color(baseline_hsv)
            color_phylogeny(daughter_id, phylogeny, daughter_color)
    return phylogeny_colors

def mutate_color(baseline_hsv):
    new_hsv = baseline_hsv[:]
    new_h = new_hsv[0] + np.random.uniform(-0.2, 0.2)
    new_hsv[0] = new_h % 1
    return new_hsv

def flourescent_color(baseline_hsv, intensity):
    # move color towards bright green (FLOURESCENT_COLOR) when intensity = 1
    new_hsv = baseline_hsv[:]
    distance = [a - b for a, b in zip(FLOURESCENT_COLOR, new_hsv)]
    new_hsv = [a + intensity*b for a, b in zip(new_hsv, distance)]
    return new_hsv

def volume_to_length(volume, cell_radius):
    """
    get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
    a=length of cylinder without rounded caps, l=total length:
    V = (4/3)*PI*r^3 + PI*r^2*a
    l = a + 2*r
    """
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
