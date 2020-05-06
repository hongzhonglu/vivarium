import copy
import csv
import os
import io

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from vivarium.compartment import emitter as emit
from vivarium.utils.dict_utils import (
    deep_merge,
    deep_merge_check,
    flatten_timeseries,
)
from vivarium.compartment.process import (
    initialize_state,
    Compartment,
    COMPARTMENT_STATE,
    Process,
    Store)
from vivarium.utils.units import units

# processes
from vivarium.processes.derive_globals import DeriveGlobals, AVOGADRO
from vivarium.processes.derive_counts import DeriveCounts
from vivarium.processes.derive_concentrations import DeriveConcs
from vivarium.processes.derive_mass import DeriveMass


REFERENCE_DATA_DIR = os.path.join('vivarium', 'reference_data')
TEST_OUT_DIR = os.path.join('out', 'tests')


deriver_library = {
    'mass': DeriveMass,
    'globals': DeriveGlobals,
    'mmol_to_counts': DeriveCounts,
    'counts_to_mmol': DeriveConcs,
}



def get_derivers(process_list, topology, deriver_config={}):
    '''
    get the derivers for a list of processes

    requires:
        - process_list: (list) with configured processes
        - topology: (dict) with topology of the processes connected to compartment ports
        - config: (dict) with deriver configurations, which are used to make deriver processes

    returns: (dict) with:
        {'deriver_processes': processes,
        'deriver_topology': topology}
    '''

    # get deriver configuration from processes
    process_derivers = get_deriver_config_from_proceses(process_list, topology)
    deriver_configs = process_derivers['deriver_configs']
    deriver_topology = process_derivers['deriver_topology']

    # update deriver_configs
    deriver_configs = deep_merge(deriver_configs, deriver_config)

    # update topology based on deriver_config
    for process_id, config in deriver_config.items():
        if process_id not in deriver_topology:
            try:
                ports = config['ports']
                deriver_topology[process_id] = ports
            except:
                print('{} deriver requires topology in deriver_config'.format(process_id))
                raise

    # configure the deriver processes
    deriver_processes = {}
    for deriver_type, deriver_config in deriver_configs.items():
        deriver_processes[deriver_type] = deriver_library[deriver_type](deriver_config)

    # TODO -- put derivers in order
    processes = [deriver_processes]

    return {
        'deriver_processes': processes,
        'deriver_topology': deriver_topology}

def get_deriver_config_from_proceses(process_list, topology):
    ''' get the deriver configuration from processes' deriver_settings'''

    deriver_configs = {}
    full_deriver_topology = {}

    for level in process_list:
        for process_id, process in level.items():
            process_settings = process.default_settings()
            deriver_setting = process_settings.get('deriver_setting', [])
            try:
                port_map = topology[process_id]
            except:
                print('{} topology port mismatch'.format(process_id))
                raise

            for setting in deriver_setting:
                deriver_type = setting['type']
                keys = setting['keys']
                source_port = setting['source_port']
                target_port = setting['derived_port']
                try:
                    source_compartment_port = port_map[source_port]
                    target_compartment_port = port_map[target_port]
                except:
                    print('source/target port mismatch for process "{}"'.format(process_id))
                    raise

                # make deriver_topology, add to full_deriver_topology
                deriver_topology = {
                    deriver_type: {
                        source_port: source_compartment_port,
                        target_port: target_compartment_port,
                        'global': 'global'}}
                deep_merge(full_deriver_topology, deriver_topology)

                # TODO -- what if multiple different source/targets?
                # TODO -- merge overwrites them. need list extend
                ports_config = {
                    'source_ports': {source_port: keys},
                    'target_ports': {target_port: keys}}

                # ports for configuration
                deriver_config = {deriver_type: ports_config}
                deep_merge(deriver_configs, deriver_config)

    return {
        'deriver_configs': deriver_configs,
        'deriver_topology': full_deriver_topology}


def get_schema(process_list, topology):
    schema = {}
    for level in process_list:
        for process_id, process in level.items():
            process_settings = process.default_settings()
            process_schema = process_settings.get('schema', {})
            try:
                port_map = topology[process_id]
            except:
                print('{} topology port mismatch'.format(process_id))
                raise

            # go through each port, and get the schema
            for process_port, settings in process_schema.items():
                compartment_port = port_map[process_port]
                compartment_schema = {
                    compartment_port: settings}

                ## TODO -- check for mismatch
                deep_merge_check(schema, compartment_schema)

    return schema

def process_in_compartment(process, settings={}):
    ''' put a process in a compartment, with all derivers added '''
    process_settings = process.default_settings()
    compartment_state_port = settings.get('compartment_state_port')
    emitter = settings.get('emitter', 'timeseries')
    deriver_config = settings.get('deriver_config', {})

    processes = [{'process': process}]
    topology = {
        'process': {
            port: port for port in process.ports
            if (not compartment_state_port
                or port != compartment_state_port)
        }
    }

    if compartment_state_port:
        topology['process'][compartment_state_port] = COMPARTMENT_STATE

    # add derivers
    derivers = get_derivers(processes, topology, deriver_config)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])

    # make the state
    state_dict = process_settings['state']
    states = initialize_state(
        all_processes,
        topology,
        state_dict)

    options = {
        'topology': topology,
        'emitter': emitter}

    return Compartment(processes, deriver_processes, states, options)

def simulate_process_with_environment(process, settings={}):
    ''' simulate a process in a compartment with an environment '''
    compartment = process_in_compartment(process, settings)
    return simulate_with_environment(compartment, settings)

def simulate_process(process, settings={}):
    ''' simulate a process in a compartment with no environment '''
    compartment = process_in_compartment(process, settings)
    return simulate_compartment(compartment, settings)

def simulate_with_environment(compartment, settings={}):
    '''
    run a compartment simulation with an environment.
    Requires:
        - a compartment with environment_port and exchange_port

    Returns:
        - a timeseries of variables from all ports.
        - if 'return_raw_data' is True, it returns the raw data instead
    '''

    # parameters
    nAvogadro = AVOGADRO

    # get environment configuration
    environment_port = settings['environment_port']
    env_volume = settings.get('environment_volume', 1e-12) * units.L
    exchange_port = settings.get('exchange_port')
    if exchange_port:
        exchange_ids = list(compartment.states[exchange_port].keys())
    else:
        print('no exchange port! simulate environment without exchange')
    environment = compartment.states.get(environment_port)
    exchange = compartment.states.get(exchange_port)

    # get timeline
    total_time = settings.get('total_time', 10)
    timeline = copy.deepcopy(settings.get('timeline', [(total_time, {})]))
    end_time = timeline[-1][0]
    timestep = compartment.time_step

    # data settings
    return_raw_data = settings.get('return_raw_data', False)

    ## run simulation
    time = 0
    while time < end_time:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for port_id, change in change_dict.items():
                    port = compartment.states.get(port_id)
                    port.assign_values(change)
                timeline.pop(0)

        # update compartment
        compartment.update(timestep)

        ## apply exchange to environment
        # get counts, convert to change in concentration
        if exchange:
            delta_counts = exchange.state_for(exchange_ids)
            mmol_to_counts = (nAvogadro.to('1/mmol') * env_volume).to('L/mmol').magnitude
            delta_concs = {mol_id: counts / mmol_to_counts for mol_id, counts in delta_counts.items()}
            environment.apply_update(delta_concs)

            # reset exchange
            reset_exchange = {key: 0 for key in exchange_ids}
            exchange.assign_values(reset_exchange)

    if return_raw_data:
        return compartment.emitter.get_data()
    else:
        return compartment.emitter.get_timeseries()



# plotting functions
def plot_compartment_topology(compartment, settings, out_dir='out', filename='topology'):
    """
    Make a plot of the topology
     - compartment: a compartment
     - settings (dict): 'network_layout' can be 'bipartite' or 'process_layers'
    """
    store_rgb = [x/255 for x in [239,131,148]]
    process_rgb = [x / 255 for x in [249, 204, 86]]
    node_size = 2500
    node_distance = 1
    layer_distance = 10

    topology = compartment.topology
    process_layers = compartment.processes

    # get figure settings
    show_ports = settings.get('show_ports', True)

    # make graph from topology
    G = nx.Graph()
    process_nodes = []
    store_nodes = []
    edges = {}
    for process_id, connections in topology.items():
        process_nodes.append(process_id)
        G.add_node(process_id)

        for port, store_id in connections.items():
            if store_id not in store_nodes:
                store_nodes.append(store_id)
            if store_id not in list(G.nodes):
                G.add_node(store_id)

            edge = (process_id, store_id)
            edges[edge]= port

            G.add_edge(process_id, store_id)

    # are there overlapping names?
    overlap = [name for name in process_nodes if name in store_nodes]
    if overlap:
        print('{} shared by processes and stores'.format(overlap))


    # get positions
    pos = {}
    n_rows = max(len(process_nodes), len(store_nodes))
    plt.figure(3, figsize=(12, 1.2 * n_rows))

    for idx, node_id in enumerate(process_nodes, 1):
        pos[node_id] = np.array([-1, -idx*node_distance])
    for idx, node_id in enumerate(store_nodes, 1):
        pos[node_id] = np.array([1, -idx*node_distance])


    # plot
    nx.draw_networkx_nodes(G, pos,
                           nodelist=process_nodes,
                           with_labels=True,
                           node_color=process_rgb,
                           node_size=node_size,
                           node_shape='o')
    nx.draw_networkx_nodes(G, pos,
                           nodelist=store_nodes,
                           with_labels=True,
                           node_color=store_rgb,
                           node_size=node_size,
                           node_shape='s')

    # edges
    colors = list(range(1,len(edges)+1))
    nx.draw_networkx_edges(G, pos,
                           edge_color=colors,
                           width=1.5)

    # labels
    nx.draw_networkx_labels(G, pos,
                            font_size=8,
                            )
    if show_ports:
        nx.draw_networkx_edge_labels(G, pos,
                                 edge_labels=edges,
                                 font_size=6,
                                 label_pos=0.85)

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.figure(3, figsize=(12, 12))
    plt.axis('off')
    plt.savefig(fig_path, bbox_inches='tight')

    plt.close()


def set_axes(ax, show_xaxis=False):
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)


def plot_simulation_output(timeseries, settings={}, out_dir='out', filename='simulation'):
    '''
    plot simulation output, with rows organized into separate columns.

    Requires:
        - timeseries (dict). This can be obtained from simulation output with convert_to_timeseries()
        - settings (dict) with:
            {
            'max_rows': (int) ports with more states than this number of states get wrapped into a new column
            'remove_zeros': (bool) if True, timeseries with all zeros get removed
            'remove_flat': (bool) if True, timeseries with all the same value get removed
            'skip_ports': (list) entire ports that won't be plotted
            'overlay': (dict) with
                {'bottom_port': 'top_port'}  ports plotted together by matching state_ids, with 'top_port' in red
            'show_state': (list) with [('port_id', 'state_id')]
                for all states that will be highlighted, even if they are otherwise to be removed
            }
    TODO -- some molecules have 'inf' concentrations for practical reasons. How should these be plotted?
    '''

    skip_keys = ['time']

    # get settings
    max_rows = settings.get('max_rows', 25)
    remove_zeros = settings.get('remove_zeros', False)
    remove_flat = settings.get('remove_flat', False)
    skip_ports = settings.get('skip_ports', [])
    overlay = settings.get('overlay', {})
    show_state = settings.get('show_state', [])
    top_ports = list(overlay.values())
    bottom_ports = list(overlay.keys())

    ports = [port for port in timeseries.keys() if port not in skip_keys + skip_ports]
    time_vec = timeseries['time']

    # remove selected states
    removed_states = []
    if remove_flat:
        # find series with all the same value
        for port in ports:
            for state_id, series in timeseries[port].items():
                if series.count(series[0]) == len(series):
                    removed_states.append((port, state_id))
    elif remove_zeros:
        # find series with all zeros
        for port in ports:
            for state_id, series in timeseries[port].items():
                if all(v == 0 for v in series):
                    removed_states.append((port, state_id))

    # if specified in show_state, keep in timeseries
    for port_state in show_state:
        if port_state in removed_states:
            removed_states.remove(port_state)

    # remove from timeseries
    for (port, state_id) in removed_states:
        del timeseries[port][state_id]

    # get the number of states in each port
    n_data = [len(timeseries[key]) for key in ports if key not in top_ports]
    if 0 in n_data:
        n_data.remove(0)

    # limit number of rows to max_rows by adding new columns
    columns = []
    for n_states in n_data:
        new_cols = n_states / max_rows
        if new_cols > 1:
            for col in range(int(new_cols)):
                columns.append(max_rows)

            mod_states = n_states % max_rows
            if mod_states > 0:
                columns.append(mod_states)
        else:
            columns.append(n_states)
    n_cols = len(columns)
    n_rows = max(columns)

    # make figure and plot
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 1.5))
    grid = plt.GridSpec(n_rows, n_cols)

    row_idx = 0
    col_idx = 0
    for port in ports:
        top_timeseries = {}

        # set up overlay
        if port in bottom_ports:
            top_port = overlay[port]
            top_timeseries = timeseries[top_port]
        elif port in top_ports + skip_ports:
            # don't give this row its own plot
            continue

        for state_id, series in sorted(timeseries[port].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            # check if series is a list of ints or floats
            if not all(isinstance(state, (int, float, np.int64, np.int32)) for state in series):
                ax.title.set_text(str(port) + ': ' + str(state_id) + ' (non numeric)')
                ax.title.set_fontsize(16)
            else:
                # plot line at zero if series crosses the zero line
                if any(x == 0.0 for x in series) or (any(x < 0.0 for x in series) and any(x > 0.0 for x in series)):
                    zero_line = [0 for t in time_vec]
                    ax.plot(time_vec, zero_line, 'k--')

                if (port, state_id) in show_state:
                    ax.plot(time_vec, series, 'indigo', linewidth=2)
                else:
                    ax.plot(time_vec, series)

                # overlay
                if state_id in top_timeseries.keys():
                    ax.plot(time_vec, top_timeseries[state_id], 'm', label=top_port)
                    ax.legend()

                ax.title.set_text(str(port) + ': ' + str(state_id))
                ax.title.set_fontsize(16)

            if row_idx == columns[col_idx]-1:
                # if last row of column
                set_axes(ax, True)
                ax.set_xlabel('time')
                row_idx = 0
                col_idx += 1
            else:
                set_axes(ax)
                row_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')


# timeseries functions
def save_timeseries(timeseries, out_dir='out'):
    '''Save a timeseries as a CSV in out_dir'''
    flattened = flatten_timeseries(timeseries)
    rows = np.transpose(list(flattened.values())).tolist()
    with open(os.path.join(out_dir, 'simulation_data.csv'), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(flattened.keys())
        writer.writerows(rows)


def load_timeseries(path_to_csv):
    '''Load a timeseries saved as a CSV using save_timeseries.

    The timeseries is returned in flattened form.
    '''
    with io.open(path_to_csv, 'r', newline='') as f:
        reader = csv.DictReader(f)
        timeseries = {}
        for row in reader:
            for header, elem in row.items():
                timeseries.setdefault(header, []).append(float(elem))
    return timeseries


def timeseries_to_ndarray(timeseries, keys=None):
    if keys is None:
        keys = timeseries.keys()
    filtered = {key: timeseries[key] for key in keys}
    array = np.array(list(filtered.values()))
    return array


def _prepare_timeseries_for_comparison(
    timeseries1, timeseries2, keys=None,
    required_frac_checked=0.9,
):
    '''Prepare two timeseries for comparison

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Returns:
        A tuple of an ndarray for each of the two timeseries and a list of
        the keys for the rows of the arrays. Each ndarray has a row for
        each key, in the order of keys. The ndarrays have only the
        columns corresponding to the timepoints common to both
        timeseries.

    Raises:
        AssertionError: If a correlation is strictly below the
            threshold or if too few timepoints are common to both
            timeseries.
    '''
    if 'time' not in timeseries1 or 'time' not in timeseries2:
        raise AssertionError('Both timeseries must have key "time"')
    if keys is None:
        keys = set(timeseries1.keys()) & set(timeseries2.keys())
    else:
        if 'time' not in keys:
            keys.append('time')
    keys = list(keys)
    time_index = keys.index('time')
    shared_times = set(timeseries1['time']) & set(timeseries2['time'])
    frac_timepoints_checked = (
        len(shared_times)
        / min(len(timeseries1), len(timeseries2))
    )
    if frac_timepoints_checked < required_frac_checked:
        raise AssertionError(
            'The timeseries share too few timepoints: '
            '{} < {}'.format(
                frac_timepoints_checked, required_frac_checked)
        )
    array1 = timeseries_to_ndarray(timeseries1, keys)
    array2 = timeseries_to_ndarray(timeseries2, keys)
    shared_times_mask1 = np.isin(array1[time_index], list(shared_times))
    shared_times_mask2 = np.isin(array2[time_index], list(shared_times))
    return (
        array1[:, shared_times_mask1],
        array2[:, shared_times_mask2],
        keys,
    )


def assert_timeseries_correlated(
    timeseries1, timeseries2, keys=None,
    default_threshold=(1 - 1e-10), thresholds={},
    required_frac_checked=0.9,
):
    '''Check that two timeseries are correlated.

    Uses a Pearson correlation coefficient. Only the data from
    timepoints common to both timeseries are compared.

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        default_threshold: The threshold correlation coefficient to use
            when a threshold is not specified in thresholds.
        thresholds: Dictionary of key-value pairs where the key is a key
            in both timeseries and the value is the threshold
            correlation coefficient to use when checking that key
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Raises:
        AssertionError: If a correlation is strictly below the
            threshold or if too few timepoints are common to both
            timeseries.
    '''
    array1, array2, keys = _prepare_timeseries_for_comparison(
        timeseries1, timeseries2, keys, required_frac_checked)
    for index, key in enumerate(keys):
        corrcoef = np.corrcoef(
            array1[index],
            array2[index],
        )[0][1]
        threshold = thresholds.get(key, default_threshold)
        if corrcoef < threshold:
            raise AssertionError(
                'The correlation coefficient for '
                '{} is too small: {} < {}'.format(
                    key, corrcoef, threshold)
            )


def assert_timeseries_close(
    timeseries1, timeseries2, keys=None,
    default_tolerance=(1 - 1e-10), tolerances={},
    required_frac_checked=0.9,
):
    '''Check that two timeseries are similar.

    Ensures that each pair of data points between the two timeseries are
    within a tolerance of each other, after filtering out timepoints not
    common to both timeseries.

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        default_tolerance: The tolerance to use when not specified in
            tolerances.
        tolerances: Dictionary of key-value pairs where the key is a key
            in both timeseries and the value is the tolerance to use
            when checking that key.
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Raises:
        AssertionError: If a pair of data points have a difference
            strictly above the tolerance threshold or if too few
            timepoints are common to both timeseries.
    '''
    array1, array2, keys = _prepare_timeseries_for_comparison(
        timeseries1, timeseries2, keys, required_frac_checked)
    for index, key in enumerate(keys):
        tolerance = tolerances.get(key, default_tolerance)
        if not np.allclose(
            array1[index], array2[index], atol=tolerance
        ):
            raise AssertionError(
                'The data for {} differed by more than {}'.format(
                    key, tolerance)
            )


# TESTS


class ToyLinearGrowthDeathProcess(Process):

    GROWTH_RATE = 1.0
    THRESHOLD = 5.0

    def __init__(self, initial_parameters={}):
        ports = {
            'compartment': ['processes'],
            'global': ['mass'],
        }
        super(ToyLinearGrowthDeathProcess, self).__init__(
            ports, initial_parameters)

    def default_settings(self):
        default_settings = {
            'emitter_keys': {
                'global': ['mass']},
            'state': {
                'global': {
                    'mass': 0.0
                }
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        mass = states['global']['mass']
        mass_grown = (
            ToyLinearGrowthDeathProcess.GROWTH_RATE * timestep)
        update = {
            'global': {'mass': mass_grown},
        }
        if mass > ToyLinearGrowthDeathProcess.THRESHOLD:
            update['compartment'] = {
                'processes': [],
            }
        return update


class TestSimulateProcess:

    def test_compartment_state_port(self):
        '''Check that compartment state ports are handled'''
        process = ToyLinearGrowthDeathProcess()
        settings = {
            'compartment_state_port': 'compartment',
        }
        timeseries = simulate_process(process, settings)
        expected_masses = [
            # Mass stops increasing the iteration after mass > 5 because
            # cell dies
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.0, 7.0]
        masses = timeseries['global']['mass']
        assert masses == expected_masses

def load_compartment(composite, boot_config={}):
    '''
    put a composite function into a compartment

    inputs:
        - composite is a function that returns a dict with 'processes', 'states', and 'options'
        for configuring a compartment
        - boot_config (dict) with specific parameters for the processes
    return:
        - a compartment object for testing
    '''

    composite_config = composite(boot_config)
    processes = composite_config['processes']
    derivers = composite_config.get('derivers', [])
    states = composite_config['states']
    options = composite_config['options']
    options['emitter'] = boot_config.get('emitter', 'timeseries')

    return Compartment(processes, derivers, states, options)


def simulate_compartment(compartment, settings={}):
    '''
    run a compartment simulation
        Requires:
        - a compartment

    Returns:
        - a timeseries of variables from all ports.
        - if 'return_raw_data' is True, it returns the raw data instead
    '''

    timestep = settings.get('timestep', 1)
    total_time = settings.get('total_time', 10)

    # data settings
    return_raw_data = settings.get('return_raw_data', False)

    # run simulation
    time = 0
    while time < total_time:
        time += timestep
        compartment.update(timestep)

    if return_raw_data:
        return compartment.emitter.get_data()
    else:
        return compartment.emitter.get_timeseries()


## functions for testing
def toy_composite(config):
    '''
    a toy composite function for testing
    returns a dictionary with 'processes', 'states', and 'options'

    '''

    # toy processes
    class ToyMetabolism(Process):
        def __init__(self, initial_parameters={}):
            ports = {'pool': ['GLC', 'MASS']}
            parameters = {'mass_conversion_rate': 1}
            parameters.update(initial_parameters)

            super(ToyMetabolism, self).__init__(ports, parameters)

        def default_settings(self):
            return {
                'emitter_keys': {
                    port_id: keys for port_id, keys in self.ports.items()}
            }

        def next_update(self, timestep, states):
            update = {}
            glucose_required = timestep / self.parameters['mass_conversion_rate']
            if states['pool']['GLC'] >= glucose_required:
                update = {
                    'pool': {
                        'GLC': -2,
                        'MASS': 1}}

            return update

    class ToyTransport(Process):
        def __init__(self, initial_parameters={}):
            ports = {
                'external': ['GLC'],
                'internal': ['GLC']}
            parameters = {'intake_rate': 2}
            parameters.update(initial_parameters)

            super(ToyTransport, self).__init__(ports, parameters)

        def default_settings(self):
            return {
                'emitter_keys': {
                    port_id: keys for port_id, keys in self.ports.items()}
            }

        def next_update(self, timestep, states):
            update = {}
            intake = timestep * self.parameters['intake_rate']
            if states['external']['GLC'] >= intake:
                update = {
                    'external': {'GLC': -2, 'MASS': 1},
                    'internal': {'GLC': 2}}

            return update

    class ToyDeriveVolume(Process):
        def __init__(self, initial_parameters={}):
            ports = {
                'compartment': ['MASS', 'DENSITY', 'VOLUME']}
            parameters = {}

            super(ToyDeriveVolume, self).__init__(ports, parameters)

        def default_settings(self):
            return {
                'emitter_keys': {
                    port_id: keys for port_id, keys in self.ports.items()}
            }

        def next_update(self, timestep, states):
            volume = states['compartment']['MASS'] / states['compartment']['DENSITY']
            update = {
                'compartment': {'VOLUME': volume}}

            return update

    class ToyDeath(Process):
        def __init__(self, initial_parameters={}):
            ports = {
                'compartment': ['VOLUME'],
                'global': ['processes']}
            super(ToyDeath, self).__init__(ports, {})

        def next_update(self, timestep, states):
            volume = states['compartment']['VOLUME']
            update = {}

            if volume > 1.0:
                # kill the cell
                update = {
                    'global': {
                        'processes': []}}

            return update


    processes = [
        {'metabolism': ToyMetabolism(
            initial_parameters={
                'mass_conversion_rate': 0.5}), # example of overriding default parameters
         'transport': ToyTransport()},
        {'death': ToyDeath()}]

    # deriver processes
    derivers_processes = [
        {'external_volume': ToyDeriveVolume(),
         'internal_volume': ToyDeriveVolume()}]

    # declare the states
    states = {
        'periplasm': Store(
            initial_state={'GLC': 20, 'MASS': 100, 'DENSITY': 10, 'VOLUME': 100/10},
            schema={
                'VOLUME': {
                    'updater': 'set'}}),
        'cytoplasm': Store(
            initial_state={'MASS': 3, 'DENSITY': 10, 'VOLUME': 3/10},
            schema={
                'VOLUME': {
                    'updater': 'set'}})}

    # hook up the ports in each process to compartment states
    topology = {
        'metabolism': {
            'pool': 'cytoplasm'},
        'transport': {
            'external': 'periplasm',
            'internal': 'cytoplasm'},
        'death': {
            'compartment': 'cytoplasm',
            'global': COMPARTMENT_STATE},
        'external_volume': {
            'compartment': 'periplasm'},
        'internal_volume': {
            'compartment': 'cytoplasm'}}

    # emitter that prints to the terminal
    emitter = emit.get_emitter({
        'type': 'print',
        'keys': {
            'periplasm': ['GLC', 'MASS'],
            'cytoplasm': ['MASS']}})

    # schema for states
    schema = {}

    options = {
        # 'environment_port': 'environment',
        # 'exchange_port': 'exchange',
        'schema': schema,
        'emitter': emitter,
        'topology': topology,
        'initial_time': 0.0}

    return {
        'processes': processes,
        'derivers': derivers_processes,
        'states': states,
        'options': options}


def test_compartment(composite=toy_composite):
    compartment = load_compartment(composite)
    settings = {
        'timestep': 1,
        'total_time': 20,
        'emit_timeseries': True,}

    return simulate_compartment(compartment, settings)


if __name__ == '__main__':
    timeseries = test_compartment()
    print(timeseries)
