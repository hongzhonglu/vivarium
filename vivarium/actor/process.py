from __future__ import absolute_import, division, print_function

import copy
import collections
import numpy as np
import vivarium.actor.emitter as emit
import random
import os
from scipy import constants
import matplotlib.pyplot as plt

from vivarium.utils.dict_utils import merge_dicts


def npize(d):
    ''' Turn a dict into an ordered set of keys and values. '''

    ordered = [[key, value] for key, value in d.items()]
    keys = [key for key, _ in ordered]
    values = np.array([value for _, value in ordered], np.float64)

    return keys, values

def update_delta(key, state_dict, current_value, new_value):
    return current_value + new_value, {}

def update_set(key, state_dict, current_value, new_value):
    return new_value, {}

# def accumulate_delta(key, state_dict, current_value, new_value):
#     new_key = key + exchange_key
#     return current_value, {new_key: state_dict[new_key] + new_value}

updater_library = {
    'delta': update_delta,
    'set': update_set,
    'accumulate': update_delta,  # TODO -- remove accumulate
}


KEY_TYPE = 'U31'


class State(object):
    ''' Represents a set of named values. '''

    def __init__(self, initial_state={}, updaters={}):
        ''' Keys and state initialize empty, with a maximum key length of 31. '''

        self.state = copy.deepcopy(initial_state)
        self.updaters = updaters

    def keys(self):
        return self.state.keys()

    def duplicate(self, initial_state={}):
        return State(
            initial_state = initial_state or self.to_dict(),
            updaters = self.dict(self.updaters))

    def merge_updaters(self, updaters):
        ''' Merge in a new set of updaters '''

        self.updaters.merge(updaters)

    def declare_state(self, keys):
        ''' Initialize values for the given keys to zero. '''

        for key in keys:
            if not key in self.state:
                self.state[key] = 0

    def assign_values(self, values):
        ''' Assign a dict of keys and values to the state. '''
        self.state.update(values)

    def apply_update(self, update):
        ''' Apply a dict of keys and values to the state using its updaters. '''

        state_dict = self.to_dict()

        for key, value in update.items():
            # updater can be a function or a key into the updater library
            updater = self.updaters.get(key, 'delta')
            if not callable(updater):
                updater = updater_library[updater]

            self.new_state[key], other_updates = updater(
                key,
                state_dict,
                self.new_state[key],
                value)

            self.new_state.update(other_updates)

    def apply_updates(self, updates):
        ''' Apply a list of updates to the state '''

        for update in updates:
            self.apply_update(update)

    def prepare(self):
        ''' Prepares for state updates by creating new copy of existing state '''

        self.new_state = copy.deepcopy(self.state)

    def proceed(self):
        ''' Once all updates are complete, swaps out state for newly calculated state '''

        self.state = self.new_state

    def state_for(self, keys):
        ''' Get the current state of these keys as a dict of values. '''

        return {
            key: self.state[key]
            for key in keys}

    def to_dict(self):
        ''' Get the current state of all keys '''

        return copy.deepcopy(self.state)


class Process(object):
    def __init__(self, roles, parameters=None):
        ''' Declare what roles this process expects. '''

        self.roles = roles
        self.parameters = parameters or {}
        self.states = None

    def default_settings(self):
        return {}

    def assign_roles(self, states):
        '''
        Provide States for some or all of the roles this Process expects.

        Roles and States must have the same keys. '''

        self.states = states
        for role, state in self.states.items():
            state.declare_state(self.roles[role])

    def update_for(self, timestep):
        ''' Called each timestep to find the next state for this process. '''

        states = {
            role: self.states[role].state_for(self.roles[role])
            for role, values in self.roles.items()}

        return self.next_update(timestep, states)

    def next_update(self, timestep, states):
        '''
        Find the next update given the current states this process cares about.

        This is the main function a new process would override.'''

        return {
            role: {}
            for role, values in self.roles.items()}


def connect_topology(process_layers, states, topology):
    ''' Given a set of processes and states, and a description of the connections
        between them, link the roles in each process to the state they refer to.'''

    for processes in process_layers:
        for name, process in processes.items():
            connections = topology[name]
            roles = {
                role: states[key]
                for role, key in connections.items()}

            process.assign_roles(roles)

def merge_default_states(processes):
    initial_state = {}
    for process_id, process in merge_dicts(processes).items():
        settings = process.default_settings()
        default_state = settings['state']
        initial_state = deep_merge(dict(initial_state), default_state)
    return initial_state

def merge_default_updaters(processes):
    updaters = {}
    for process_id, process in merge_dicts(processes).items():
        settings = process.default_settings()
        default_updaters = settings['updaters']
        updaters = deep_merge(dict(updaters), default_updaters)
    return updaters

def deep_merge(dct, merge_dct):
    '''
    Recursive dict merge
    This mutates dct - the contents of merge_dct are added to dct (which is also returned).
    If you want to keep dct you could call it like deep_merge(dict(dct), merge_dct)'''
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            deep_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct

def default_divide_condition(compartment):
    return False

def default_divide_state(compartment):
    divided = [{}, {}]
    for state_key, state in compartment.states.items():
        left = random.randint(0, 1)
        for index in range(2):
            divided[index][state_key] = {}
            for key, value in state.to_dict().items():
                divided[index][state_key][key] = value // 2 + (
                    value % 2 if index == left else 0)

    print('divided {}'.format(divided))
    return divided

def get_compartment_timestep(process_layers):
    # get the minimum time_step from all processes
    processes = merge_dicts(process_layers)
    minimum_step = 10

    for process_id, process_object in processes.items():
        settings = process_object.default_settings()
        time_step = settings.get('time_step', 1.0)
        minimum_step = min(time_step, minimum_step)

    return minimum_step

def initialize_state(process_layers, topology, initial_state):
    processes = merge_dicts(process_layers)

    # make a dict with the compartment's default states {roles: states}
    compartment_states = {}
    compartment_updaters = {}
    for process_id, roles_map in topology.items():
        process_roles = processes[process_id].roles

        settings = processes[process_id].default_settings()
        default_process_states = settings['state']
        default_process_updaters = settings['updaters']

        for process_role, states in process_roles.items():
            compartment_role = topology[process_id][process_role]

            # initialize the default states
            default_states = default_process_states.get(process_role, {})

            # initialize the default updaters
            default_updaters = default_process_updaters.get(process_role, {})

            # update the states
            c_states = deep_merge(default_states, compartment_states.get(compartment_role, {}))
            compartment_states[compartment_role] = c_states

            # update the updaters
            c_updaters = deep_merge(default_updaters, compartment_updaters.get(compartment_role, {}))
            compartment_updaters[compartment_role] = c_updaters

    # initialize state for each compartment role
    initialized_state = {}
    for compartment_role, states in compartment_states.items():
        updaters = compartment_updaters[compartment_role]
        make_state = State(
            initial_state=deep_merge(states, dict(initial_state.get(compartment_role, {}))),
            updaters=updaters)
        initialized_state[compartment_role] = make_state

    return initialized_state

class Compartment(object):
    ''' Track a set of processes and states and the connections between them. '''

    def __init__(self, processes, states, configuration):
        ''' Given a set of processes and states, and a topology describing their
            connections, perform those connections. '''

        self.initial_time = configuration.get('initial_time', 0.0)
        self.local_time = 0.0
        self.time_step = configuration.get('time_step', 1.0)

        self.processes = processes
        self.states = states
        self.topology = configuration['topology']

        self.divide_condition = configuration.get('divide_condition', default_divide_condition)
        self.divide_state = configuration.get('divide_state', default_divide_state)

        # emitter
        if configuration.get('emitter'):
            self.emitter_keys = configuration['emitter'].get('keys')
            self.emitter = configuration['emitter'].get('object')
        else:
            emitter = emit.get_emitter({})
            self.emitter_keys = emitter.get('keys')
            self.emitter = emitter.get('object')

        connect_topology(processes, self.states, self.topology)

        # log experiment configuration
        emit_config = {
            'table': 'configuration',
            'data': {'topology': self.topology}}
        self.emitter.emit(emit_config)

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        for processes in self.processes:
            updates = {}
            for name, process in processes.items():
                update = process.update_for(timestep)
                for role, update_dict in update.items():
                    key = self.topology[name][role]
                    if not updates.get(key):
                        updates[key] = []
                    updates[key].append(update_dict)

            for key in self.states.keys():
                self.states[key].prepare()

            for key, update in updates.items():
                self.states[key].apply_updates(update)

            for key in self.states.keys():
                self.states[key].proceed()

        self.local_time += timestep

        # run emitters
        self.emit_data()

    def current_state(self):
        ''' Construct the total current state from the existing substates. '''

        return {
            key: state.to_dict()
            for key, state in self.states.items()}

    def current_parameters(self):
        return {
            name: process.parameters
            for name, process in merge_dicts(self.processes).items()}

    def time(self):
        return self.initial_time + self.local_time

    def emit_data(self):
        data = {}
        for role_key, emit_keys in self.emitter_keys.items():
            data[role_key] = self.states[role_key].state_for(emit_keys)

        data.update({
            'type': 'compartment',
            'time': self.time()})

        emit_config = {
            'table': 'history',
            'data': data}

        self.emitter.emit(emit_config)




## functions for testing
def toy_composite(config):
    '''
    a toy composite function for testing
    returns a dictionary with 'processes', 'states', and 'options'

    '''

    # toy processes
    class ToyMetabolism(Process):
        def __init__(self, initial_parameters={}):
            roles = {'pool': ['GLC', 'MASS']}
            parameters = {'mass_conversion_rate': 1}
            parameters.update(initial_parameters)

            super(ToyMetabolism, self).__init__(roles, parameters)

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
            roles = {
                'external': ['GLC'],
                'internal': ['GLC']}
            parameters = {'intake_rate': 2}
            parameters.update(initial_parameters)

            super(ToyTransport, self).__init__(roles, parameters)

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
            roles = {
                'compartment': ['MASS', 'DENSITY', 'VOLUME']}
            parameters = {}

            super(ToyDeriveVolume, self).__init__(roles, parameters)

        def next_update(self, timestep, states):
            volume = states['compartment']['MASS'] / states['compartment']['DENSITY']
            update = {
                'compartment': {'VOLUME': volume}}

            return update

    processes = [
        {'metabolism': ToyMetabolism(
            initial_parameters={
                'mass_conversion_rate': 0.5}), # example of overriding default parameters
         'transport': ToyTransport()},
        {'external_volume': ToyDeriveVolume(),
         'internal_volume': ToyDeriveVolume()}]

    def update_mass(key, state, current, new):
        return current / (current + new), {}

    # declare the states
    states = {
        'periplasm': State(
            initial_state={'GLC': 20, 'MASS': 100, 'DENSITY': 10},
            updaters={'MASS': update_mass, 'VOLUME': 'set'}),
        'cytoplasm': State(
            initial_state={'MASS': 3, 'DENSITY': 10},
            updaters={'VOLUME': 'set'})}

    # hook up the roles in each process to compartment states
    topology = {
        'metabolism': {
            'pool': 'cytoplasm'},
        'transport': {
            'external': 'periplasm',
            'internal': 'cytoplasm'},
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

    options = {
        # 'environment_role': 'environment',
        # 'exchange_role': 'exchange',
        'emitter': emitter,
        'topology': topology,
        'initial_time': 0.0}

    return {
        'processes': processes,
        'states': states,
        'options': options}

def test_compartment():
    out_dir = os.path.join('out', 'tests', 'toy_compartment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(toy_composite)

    settings = {
        'timestep': 1,
        'total_time': 10}
    plot_settings = {}

    saved_state = simulate_compartment(compartment, settings)
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)

def load_compartment(composite=toy_composite, boot_config={}):
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
    states = composite_config['states']
    options = composite_config['options']

    compartment = Compartment(processes, states, options)
    # print('current_parameters: {}'.format(compartment.current_parameters()))
    # print('current_state: {}'.format(compartment.current_state()))

    return compartment

def simulate_compartment(compartment, settings={}):
    '''
    run a compartment simulation
    '''

    timestep = settings.get('timestep', 1)
    total_time = settings.get('total_time', 10)

    # save initial state
    time = 0
    saved_state = {}
    saved_state[time] = compartment.current_state()

    # run simulation
    while time < total_time:
        time += timestep
        compartment.update(timestep)
        saved_state[time] = compartment.current_state()

    return saved_state

def simulate_with_environment(compartment, settings={}):
    '''
    run a compartment simulation with an environment.
    requires processes made for LatticeCompartment, with environment_role and exchange_role
    '''

    environment_role = settings['environment_role']
    exchange_role = settings['exchange_role']
    exchange_ids = list(compartment.states[exchange_role].keys())
    env_volume = settings['environment_volume']  # (L)
    nAvogadro = constants.N_A

    timestep = settings.get('timestep', 1)
    total_time = settings.get('total_time', 10)
    timeline = settings.get('timeline', [(total_time, {})])
    end_time = timeline[-1][0]

    environment = compartment.states.get(environment_role)
    exchange = compartment.states.get(exchange_role)

    # initialize saved_state
    saved_state = {}

    ## run simulation
    time = 0
    saved_state[time] = compartment.current_state()
    while time < end_time:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for role_id, change in change_dict.items():
                    role = compartment.states.get(role_id)
                    role.assign_values(change)
                timeline.pop(0)

        # update compartment
        compartment.update(timestep)

        ## apply exchange to environment
        delta_counts = exchange.state_for(exchange_ids)

        # convert counts to change in concentration in environemnt
        mmol_to_count = nAvogadro * env_volume * 1e-3  # (L/mmol)
        delta_concs = {mol_id: counts / mmol_to_count  for mol_id, counts in delta_counts.items()}
        environment.apply_update(delta_concs)
        # if state[environment_role][mol_id] < 0.0:  # this shouldn't be needed
        #     state[environment_role][mol_id] = 0.0

        # reset exchange
        reset_exchange = {key: 0 for key in exchange_ids}
        exchange.assign_values(reset_exchange)

        saved_state[time] = compartment.current_state()

    return saved_state

def convert_to_timeseries(sim_output):
    '''
    input:
        - saved_states (dict) with {timestep: state_dict}
    returns:
        - timeseries (dict) with timeseries in lists {'time': [], 'role1': {'state': []}}
    TODO --  currently assumes state is 1 dictionary deep. make a more general state embedding
    '''

    time_vec = list(sim_output.keys())
    initial_state = sim_output[time_vec[0]]
    timeseries = {role: {state: []
        for state, initial in states.items()}
        for role, states in initial_state.items()}
    timeseries['time'] = time_vec

    for time, all_states in sim_output.items():
        for role, states in all_states.items():
            for state_id, state in states.items():
                timeseries[role][state_id].append(state)

    return timeseries

def set_axes(ax, show_xaxis=False):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)

def plot_simulation_output(timeseries, settings={}, out_dir='out'):
    '''
    plot simulation output
        args:
        - timeseries (dict). This can be obtained from simulation output with convert_to_timeseries()
        - settings (dict) with:
            {
            'max_rows': (int) roles with more states than this number of states get wrapped into a new column
            'remove_zeros': (bool) if True, timeseries with all 0's get removed
            'remove_flat': (bool) if True, timeseries with all the same value get removed
            'skip_roles': (list) roles that won't be plotted
            'overlay': (dict) with
                {'bottom_role': 'top_role'}  roles plotted together by matching state_ids, with 'top_role' in red
            }
    TODO -- some molecules have 'inf' concentrations for practical reasons. How should these be plotted?
    '''

    skip_keys = ['time']

    # get settings
    max_rows = settings.get('max_rows', 25)
    overlay = settings.get('overlay', {})
    skip_roles = settings.get('skip_roles', [])
    remove_flat = settings.get('remove_flat', False)
    remove_zeros = settings.get('remove_zeros', False)
    show_state = settings.get('show_state', [])
    top_roles = list(overlay.values())
    bottom_roles = list(overlay.keys())

    roles = [role for role in timeseries.keys() if role not in skip_keys + skip_roles]
    time_vec = timeseries['time']

    # remove selected states
    # TODO -- plot removed_states as text
    removed_states = []
    if remove_flat:
        # find series with all the same value
        for role in roles:
            for state_id, series in timeseries[role].items():
                if series.count(series[0]) == len(series):
                    removed_states.append((role, state_id))
    elif remove_zeros:
        # find series with all zeros
        for role in roles:
            for state_id, series in timeseries[role].items():
                if all(v == 0 for v in series):
                    removed_states.append((role, state_id))

    # if specified in show_state, keep in timeseries
    for role_state in show_state:
        if role_state in removed_states:
            removed_states.remove(role_state)

    # remove from timeseries
    for (role, state_id) in removed_states:
        del timeseries[role][state_id]

    # get the number of states in each role
    n_data = [len(timeseries[key]) for key in roles if key not in top_roles]
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
    for role in roles:
        top_timeseries = {}
        if role in bottom_roles:
            # get overlay
            top_role = overlay[role]
            top_timeseries = timeseries[top_role]
        elif role in top_roles + skip_roles:
            # don't give this row its own plot
            continue

        for state_id, series in sorted(timeseries[role].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            # plot line at zero if series crosses the zero line
            if any(x == 0.0 for x in series) or \
                    (any(x < 0.0 for x in series) and any(x > 0.0 for x in series)):
                zero_line = [0 for t in time_vec]
                ax.plot(time_vec, zero_line, 'k--')

            if (role, state_id) in show_state:
                ax.plot(time_vec, series, 'indigo')
            else:
                ax.plot(time_vec, series)

            # overlay
            if state_id in top_timeseries.keys():
                ax.plot(time_vec, top_timeseries[state_id], 'm', label=top_role)
                ax.legend()

            ax.title.set_text(str(role) + ': ' + str(state_id))
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
    fig_path = os.path.join(out_dir, 'simulation')
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    test_compartment()
