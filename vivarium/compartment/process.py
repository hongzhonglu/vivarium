from __future__ import absolute_import, division, print_function

import copy
import random

import numpy as np
import logging as log

import vivarium.compartment.emitter as emit
from vivarium.utils.dict_utils import merge_dicts, deep_merge, deep_merge_check


COMPARTMENT_STATE = '__compartment_state__'

DEFAULT_TIMESTEP = 1.0

INFINITY = float('inf')
VERBOSE = False

class topologyError(Exception):
    pass

def npize(d):
    ''' Turn a dict into an ordered set of keys and values. '''

    ordered = [[key, value] for key, value in d.items()]
    keys = [key for key, _ in ordered]
    values = np.array([value for _, value in ordered], np.float64)

    return keys, values


## updater functions
# these function take in a variable key, the entire store's dict,
# the variable's current value, the variable's current update,
# and returns a new value, and other updates
def update_accumulate(key, state_dict, current_value, new_value):
    return current_value + new_value, {}

def update_set(key, state_dict, current_value, new_value):
    return new_value, {}

def update_merge(key, state_dict, current_value, new_value):
    # merge dicts, with new_value replacing any shared keys with current_value
    update = current_value.copy()
    for k, v in current_value.items():
        new = new_value.get(k)
        if isinstance(new, dict):
            update[k] = deep_merge(dict(v), new)
        else:
            update[k] = new
    return update, {}

updater_library = {
    'accumulate': update_accumulate,
    'set': update_set,
    'merge': update_merge}


## divider functions
# these functions take in a value, are return two values for each daughter
def default_divide_condition(compartment):
    return False

def divide_set(state):
    return [state, state]

def divide_split(state):
    if isinstance(state, int):
        remainder = state % 2
        half = int(state / 2)
        if random.choice([True, False]):
            return [half + remainder, half]
        else:
            return [half, half + remainder]
    elif state == float('inf') or state == 'Infinity':
        # some concentrations are considered infinite in the environment
        # an alternative option is to not divide the local environment state
        return [state, state]
    elif isinstance(state, float):
        half = state/2
        return [half, half]
    else:
        raise Exception('can not divide state {} of type {}'.format(state, type(state)))

def divide_zero(state):
    return [0, 0]

def divide_split_dict(state):
    d1 = dict(list(state.items())[len(state) // 2:])
    d2 = dict(list(state.items())[:len(state) // 2])
    return [d1, d2]

divider_library = {
    'set': divide_set,
    'split': divide_split,
    'split_dict': divide_split_dict,
    'zero': divide_zero}



KEY_TYPE = 'U31'

def keys_list(d):
    return list(d.keys())

class Store(object):
    ''' Represents a set of named values. '''

    def __init__(self, initial_state={}, schema={}):
        ''' Keys and state initialize empty, with a maximum key length of 31. '''

        self.state = copy.deepcopy(initial_state)
        self.schema = schema

        # get updaters from schema
        updaters = {}
        for state, state_schema in schema.items():
            updater = state_schema.get('updater')
            if updater:
                updaters.update({state: updater})
        self.updaters = updaters

    def schema_properties(self, keys, schema_type):
        return {key: self.schema[key].get(schema_type) for key in keys}

    def keys(self):
        return self.state.keys()

    def duplicate(self, initial_state={}):
        return Store(
            initial_state = initial_state or self.to_dict(),
            schema = self.dict(self.schema))

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

        for key, value in update.items():
            # updater can be a function or a key into the updater library
            updater = self.updaters.get(key, 'accumulate')
            if not callable(updater):
                updater = updater_library[updater]

            try:
                self.new_state[key], other_updates = updater(
                    key,
                    self.state,
                    self.new_state[key],
                    value)
            except TypeError as e:
                log.error(
                    'bad update - {}: {}, {}'.format(
                        key, value, self.new_state[key]))
                raise e

            self.new_state.update(other_updates)

    def apply_updates(self, updates):
        ''' Apply a list of updates to the state '''

        for update in updates:
            self.apply_update(update)

    def prepare(self):
        ''' Prepares for state updates by creating new copy of existing state '''

        self.new_state = self.state
        # self.new_state = copy.deepcopy(self.state)

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

        # return copy.deepcopy(self.state)
        return {
            key: value
            for key, value in self.state.items()}


class Process(object):
    def __init__(self, ports, parameters=None):
        ''' Declare what ports this process expects. '''

        self.ports = ports
        self.parameters = parameters or {}
        self.states = None

        default_timestep = self.default_settings().get('time_step', 1.0)
        self.time_step = self.parameters.get('time_step', default_timestep)

    def local_timestep(self):
        '''
        Returns the favored timestep for this process.
        Meant to be overridden in subclasses, unless 1.0 is a happy value. 
        '''

        return self.time_step

    def default_settings(self):
        return {}


    def schema_properties(self, states, schema_type):
        '''
        Requires:
            - states (dict)
            - schema_type (str)

        Returns a dictionary with {store_id: {key: schema_value}}
        for all store_ids and list of keys in states,
        with schema_value specified by schema_type
        '''

        schema = {}
        for store_id, keys in states.items():
            schema[store_id] = self.states[store_id].schema_properties(keys, schema_type)

        return schema

    def assign_ports(self, states):
        '''
        Provide States for some or all of the ports this Process expects.

        Roles and States must have the same keys. '''

        self.states = states
        for port, state in self.states.items():
            state.declare_state(self.ports[port])

    def update_for(self, timestep):
        ''' Called each timestep to find the next state for this process. '''

        states = {
            port: self.states[port].state_for(values)
            for port, values in self.ports.items()}

        return self.next_update(timestep, states)

    def parameters_for(self, parameters, key):
        ''' Return key in parameters or from self.default_parameters if not present. '''

        return parameters.get(key, self.default_parameters[key])

    def derive_defaults(self, parameters, original_key, derived_key, f):
        present = self.parameters_for(parameters, original_key)
        self.default_parameters[derived_key] = f(present)
        return self.default_parameters[derived_key]

    def next_update(self, timestep, states):
        '''
        Find the next update given the current states this process cares about.

        This is the main function a new process would override.'''

        return {
            port: {}
            for port, values in self.ports.items()}


def connect_topology(process, derivers, states, topology):
    ''' Given a set of processes and states, and a description of the connections
        between them, link the ports in each process to the state they refer to.'''
    process_layers = process + derivers
    for processes in process_layers:
        for name, process in processes.items():
            connections = topology[name]
            ports = {
                port: states[key]
                for port, key in connections.items()}
            try:
                process.assign_ports(ports)
            except:
                print('{} mismatched ports'.format(name))

def get_minimum_timestep(process_layers):
    # get the minimum time_step from all processes
    processes = merge_dicts(process_layers)
    minimum_step = 10

    for process_id, process_object in processes.items():
        settings = process_object.default_settings()
        time_step = settings.get('time_step', DEFAULT_TIMESTEP)
        minimum_step = min(time_step, minimum_step)

    return minimum_step

def get_maximum_timestep(process_layers):
    # get the minimum time_step from all processes
    processes = merge_dicts(process_layers)
    maximum_step = 0.0

    for process_id, process_object in processes.items():
        settings = process_object.default_settings()
        time_step = settings.get('time_step', DEFAULT_TIMESTEP)
        maximum_step = max(time_step, maximum_step)

    return maximum_step

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

def initialize_state(process_layers, topology, initial_state):
    schema = get_schema(process_layers, topology)
    processes = merge_dicts(process_layers)

    # make a dict with the compartment's default states {ports: states}
    compartment_states = {}
    for process_id, ports_map in topology.items():
        process_ports = processes[process_id].ports

        settings = processes[process_id].default_settings()
        default_process_states = settings['state']

        for process_port, states in process_ports.items():
            try:
                compartment_port = topology[process_id][process_port]
            except:
                raise topologyError(
                    'no "{}" port assigned to "{}" process in topology'.format(process_port, process_id))

            # initialize the default states
            default_states = default_process_states.get(process_port, {})

            # update the states
            # TODO -- make this a deep_merge_check, requires better handling of initial state conflicts
            c_states = deep_merge(default_states, compartment_states.get(compartment_port, {}))
            compartment_states[compartment_port] = c_states

    # initialize state for each compartment port
    initialized_state = {}
    for compartment_port, states in compartment_states.items():
        state_schema = schema.get(compartment_port, {})

        make_state = Store(
            initial_state=deep_merge(states, dict(initial_state.get(compartment_port, {}))),
            schema=state_schema)
        initialized_state[compartment_port] = make_state

    return initialized_state

def flatten_process_layers(process_layers):
    processes = {}
    for layer in process_layers:
        processes.update(layer)
    return processes

class Compartment(Store):
    ''' Track a set of processes and states and the connections between them. '''

    def __init__(self, processes, derivers, states, configuration):
        ''' Given a set of processes and states, and a topology describing their
            connections, perform those connections. '''

        self.processes = processes
        self.states = states
        self.derivers = derivers
        self.configuration = configuration

        self.topology = configuration['topology']
        self.initial_time = configuration.get('initial_time', 0.0)
        self.local_time = 0.0
        self.time_step = configuration.get('time_step', get_maximum_timestep(processes))

        # configure compartment state
        self.states[COMPARTMENT_STATE] = self
        self.state = {
            'processes': self.processes,
            'derivers': self.derivers}
        self.updaters = {'processes': 'set'}

        # divide condition
        self.divide_condition = configuration.get('divide_condition', default_divide_condition)

        # emitter
        emitter_config = configuration.get('emitter')
        if emitter_config is None:
            emitter = emit.get_emitter({})
            self.emitter_keys = emitter.get('keys')
            self.emitter = emitter.get('object')
        elif isinstance(emitter_config, str):
            emitter = emit.configure_emitter(
                {'emitter': {'type': emitter_config}},
                self.processes + self.derivers,
                self.topology)
            self.emitter_keys = emitter.get('keys')
            self.emitter = emitter.get('object')
        else:
            self.emitter_keys = configuration['emitter'].get('keys')
            self.emitter = configuration['emitter'].get('object')

        connect_topology(processes, derivers, self.states, self.topology)

        # log experiment configuration
        data = {
            'type': 'compartment',
            'name': configuration.get('name', 'compartment'),
            'topology': self.topology}

        emit_config = {
            'table': 'configuration',
            'data': data}
        self.emitter.emit(emit_config)

    def divide_state(self):
        daughter_states = [{}, {}]
        for port_id, state in self.states.items():
            if port_id == COMPARTMENT_STATE:
                # TODO -- copy compartment_state to each daughter???
                break
            else:
                schema = state.schema

            for state_id, value in state.to_dict().items():
                state_schema = schema.get(state_id, {})
                divide_type = state_schema.get('divide', 'split')  # default divider is 'split'
                divider = divider_library[divide_type]

                # divide the state
                divided_state = divider(value)

                for index in range(2):
                    new_state = {
                        port_id: {
                            state_id: divided_state[index]}}
                    deep_merge(daughter_states[index], new_state)

        print('divided {}'.format(daughter_states))
        return daughter_states

    def prepare(self):
        ''' Avoid creating a copy of the process objects. '''
        self.new_state = self.state

    def to_dict(self):
        ''' overriding from State '''
        return self.current_state()

    def collect_updates(self, updates, process_name, new_update):
        for port, update_store in new_update.items():
            key = self.topology[process_name][port]
            if not updates.get(key):
                updates[key] = []
            updates[key].append(update_store)
        return updates

    def run_derivers(self):

        derivers = flatten_process_layers(self.state['derivers'])

        updates = {}
        for name, process in derivers.items():
            new_update = process.update_for(0)  # timestep shouldn't influence derivers
            updates = self.collect_updates(updates, name, new_update)

        for key, update in updates.items():
            self.states[key].apply_updates(update)

    def send_updates(self, updates):
        ''' Prepare the states, apply the updates, run derivers, and proceed'''

        for key in self.states.keys():
            self.states[key].prepare()

        for key, update in updates.items():
            self.states[key].apply_updates(update)

        # run derivers after every update
        # TODO -- only run derivers if any of their states have been updated
        self.run_derivers()

        for key in self.states.keys():
            self.states[key].proceed()

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        time = 0

        # flatten all process layers into a single process dict
        processes = flatten_process_layers(self.state['processes'])

        # keep track of which processes have simulated until when
        front = {
            process_name: {
                'time': 0,
                'update': {}}
            for process_name in processes.keys()}

        while time < timestep:
            step = INFINITY

            if VERBOSE:
                for state_id in self.states:
                    print('{}: {}'.format(time, self.states[state_id].to_dict()))

            for process_name, process in processes.items():
                process_time = front[process_name]['time']

                if process_time <= time:
                    future = min(process_time + process.local_timestep(), timestep)
                    interval = future - process_time
                    update = process.update_for(interval)

                    if interval < step:
                        step = interval
                    front[process_name]['time'] = future
                    front[process_name]['update'] = update

            if step == INFINITY:
                # no processes ran, jump to next process
                next_event = timestep
                for process_name in front.keys():
                    if front[process_name]['time'] < next_event:
                        next_event = front[process_name]['time']
                time = next_event
            else:
                # at least one process ran, apply updates and continue
                future = time + step

                updates = {}
                for process_name, advance in front.items():
                    if advance['time'] <= future:
                        updates = self.collect_updates(updates, process_name, advance['update'])
                        advance['update'] = {}

                self.send_updates(updates)

                time = future

        for process_name, advance in front.items():
            assert advance['time'] == time == timestep
            assert len(advance['update']) == 0

        self.local_time += timestep

        # run emitters
        self.emit_data()

    def current_state(self):
        ''' Construct the total current state from the existing substates. '''

        return {
            key: state.to_dict()
            for key, state in self.states.items()
            if key != COMPARTMENT_STATE}

    def current_parameters(self):
        return {
            name: process.parameters
            for name, process in merge_dicts(self.state['processes']).items()}

    def time(self):
        return self.initial_time + self.local_time

    def emit_data(self):
        data = {}
        for port_key, emit_keys in self.emitter_keys.items():
            data[port_key] = self.states[port_key].state_for(emit_keys)

        data.update({
            'type': 'compartment',
            'time': self.time()})

        emit_config = {
            'table': 'history',
            'data': data}

        self.emitter.emit(emit_config)


def test_timescales():
    class Slow(Process):
        def __init__(self):
            self.timestep = 3.0
            self.ports = {
                'state': ['base']}

        def local_timestep(self):
            return self.timestep

        def next_update(self, timestep, states):
            base = states['state']['base']
            next_base = timestep * base * 0.1

            return {
                'state': {'base': next_base}}

    class Fast(Process):
        def __init__(self):
            self.timestep = 0.1
            self.ports = {
                'state': ['base', 'motion']}

        def local_timestep(self):
            return self.timestep

        def next_update(self, timestep, states):
            base = states['state']['base']
            motion = timestep * base * 0.001

            return {
                'state': {'motion': motion}}

    processes = [{
        'slow': Slow(),
        'fast': Fast()}]

    derivers = []

    states = {
        'state': Store({
            'base': 1.0,
            'motion': 0.0})}

    topology = {
        'slow': {'state': 'state'},
        'fast': {'state': 'state'}}

    emitter_config = {
            'type': 'print',
            'keys': {
                'state': ['base', 'motion']}}

    configuration = {
        'topology': topology,
        'emitter': emit.get_emitter(emitter_config)}

    compartment = Compartment(
        processes,
        derivers,
        states,
        configuration)

    compartment.update(10.0)

if __name__ == '__main__':
    test_timescales()
