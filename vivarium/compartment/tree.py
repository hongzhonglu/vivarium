from __future__ import absolute_import, division, print_function

import copy
import random

import numpy as np
import logging as log

def get_in(d, path):
    if path:
        head = path[0]
        if head in d:
            return get_in(d[head], path[1:])
    else:
        return d


def assoc_in(d, path, value):
    if path:
        head = path[0]
        if len(path) == 1:
            d[head] = value
        else:
            if head not in d:
                d[head] = {}
            assoc_in(d[head], path[1:], value)
    else:
        value


def update_in(d, path, f):
    if path:
        head = path[0]
        if len(path) == 1:
            d[head] = f(d.get(head, None))
        else:
            if not head in d:
                d[head] = {}
            update_in(d[head], path[1:], f)


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
def update_accumulate(current_value, new_value):
    return current_value + new_value

def update_set(current_value, new_value):
    return new_value

def update_merge(current_value, new_value):
    # merge dicts, with new_value replacing any shared keys with current_value
    update = current_value.copy()
    for k, v in current_value.items():
        new = new_value.get(k)
        if isinstance(new, dict):
            update[k] = deep_merge(dict(v), new)
        else:
            update[k] = new
    return update

updater_library = {
    'accumulate': update_accumulate,
    'set': update_set,
    'merge': update_merge}



DEFAULT = '_default'

class State(object):
    def __init__(self, config, parent=None):
        self.config = config
        self.parent = parent

        if DEFAULT in self.config:
            self.default = self.config.get(DEFAULT)
            self.updater = self.config.get('updater', 'delta')
            if isinstance(self.updater, str):
                self.updater = updater_library[self.updater]
            self.value = self.default
            self.properties = self.config.get('properties', {})
            self.units = None
            self.children = None
        else:
            self.children = {
                key: State(child, self)
                for key, child in self.config.items()}

    def get_value(self):
        if self.children:
            return {
                key: child.get_value()
                for key, child in self.children.items()}
        else:
            return self.value

    def get_path(self, path):
        if len(path) > 0:
            step = path[0]
            if step == '..':
                child = self.parent
            else:
                child = self.children.get(step)

            if child:
                return child.get_in(path[1:])
            else:
                # TODO: more handling for bad paths?
                return None
        else:
            return self

    def get_in(self, path):
        return self.get_path(path).get_value()

    def apply_update(self, update):
        if self.children:
            for key, value in update.items():
                child = self.children[key]
                child.apply_update(value)
        else:
            self.value = self.updater(self.value, update)

    def get_template(self, template):
        '''
        Pass in a template dict with None for each value you want to
        retrieve from the tree!
        '''

        state = {}
        for key, value in template.items():
            child = self.children[key]
            if value is None:
                state[key] = child.get_value()
            else:
                state[key] = child.get_template(value)
        return state

    def state_for(self, path, keys):
        state = self.get_path(path)
        return {
            key: state.children[key].get_value()
            for key in keys}

    def depth(self, path=()):
        base = [(path, self)]
        for key, child in self.children.items():
            down = tuple(path + (key,))
            base += child.depth(down)
        return base

    def processes(self, path=()):
        return {
            path: state
            for path, state in self.depth()
            if state.value and isinstance(state.value, Process)}

class Process(object):
    def __init__(
            self,
            ports,
            parameters=None):
        ''' Declare what ports this process expects. '''

        self.ports = ports
        self.parameters = parameters or {}
        self.local_time = self.parameters.get('time', 0)
        self.states = None

        default_timestep = self.default_settings().get('time_step', 1.0)
        self.time_step = self.parameters.get('time_step', default_timestep)

        # set agent_id
        self.agent_id = parameters.get('agent_id')

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

    def find_states(self, tree, topology):
        return {
            port: tree.state_for(topology[port], keys)
            for port, keys in self.ports}

    def next_update(self, timestep, states):
        '''
        Find the next update given the current states this process cares about.

        This is the main function a new process would override.'''

        return {
            port: {}
            for port, values in self.ports.items()}


def connect_topology(processes, derivers, states, topology):
    '''
    Given a set of processes and states, and a description of the connections
    between them, link the ports in each process to the state they refer to.
    '''

    all_processes = processes.copy()
    all_processes.update(derivers)
    for name, process in all_processes.items():
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

def get_maximum_timestep(processes):
    # get the minimum time_step from all processes
    maximum_step = 0.0

    for process_id, process_object in processes.items():
        settings = process_object.default_settings()
        time_step = settings.get('time_step', DEFAULT_TIMESTEP)
        maximum_step = max(time_step, maximum_step)

    return maximum_step

def get_schema(processes, topology):
    schema = {}
    for process_id, process in processes.items():
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

def initialize_state(processes, topology, initial_state):
    schema = get_schema(processes, topology)

    # make a dict with the compartment's default states {ports: states}
    compartment_states = {}
    for process_id, ports_map in topology.items():
        process_ports = processes[process_id].ports

        settings = processes[process_id].default_settings()
        default_process_states = settings['state']

        for process_port, states in process_ports.items():
            try:
                store_id = topology[process_id][process_port]
            except:
                raise topologyError(
                    'no "{}" port assigned to "{}" process in topology'.format(process_port, process_id))

            # initialize the default states
            default_states = default_process_states.get(process_port, {})

            # update the states
            # TODO -- make this a deep_merge_check, requires better handling of initial state conflicts
            c_states = deep_merge(default_states, compartment_states.get(store_id, {}))
            compartment_states[store_id] = c_states

    # initialize state for each compartment port
    initialized_state = {}
    for store_id, states in compartment_states.items():
        state_schema = schema.get(store_id, {})

        make_state = Store(
            initial_state=deep_merge(states, dict(initial_state.get(store_id, {}))),
            schema=state_schema)
        initialized_state[store_id] = make_state

    return initialized_state


def append_update(existing_updates, new_update):
    if existing_updates is None:
        existing_updates = []
    existing_updates.append(new_update)


class Experiment(object):
    def __init__(self, config):
        self.processes = config['processes']
        self.topology = config['topology']
        self.state = generate_state(self.processes, self.topology)

    def collect_updates(self, updates, path, new_update):
        for port, update in new_update.items():
            get_in(self.topology, path)
            state_path = self.topology[port]
            update_in(updates, state_path, lambda x: append_update(x, update))
        return updates

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        time = 0

        def empty_front():
            return {
                'time': 0,
                'update': {}}

        # keep track of which processes have simulated until when
        front = {
            process_name: empty_front()
            for process_name in processes.keys()}

        while time < timestep:
            step = INFINITY

            if VERBOSE:
                for state_id in self.states:
                    print('{}: {}'.format(time, self.states[state_id].to_dict()))

            processes = {
                path: state
                for path, state in states.depth()
                if state.value and isinstance(state.value, Process)}

            for path, state in processes:
                if not path in front:
                    front[path] = empty_front()
                process_time = front[path]['time']

                if process_time <= time:
                    future = min(process_time + process.local_timestep(), timestep)
                    interval = future - process_time
                    process = state.value
                    process_topology = get_in(self.topology, path)
                    ports = process.find_states(state.parent, process_topology)
                    update = process.update_for(interval, ports)

                    if interval < step:
                        step = interval
                    front[path]['time'] = future
                    front[path]['update'] = update

            if step == INFINITY:
                # no processes ran, jump to next process
                next_event = timestep
                for process_name in front.keys():
                    if front[path]['time'] < next_event:
                        next_event = front[path]['time']
                time = next_event
            else:
                # at least one process ran, apply updates and continue
                future = time + step

                updates = {}
                for path, advance in front.items():
                    if advance['time'] <= future:
                        updates = self.collect_updates(updates, path, advance['update'])
                        advance['update'] = {}

                self.send_updates(updates)

                time = future

        for process_name, advance in front.items():
            assert advance['time'] == time == timestep
            assert len(advance['update']) == 0

        self.local_time += timestep

        # # run emitters
        # self.emit_data()
def test_recursive_store():
    environment_config = {
        'environment': {
            'temperature': {
                DEFAULT: 0.0,
                'updater': 'accumulate'},
            'fields': {
                (0, 1): {
                    'enzymeX': {
                        DEFAULT: 0.0,
                        'updater': 'set'},
                    'enzymeY': {
                        DEFAULT: 0.0,
                        'updater': 'set'}},
                (0, 2): {
                    'enzymeX': {
                        DEFAULT: 0.0,
                        'updater': 'set'},
                    'enzymeY': {
                        DEFAULT: 0.0,
                        'updater': 'set'}}},
            'agents': {
                '1': {
                    'location': {
                        DEFAULT: (0, 0),
                        'updater': 'set'},
                    'boundary': {
                        'external': {
                            DEFAULT: 0.0,
                            'updater': 'set'},
                        'internal': {
                            DEFAULT: 0.0,
                            'updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            DEFAULT: 0,
                            'updater': 'accumulate'},
                        'fliA': {
                            DEFAULT: 0,
                            'updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            DEFAULT: 0,
                            'updater': 'set'},
                        'flagella': {
                            DEFAULT: 0,
                            'updater': 'accumulate'}}},
                '2': {
                    'location': {
                        DEFAULT: (0, 0),
                        'updater': 'set'},
                    'boundary': {
                        'external': {
                            DEFAULT: 0.0,
                            'updater': 'set'},
                        'internal': {
                            DEFAULT: 0.0,
                            'updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            DEFAULT: 0,
                            'updater': 'accumulate'},
                        'fliA': {
                            DEFAULT: 0,
                            'updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            DEFAULT: 0,
                            'updater': 'set'},
                        'flagella': {
                            DEFAULT: 0,
                            'updater': 'accumulate'}}}}}}

    state = State(environment_config)

    import ipdb; ipdb.set_trace()

    state.apply_updates({})
    state.state_for({})

def test_in():
    blank = {}
    path = ['where', 'are', 'we']
    assoc_in(blank, path, 5)
    print(blank)
    print(get_in(blank, path))
    update_in(blank, path, lambda x: x + 6)
    print(blank)


if __name__ == '__main__':
    test_recursive_store()
    test_in()
