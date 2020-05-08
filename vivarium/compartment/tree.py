from __future__ import absolute_import, division, print_function

import copy
import random

import numpy as np
import logging as log

from vivarium.compartment.process import Process

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



class State(object):
    schema_keys = set([
        '_default',
        '_updater',
        '_value',
        '_properties',
        '_units'])

    def __init__(self, config, parent=None):
        self.config = {}
        self.parent = parent
        self.apply_config(config)

    def apply_config(self, config):
        self.config.update(config)

        if self.schema_keys & self.config.keys():
            self.default = self.config.get('_default')
            self.updater = self.config.get('_updater', 'accumulate')
            if isinstance(self.updater, str):
                self.updater = updater_library[self.updater]
            self.value = self.config.get('_value', self.default)
            self.properties = self.config.get('_properties', {})
            self.units = self.config.get('_units')
            self.children = {}
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
                return child.get_path(path[1:])
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

    def establish_path(self, path, config, child_key='child'):
        if len(path) > 0:
            step = path[0]
            remaining = path[1:]

            if step == '..':
                if not self.parent:
                    self.parent = State({})
                self.parent.children[child_key] = self
                self.parent.establish_path(remaining, config, child_key)
            else:
                if not step in self.children:
                    self.children[step] = State({}, self)
                self.children[step].establish_path(remaining, config, child_key)
        else:
            self.apply_config(config)

    def generate_paths(self, processes, topology, initial_state):
        for key, subprocess in processes.items():
            subtopology = topology[key]
            if isinstance(subprocess, Process):
                for port, targets in subprocess.ports_schema().items():
                    path = subtopology[port]
                    if path:
                        initial = initial_state.get(port, {})
                        for target, schema in targets.items():
                            if target in initial:
                                schema = dict(
                                    schema,
                                    _value=initial[target])
                            subpath = path + [target]
                            # if subpath[0] == '..':
                            #     import ipdb; ipdb.set_trace()
                            self.establish_path(subpath, schema)
            else:
                if not key in self.children:
                    self.children[key] = State({}, self)
                substate = initial_state.get(key, {})
                self.children[key].generate_paths(
                    subprocess,
                    subtopology,
                    substate)


def initialize_state(processes, topology, initial_state):
    state = State({})
    state.generate_paths(processes, topology, initial_state)
    return state


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
                '_default': 0.0,
                '_updater': 'accumulate'},
            'fields': {
                (0, 1): {
                    'enzymeX': {
                        '_default': 0.0,
                        '_updater': 'set'},
                    'enzymeY': {
                        '_default': 0.0,
                        '_updater': 'set'}},
                (0, 2): {
                    'enzymeX': {
                        '_default': 0.0,
                        '_updater': 'set'},
                    'enzymeY': {
                        '_default': 0.0,
                        '_updater': 'set'}}},
            'agents': {
                '1': {
                    'location': {
                        '_default': (0, 0),
                        '_updater': 'set'},
                    'boundary': {
                        'external': {
                            '_default': 0.0,
                            '_updater': 'set'},
                        'internal': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            '_default': 0,
                            '_updater': 'accumulate'},
                        'fliA': {
                            '_default': 0,
                            '_updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            '_default': 0,
                            '_updater': 'set'},
                        'flagella': {
                            '_default': 0,
                            '_updater': 'accumulate'}}},
                '2': {
                    'location': {
                        '_default': (0, 0),
                        '_updater': 'set'},
                    'boundary': {
                        'external': {
                            '_default': 0.0,
                            '_updater': 'set'},
                        'internal': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            '_default': 0,
                            '_updater': 'accumulate'},
                        'fliA': {
                            '_default': 0,
                            '_updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            '_default': 0,
                            '_updater': 'set'},
                        'flagella': {
                            '_default': 0,
                            '_updater': 'accumulate'}}}}}}

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
