from __future__ import absolute_import, division, print_function

import copy
import random

import numpy as np
import logging as log

from vivarium.compartment.process import Process
from vivarium.compartment.composition import deriver_library
from vivarium.utils.dict_utils import merge_dicts, deep_merge, deep_merge_check

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



def schema_for(port, keys, initial_state, default=0.0, updater='accumulate'):
    return {
        key: {
            '_default': initial_state.get(
                port, {}).get(key, default),
            '_updater': updater}
        for key in keys}

class State(object):
    schema_keys = set([
        '_default',
        '_updater',
        '_value',
        '_properties',
        '_units',
        '_subschema'])

    def __init__(self, config, parent=None):
        self.config = {}
        self.parent = parent
        self.children = {}
        self.subschema = {}
        self.properties = {}

        self.apply_config(config)

    def apply_config(self, config):
        self.config.update(config)

        if self.schema_keys & self.config.keys():
            self.default = self.config.get('_default')
            self.updater = self.config.get('_updater', 'accumulate')
            if isinstance(self.updater, str):
                self.updater = updater_library[self.updater]
            self.value = self.config.get('_value', self.default)
            self.properties = deep_merge(
                self.properties,
                self.config.get('_properties', {}))
            self.units = self.config.get('_units')
            self.subschema = deep_merge(
                self.subschema,
                self.config.get('_subschema', {}))
        else:
            self.value = None
            for key, child in self.config.items():
                if not key in self.children:
                    self.children[key] = State(child, self)
                else:
                    self.children[key].apply_config(child)

    def get_value(self):
        if self.children:
            return {
                key: child.get_value()
                for key, child in self.children.items()}
        else:
            return self.value

    def get_path(self, path):
        if path:
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

    def child_value(self, key):
        if key in self.children:
            return self.children[key].get_value()

    def state_for(self, path, keys):
        state = self.get_path(path)
        if state is None:
            print('nil path: {} -- {}'.format(path, keys))
            return {}
        elif keys and keys[0] == '*':
            return state.get_value()
        else:
            return {
                key: state.child_value(key)
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
                return self.parent.establish_path(remaining, config, child_key)
            else:
                if not step in self.children:
                    self.children[step] = State({}, self)
                return self.children[step].establish_path(remaining, config, child_key)
        else:
            self.apply_config(config)
            return self

    def apply_subschema(self, subschema=None):
        if subschema is None:
            subschema = self.subschema
        for child_key, child in self.children.items():
            child.apply_config(subschema)

    def apply_subschemas(self):
        if self.subschema:
            self.apply_subschema()
        for child in self.children.values():
            child.apply_subschemas()

    def generate_paths(self, processes, topology, initial_state):
        for key, subprocess in processes.items():
            subtopology = topology[key]
            if isinstance(subprocess, Process):
                process_state = State({'_value': subprocess}, self)
                self.children[key] = process_state
                for port, targets in subprocess.ports_schema().items():
                    path = subtopology[port]
                    if path:
                        initial = initial_state.get(port, {})
                        for target, schema in targets.items():
                            if target == '*':
                                glob = self.establish_path(path, {
                                    '_subschema': schema})
                                glob.apply_subschema()
                            else:
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


def process_derivers_config(processes, topology):
    ''' get the deriver configuration from processes' deriver_settings'''

    deriver_configs = {}
    full_deriver_topology = {}

    for process_id, process in processes.items():
        if isinstance(process, Process):
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
                global_port = setting.get('global_port', ['global'])

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
                        # TODO
                        'global': ['..', 'global']}}
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

def process_derivers(processes, topology, deriver_config={}):
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
    process_config = process_derivers_config(processes, topology)

    deriver_configs = process_config['deriver_configs']
    deriver_topology = process_config['deriver_topology']

    # update deriver_configs
    deriver_configs = deep_merge(deriver_configs, deriver_config)

    # update topology based on deriver_config
    for process_id, config in deriver_configs.items():
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

    return {
        'processes': deriver_processes,
        'topology': deriver_topology}



def generate_state(processes, topology, initial_state):
    state = State({})
    state.generate_paths(processes, topology, initial_state)
    state.apply_subschemas()
    return state


def append_update(existing_updates, new_update):
    if existing_updates is None:
        existing_updates = []
    existing_updates.append(new_update)
    return existing_updates

def normalize_path():
    pass


class Experiment(object):
    def __init__(self, config):
        self.processes = config['processes']
        self.topology = config['topology']
        self.initial_state = config['initial_state']
        self.state = generate_state(
            self.processes,
            self.topology,
            self.initial_state)

    def collect_updates(self, updates, path, new_update):
        for port, update in new_update.items():
            topology = get_in(self.topology, path)
            state_path = list(path[:-1]) + topology[port]
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
        front = {}
        # front = {
        #     process_name: empty_front()
        #     for process_name in self.processes.keys()}

        while time < timestep:
            step = INFINITY

            if VERBOSE:
                for state_id in self.states:
                    print('{}: {}'.format(time, self.states[state_id].to_dict()))

            processes = {
                path: state
                for path, state in self.state.depth()
                if state.value is not None and isinstance(state.value, Process)}

            for path, state in processes.items():
                if not path in front:
                    front[path] = empty_front()
                process_time = front[path]['time']

                if process_time <= time:
                    process = state.value
                    future = min(process_time + process.local_timestep(), timestep)
                    interval = future - process_time
                    process_topology = get_in(self.topology, path)
                    ports = process.find_states(state.parent, process_topology)
                    update = process.next_update(interval, ports)

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

                import ipdb; ipdb.set_trace()

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
