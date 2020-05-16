from __future__ import absolute_import, division, print_function

import copy
import random

import numpy as np
import logging as log

from vivarium.compartment.process import (
    Process,
    divider_library)
from vivarium.utils.dict_utils import merge_dicts, deep_merge, deep_merge_check

# processes
from vivarium.processes.derive_globals import DeriveGlobals, AVOGADRO
from vivarium.processes.derive_counts import DeriveCounts
from vivarium.processes.derive_concentrations import DeriveConcs
from vivarium.processes.derive_mass import DeriveMass


deriver_library = {
    'mass': DeriveMass,
    'globals': DeriveGlobals,
    'mmol_to_counts': DeriveCounts,
    'counts_to_mmol': DeriveConcs,
}


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


def dissoc(d, removing):
    return {
        key: value
        for key, value in d.items()
        if not key in removing}
        

def select_keys(d, keys):
    return {
        key: d.get(key)
        for key in keys}


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

def always_true(x):
    return True

def identity(y):
    return y

class Store(object):
    schema_keys = set([
        '_default',
        '_updater',
        '_divider',
        '_value',
        '_divider',
        '_properties',
        '_units'])

    def __init__(self, config, parent=None):
        self.parent = parent
        self.children = {}
        self.subschema = {}
        self.properties = {}
        self.default = None
        self.updater = None
        self.value = None
        self.units = None
        self.divider = None

        self.apply_config(config)

    def apply_config(self, config):
        if '_subschema' in config:
            self.subschema = deep_merge(
                self.subschema,
                config.get('_subschema', {}))
            config = {
                key: value
                for key, value in config.items()
                if key != '_subschema'}

        if self.schema_keys & config.keys():
            # TODO (Eran) -- check if config schema match existing schema, make exceptions if mismatched
            self.default = config.get('_default', self.default)
            self.value = config.get('_value', self.value or self.default)

            self.updater = config.get('_updater', self.updater or 'accumulate')
            if isinstance(self.updater, str):
                self.updater = updater_library[self.updater]

            self.divider = config.get('_divider', self.divider)
            if isinstance(self.divider, str):
                self.divider = divider_library[self.divider]
            if isinstance(self.divider, dict) and isinstance(self.divider['divider'], str):
                self.divider['divider'] = divider_library[self.divider['divider']]

            self.properties = deep_merge(
                self.properties,
                config.get('_properties', {}))

            self.units = config.get('_units', self.units)
        else:
            self.value = None

            for key, child in config.items():
                if not key in self.children:
                    self.children[key] = Store(child, parent=self)
                else:
                    self.children[key].apply_config(child)

    def get_config(self):
        config = {}
        if self.properties:
            config['_properties'] = self.properties
        if self.subschema:
            config['_subschema'] = self.subschema

        if self.children:
            child_config = {
                key: child.get_config()
                for key, child in self.children.items()}
            config.update(child_config)
        else:
            config.update({
                '_default': self.default,
                '_value': self.value})
            if self.updater:
                config['_updater'] = self.updater
            if self.divider:
                config['_divider'] = self.divider
            if self.units:
                config['_units'] = self.units

        return config

    def get_value(self, condition=None, f=None):
        if self.children:
            if condition is None:
                condition = always_true

            if f is None:
                f = identity

            return {
                key: f(child.get_value(condition, f))
                for key, child in self.children.items()
                if condition(child)}
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

    def get_paths(self, paths):
        return {
            key: self.get_path(path)
            for key, path in paths.items()}

    def get_values(self, paths):
        return {
            key: self.get_in(path)
            for key, path in paths.items()}

    def get_in(self, path):
        return self.get_path(path).get_value()

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

    def delete_path(self, path):
        if not path:
            self.children = {}
            self.value = None
            return self
        else:
            target = self.get_path(path[:-1])
            remove = path[-1]
            if remove in target.children:
                lost = target.children[remove]
                del target.children[remove]
                return lost

    def divide_value(self):
        if self.divider:
            # divider is either a function or a dict with topology
            if isinstance(self.divider, dict):
                divider = self.divider['divider']
                topology = self.divider['topology']
                state = self.parent.get_values(topology)
                return divider(self.value, state)
            else:
                return self.divider(self.value)
        elif self.children:
            daughters = [{}, {}]
            for key, child in self.children.items():
                division = child.divide_value()
                if division:
                    for daughter, divide in zip(daughters, division):
                        daughter[key] = divide
            return daughters

    def set_value(self, value):
        if self.children:
            for child, child_value in value.items():
                if child in self.children:
                    self.children[child].set_value(child_value)
        else:
            self.value = value

    def apply_update(self, update):
        if self.children:
            topology_updates = {}

            if '_delete' in update:
                # delete a list of paths
                for path in update['_delete']:
                    self.delete_path(path)

                update = dissoc(update, ['_delete'])

            if '_generate' in update:
                # generate a list of new compartments
                for generate in update['_generate']:
                    self.generate(
                        generate['path'],
                        generate['processes'],
                        generate['topology'],
                        generate['initial_state'])
                    assoc_in(
                        topology_updates,
                        generate['path'],
                        generate['topology'])
                self.apply_subschemas()

                update = dissoc(update, '_generate')

            if '_divide' in update:
                # use dividers to find initial states for daughters
                divide = update['_divide']
                mother = divide['mother']
                daughters = divide['daughters']
                initial_state = self.children[mother].get_value(
                    condition=lambda child: not(isinstance(child.value, Process)),
                    f=lambda child: copy.deepcopy(child))
                states = self.children[mother].divide_value()
                for daughter, state in zip(daughters, states):
                    daughter_id = daughter['daughter']
                    initial_state = deep_merge(
                        initial_state,
                        state)

                    self.generate(
                        daughter['path'],
                        daughter['processes'],
                        daughter['topology'],
                        daughter['initial_state'])
                    assoc_in(
                        topology_updates,
                        daughter['path'],
                        daughter['topology'])

                    self.children[daughter_id].set_value(initial_state)

                self.delete_path((mother,))
                self.apply_subschemas()

                update = dissoc(update, '_divide')

            for key, value in update.items():
                if key in self.children:
                    child = self.children[key]
                    topology_updates = deep_merge(
                        topology_updates,
                        {key: child.apply_update(value)})

            return topology_updates

        else:
            self.value = self.updater(self.value, update)

    def child_value(self, key):
        if key in self.children:
            return self.children[key].get_value()

    def state_for(self, path, keys):
        state = self.get_path(path)
        if state is None:
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
                    self.parent = Store({})
                    self.parent.children[child_key] = self
                return self.parent.establish_path(remaining, config, child_key)
            else:
                if not step in self.children:
                    self.children[step] = Store({}, self)
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
                process_state = Store({
                    '_value': subprocess,
                    '_updater': 'set'}, self)
                self.children[key] = process_state
                for port, targets in subprocess.ports_schema().items():
                    path = subtopology[port]
                    if path:
                        initial = get_in(initial_state, path)
                        for target, schema in targets.items():
                            if target == '*':
                                glob = self.establish_path(path, {
                                    '_subschema': schema})
                                glob.apply_subschema()
                            else:
                                if initial and target in initial:
                                    schema = dict(
                                        schema,
                                        _value=initial[target])
                                subpath = tuple(path) + (target,)
                                self.establish_path(subpath, schema)
            else:
                if not key in self.children:
                    self.children[key] = Store({}, self)
                substate = initial_state.get(key, {})
                self.children[key].generate_paths(
                    subprocess,
                    subtopology,
                    substate)

    def generate(self, path, processes, topology, initial_state):
        target = self.establish_path(path, {})
        target.generate_paths(processes, topology, initial_state)


def process_derivers_config(processes, topology):
    ''' get the deriver configuration from processes' deriver_settings'''

    deriver_configs = {}
    full_deriver_topology = {}

    for process_id, process in processes.items():
        if isinstance(process, Process):
            process_settings = process.default_settings()
            deriver_setting = process_settings.get('deriver_setting', [])
            if process_id in topology:
                port_map = topology[process_id]
            else:
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
                        'global': ('..', 'global')}}
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


class Compartment(object):
    def __init__(self, config):
        self.config = config

    def generate_processes(self, config):
        return {}

    def generate_topology(self, config):
        return {}

    def generate(self, config):
        processes = self.generate_processes(config)
        topology = self.generate_topology(config)

        # add derivers
        derivers = process_derivers(processes, topology)
        processes.update(derivers['processes'])
        topology.update(derivers['topology'])

        return {
            'processes': processes,
            'topology': topology}


def generate_state(processes, topology, initial_state):
    state = Store({})
    state.generate_paths(processes, topology, initial_state)
    state.apply_subschemas()
    return state


def normalize_path(path):
    progress = []
    for step in path:
        if step == '..' and len(progress) > 0:
            progress = progress[:-1]
        else:
            progress.append(step)
    return progress


class Experiment(object):
    def __init__(self, config):
        self.processes = config['processes']
        self.topology = config['topology']
        self.initial_state = config['initial_state']

        self.state = generate_state(
            self.processes,
            self.topology,
            self.initial_state)

        self.local_time = 0.0

    def absolute_update(self, path, new_update):
        absolute = {}
        for port, update in new_update.items():
            topology = get_in(self.topology, path + (port,))
            if topology is not None:
                state_path = path[:-1] + topology
                normal_path = normalize_path(state_path)
                assoc_in(absolute, normal_path, update)
        return absolute

    def process_update(self, path, state, interval):
        process = state.value
        process_topology = get_in(self.topology, path)
        ports = process.find_states(state.parent, process_topology)
        update = process.next_update(interval, ports)
        absolute = self.absolute_update(path, update)
        return absolute

    def apply_update(self, update):
        topology_updates = self.state.apply_update(update)
        if topology_updates:
            self.topology = deep_merge(self.topology, topology_updates)

    def run_derivers(self, derivers):
        for path, deriver in derivers.items():
            # timestep shouldn't influence derivers
            update = self.process_update(path, deriver, 0)
            self.apply_update(update)

    def send_updates(self, updates, derivers=None):
        for update in updates:
            self.apply_update(update)
        if derivers is None:
            derivers = {
                path: state
                for path, state in self.state.depth()
                if state.value is not None and isinstance(state.value, Process) and state.value.is_deriver()}
        self.run_derivers(derivers)

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        time = 0

        def empty_front():
            return {
                'time': 0,
                'update': {}}

        # keep track of which processes have simulated until when
        front = {}

        while time < timestep:
            step = INFINITY

            if VERBOSE:
                for state_id in self.states:
                    print('{}: {}'.format(time, self.states[state_id].to_dict()))

            all_processes = {
                path: state
                for path, state in self.state.depth()
                if state.value is not None and isinstance(state.value, Process)}

            processes = {
                path: state
                for path, state in all_processes.items()
                if not state.value.is_deriver()}

            derivers = {
                path: state
                for path, state in all_processes.items()
                if state.value.is_deriver()}

            for path, state in processes.items():
                if not path in front:
                    front[path] = empty_front()
                process_time = front[path]['time']

                if process_time <= time:
                    process = state.value
                    future = min(process_time + process.local_timestep(), timestep)
                    interval = future - process_time
                    update = self.process_update(path, state, interval)

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

                updates = []
                for path, advance in front.items():
                    if advance['time'] <= future:
                        new_update = advance['update']
                        new_update['_path'] = path
                        updates.append(new_update)
                        advance['update'] = {}

                self.send_updates(updates, derivers)

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

    state = Store(environment_config)

    state.apply_update({})
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
