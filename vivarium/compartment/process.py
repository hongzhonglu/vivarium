from __future__ import absolute_import, division, print_function

import copy
import os

import numpy as np

import vivarium.compartment.emitter as emit
from vivarium.utils.dict_utils import merge_dicts, deep_merge, deep_merge_check


COMPARTMENT_STATE = '__compartment_state__'

class topologyError(Exception):
    pass

def npize(d):
    ''' Turn a dict into an ordered set of keys and values. '''

    ordered = [[key, value] for key, value in d.items()]
    keys = [key for key, _ in ordered]
    values = np.array([value for _, value in ordered], np.float64)

    return keys, values


# updater functions
def update_accumulate(key, state_dict, current_value, new_value):
    return current_value + new_value, {}

def update_set(key, state_dict, current_value, new_value):
    return new_value, {}

updater_library = {
    'accumulate': update_accumulate,
    'set': update_set}


# divider functions
def default_divide_condition(compartment):
    return False

def divide_set(state):
    return [state, state]

def divide_split(state):
    if isinstance(state, int):
        remainder = state % 2
        half = int(state / 2)
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

    def __init__(self, initial_state={}, updaters={}):
        ''' Keys and state initialize empty, with a maximum key length of 31. '''

        self.state = copy.deepcopy(initial_state)
        self.updaters = updaters

    def keys(self):
        return self.state.keys()

    def duplicate(self, initial_state={}):
        return Store(
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
            updater = self.updaters.get(key, 'accumulate')
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
            role: self.states[role].state_for(values)
            for role, values in self.roles.items()}

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
            try:
                process.assign_roles(roles)
            except:
                print('{} mismatched roles'.format(name))

def get_compartment_timestep(process_layers):
    # get the minimum time_step from all processes
    processes = merge_dicts(process_layers)
    minimum_step = 10

    for process_id, process_object in processes.items():
        settings = process_object.default_settings()
        time_step = settings.get('time_step', 1.0)
        minimum_step = min(time_step, minimum_step)

    return minimum_step

def initialize_state(process_layers, topology, schema, initial_state):
    processes = merge_dicts(process_layers)

    # make a dict with the compartment's default states {roles: states}
    compartment_states = {}
    compartment_updaters = {}
    for process_id, roles_map in topology.items():
        process_ports = processes[process_id].roles

        settings = processes[process_id].default_settings()
        default_process_states = settings['state']

        for process_port, states in process_ports.items():
            try:
                compartment_port = topology[process_id][process_port]
            except:
                raise topologyError(
                    'no "{}" role assigned to "{}" process in topology'.format(process_port, process_id))

            # initialize the default states
            default_states = default_process_states.get(process_port, {})

            # get updater from schema
            updaters = {}
            role_schema = schema.get(compartment_port, {})
            for state in states:
                if state in role_schema:
                    updater = role_schema[state].get('updater')
                    updaters.update({state: updater})

            # update the states
            # TODO -- make this a deep_merge_check, requires better handling of initial state conflicts
            c_states = deep_merge(default_states, compartment_states.get(compartment_port, {}))
            compartment_states[compartment_port] = c_states

            # update the updaters
            c_updaters = deep_merge_check(updaters, compartment_updaters.get(compartment_port, {}))
            compartment_updaters[compartment_port] = c_updaters

    # initialize state for each compartment role
    initialized_state = {}
    for compartment_port, states in compartment_states.items():
        updaters = compartment_updaters[compartment_port]
        make_state = Store(
            initial_state=deep_merge(states, dict(initial_state.get(compartment_port, {}))),
            updaters=updaters)
        initialized_state[compartment_port] = make_state

    return initialized_state

class Compartment(Store):
    ''' Track a set of processes and states and the connections between them. '''

    def __init__(self, processes, states, configuration):
        ''' Given a set of processes and states, and a topology describing their
            connections, perform those connections. '''

        self.initial_time = configuration.get('initial_time', 0.0)
        self.local_time = 0.0
        self.time_step = min(configuration.get('time_step', 1.0), get_compartment_timestep(processes))

        self.processes = processes
        self.states = states

        # configure compartment state
        self.states[COMPARTMENT_STATE] = self
        self.state = {'processes': self.processes}
        self.updaters = {'processes': 'set'}

        self.configuration = configuration
        self.topology = configuration['topology']
        self.schema = configuration['schema']

        self.divide_condition = configuration.get('divide_condition', default_divide_condition)

        # emitter
        emitter_type = configuration.get('emitter')
        if emitter_type is None:
            emitter = emit.get_emitter({})
            self.emitter_keys = emitter.get('keys')
            self.emitter = emitter.get('object')
        elif emitter_type == 'null':
            emitter = emit.get_emitter({'type': 'null'})
            self.emitter_keys = emitter.get('keys')
            self.emitter = emitter.get('object')
        else:
            self.emitter_keys = configuration['emitter'].get('keys')
            self.emitter = configuration['emitter'].get('object')

        connect_topology(processes, self.states, self.topology)

        # log experiment configuration
        data = {
            'type': 'compartment',
            'name': configuration.get('name', 'compartment'),
            'topology': self.topology,
            'schema': self.schema}

        emit_config = {
            'table': 'configuration',
            'data': data}
        self.emitter.emit(emit_config)


    def divide_state(self):
        daughter_states = [{}, {}]
        for role_id, state in self.states.items():
            if role_id == COMPARTMENT_STATE:
                # TODO -- copy compartment_state to each daughter???
                break

            for state_id, value in state.to_dict().items():
                if role_id in self.schema:
                    state_schema = self.schema[role_id].get(state_id, {})
                    divide_type = state_schema.get('divide', 'split')
                    divider = divider_library[divide_type]
                else:
                    # default divider is 'split'
                    divider = divider_library['split']

                # divide the state
                divided_state = divider(value)

                for index in range(2):
                    new_state = {
                        role_id: {
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

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        # use processes from state
        for processes in self.state['processes']:
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

    class ToyDeath(Process):
        def __init__(self, initial_parameters={}):
            roles = {
                'compartment': ['VOLUME'],
                'global': ['processes']}
            super(ToyDeath, self).__init__(roles, {})

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
        {'external_volume': ToyDeriveVolume(),
         'internal_volume': ToyDeriveVolume()},
        {'death': ToyDeath()}]

    def update_mass(key, state, current, new):
        return current / (current + new), {}

    # declare the states
    states = {
        'periplasm': Store(
            initial_state={'GLC': 20, 'MASS': 100, 'DENSITY': 10},
            updaters={'MASS': update_mass, 'VOLUME': 'set'}),
        'cytoplasm': Store(
            initial_state={'MASS': 3, 'DENSITY': 10},
            updaters={'VOLUME': 'set'})}

    # hook up the roles in each process to compartment states
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

    saved_state = simulate_compartment(compartment, settings)

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
    options.update({'emitter': boot_config.get('emitter')})

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


if __name__ == '__main__':
    test_compartment()
