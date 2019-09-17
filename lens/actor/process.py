from __future__ import absolute_import, division, print_function

import numpy as np

def npize(d):
    ''' Turn a dict into an ordered set of keys and values. '''

    ordered = [[key, value] for key, value in d.iteritems()]
    keys = [key for key, _ in ordered]
    values = np.array([value for _, value in ordered], np.float64)

    return keys, values


KEY_TYPE = 'U31'


class State(object):
    ''' Represents a set of named values. '''

    def __init__(self, initial_state={}):
        ''' Keys and state initialize empty, with a maximum key length of 31. '''

        self.keys = np.array([], dtype=KEY_TYPE) # maximum key length
        self.state = np.array([], dtype=np.float64)

        self.initialize_state(initial_state)

    def initialize_state(self, initial):
        ''' Provide an initial value for any keys in this dict. '''

        self.declare_state(initial.keys())
        self.apply_delta(initial)

    def sort_keys(self):
        ''' re-sort keys after adding some '''

        sort = np.argsort(self.keys)
        self.keys = self.keys[sort]
        self.state = self.state[sort]

    def declare_state(self, keys):
        ''' Initialize values for the given keys to zero. '''

        existing_keys = np.isin(keys, self.keys)
        novel = np.array(keys, dtype=KEY_TYPE)[~existing_keys]

        self.keys = np.concatenate([self.keys, novel])
        self.state = np.concatenate([self.state, np.zeros(novel.shape)])
        self.sort_keys()

    def index_for(self, keys):
        return np.searchsorted(self.keys, keys)

    def assign_values(self, values):
        ''' Assign a dict of keys and values to the state. '''

        keys, values = npize(values)
        index = self.index_for(keys)
        self.state[index] = values

    def apply_delta(self, delta):
        ''' Apply a dict of keys and deltas to the state. '''

        keys, values = npize(delta)
        index = self.index_for(keys)
        self.state[index] += values

    def apply_deltas(self, deltas):
        ''' Apply a list of deltas to the state. '''

        for delta in deltas:
            self.apply_delta(delta)

    def state_for(self, keys):
        ''' Get the current state of these keys as a dict of values. '''

        index = self.index_for(keys)
        return dict(zip(keys, self.state[index]))

    def to_dict(self):
        ''' Get the current state of all keys '''

        return {
            self.keys[index]: self.state[index]
            for index in range(self.keys.shape[0])}


class Process(object):
    def __init__(self, roles, parameters=None, deriver=False):
        ''' Declare what roles this process expects. '''

        self.roles = roles
        self.parameters = parameters or {}
        self.states = None
        self.deriver = deriver

    def default_state(self):
        return {}

    def default_emitter_keys(self):
        return {'external': [], 'internal': []}

    # def default_parameters(self):
    #     return {}

    def assign_roles(self, states):
        '''
        Provide States for some or all of the roles this Process expects.

        Roles and States must have the same keys. '''

        self.states = states
        for role, state in self.states.iteritems():
            state.declare_state(self.roles[role])

    def update_for(self, timestep):
        ''' Called each timestep to find the next state for this process. '''

        states = {
            role: self.states[role].state_for(self.roles[role])
            for role, values in self.roles.iteritems()}

        return self.next_update(timestep, states)

    def next_update(self, timestep, states):
        '''
        Find the next update given the current states this process cares about.

        This is the main function a new process would override.'''

        return {
            role: {}
            for role, values in self.roles.iteritems()}


def connect_topology(processes, states, topology):
    ''' Given a set of processes and states, and a description of the connections
        between them, link the roles in each process to the state they refer to.'''

    for name, process in processes.iteritems():
        connections = topology[name]
        roles = {
            role: states[key]
            for role, key in connections.iteritems()}

        process.assign_roles(roles)


class Compartment(object):
    ''' Track a set of processes and states and the connections between them. '''

    def __init__(self, processes, states, configuration):
        ''' Given a set of processes and states, and a topology describing their
            connections, perform those connections. '''

        self.local_time = 0.0
        self.time_step = configuration.get('time_step', 1.0)

        self.derivers = {
            name: process
            for name, process in processes.iteritems()
            if process.deriver}

        self.processes = {
            name: process
            for name, process in processes.iteritems()
            if not process.deriver}

        self.states = states
        self.topology = configuration['topology']

        # emitter
        self.emitter_keys = configuration['emitter'].get('keys')
        self.emitter = configuration['emitter'].get('object')

        connect_topology(processes, self.states, self.topology)
        self.run_derivers(0)

    def run_derivers(self, timestep):
        ''' Run each deriver process to set up state for subsequent processes. '''

        for name, deriver in self.derivers.iteritems():
            all_updates = deriver.update_for(timestep)
            for role, update in all_updates.iteritems():
                key = self.topology[name][role]
                self.states[key].assign_values(update)

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        updates = {}
        for name, process in self.processes.iteritems():
            update = process.update_for(timestep)
            for role, deltas in update.iteritems():
                key = self.topology[name][role]

                if not updates.get(key):
                    updates[key] = []
                updates[key].append(deltas)

        for key, update in updates.iteritems():
            self.states[key].apply_deltas(update)

        self.run_derivers(timestep)

        self.local_time += timestep

        # run emitters
        self.emit_data()

    def current_state(self):
        ''' Construct the total current state from the existing substates. '''

        return {
            key: state.to_dict()
            for key, state in self.states.iteritems()}

    def current_parameters(self):
        return {
            name: process.parameters
            for name, process in self.processes.iteritems()}

    def time(self):
        return self.local_time

    def emit_data(self):
        data = {}
        for role_key, emit_keys in self.emitter_keys.iteritems():
            data[role_key] = self.states[role_key].state_for(emit_keys)
        data['time'] = self.time()

        self.emitter.emit(data)

def test_compartment():
    # simplest possible metabolism
    class Metabolism(Process):
        def __init__(self, initial_parameters={}):
            roles = {'pool': ['GLC', 'MASS']}
            parameters = {'mass_conversion_rate': 1}
            parameters.update(initial_parameters)

            super(Metabolism, self).__init__(roles, parameters)

        def next_update(self, timestep, states):
            update = {}
            glucose_required = timestep / self.parameters['mass_conversion_rate']
            if states['pool']['GLC'] >= glucose_required:
                update = {
                    'pool': {
                        'GLC': -2,
                        'MASS': 1}}

            return update

    # simplest possible transport
    class Transport(Process):
        def __init__(self, initial_parameters={}):
            roles = {
                'external': ['GLC'],
                'internal': ['GLC']}
            parameters = {'intake_rate': 2}
            parameters.update(initial_parameters)

            super(Transport, self).__init__(roles, parameters)

        def next_update(self, timestep, states):
            update = {}
            intake = timestep * self.parameters['intake_rate']
            if states['external']['GLC'] >= intake:
                update = {
                    'external': {'GLC': -2},
                    'internal': {'GLC': 2}}

            return update

    class DeriveVolume(Process):
        def __init__(self, initial_parameters={}):
            roles = {
                'compartment': ['MASS', 'DENSITY', 'VOLUME']}
            parameters = {}

            super(DeriveVolume, self).__init__(roles, parameters, deriver=True)

        def next_update(self, timestep, states):
            volume = states['compartment']['MASS'] / states['compartment']['DENSITY']
            update = {
                'compartment': {'VOLUME': volume}}

            return update

    # declare the processes
    processes = {
        'metabolism': Metabolism(initial_parameters={'mass_conversion_rate': 0.5}), # example of overriding default parameters
        'transport': Transport(),
        'external_volume': DeriveVolume(),
        'internal_volume': DeriveVolume()}

    # declare the states
    states = {
        'periplasm': State({'GLC': 20, 'MASS': 100, 'DENSITY': 10}),
        'cytoplasm': State({'MASS': 3, 'DENSITY': 10})}

    # hook up the states to the roles in each process
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

    # create the compartment (which automatically hooks everything up)
    compartment = Compartment(processes, states, topology)
    timestep = 1

    # print initial parameters and state
    print(compartment.current_parameters())
    print(compartment.current_state())

    for steps in np.arange(13):
        # run the simulation
        compartment.update(timestep)
        print(compartment.current_state())


if __name__ == '__main__':
    test_compartment()
