from __future__ import absolute_import, division, print_function

import uuid

from lens.actor.process import Compartment, State, deep_merge
from lens.actor.emitter import get_emitter
from lens.actor.inner import Simulation

DEFAULT_COLOR = [color/255 for color in [153, 204, 255]]

exchange_key = '__exchange'  # TODO -- this is declared in multiple locations


def add_str_in_list(molecule_ids, key_str):
   return [mol_id + key_str for mol_id in molecule_ids]

def remove_str_in_list(molecule_ids, key_str):
   return [mol_id.replace(key_str, '') for mol_id in molecule_ids]

def add_str_to_keys(dct, key_str):
   ''' convert dictionary keys by adding key_str'''
   new_dct = {}
   for key, value in dct.iteritems():
       if key_str in key:
           new_dct[key] = value
       else:
           new_dct[key + key_str] = value
   return new_dct

def remove_str_from_keys(dct, key_str):
    ''' convert dictionary keys by removing key_str'''
    new_dct = {}
    for key, value in dct.iteritems():
        if key_str in key:
            new_key = key.replace(key_str, '')
            new_dct[new_key] = value
        else:
            new_dct[key] = value
    return new_dct



class LatticeCompartment(Compartment, Simulation):
    def __init__(self, processes, states, configuration):
        self.environment = configuration['environment']
        self.compartment = configuration['compartment']
        self.exchange_key = configuration['exchange_key']
        self.environment_ids = configuration['environment_ids']
        self.exchange_ids = [state_id + self.exchange_key for state_id in self.environment_ids]
        self.configuration = configuration
        self.color = DEFAULT_COLOR

        super(LatticeCompartment, self).__init__(processes, states, configuration)

    def run_incremental(self, run_until):
        while self.time() < run_until:
            self.update(self.time_step)

    def apply_outer_update(self, update):
        self.last_update = update
        environment = self.states.get(self.environment)
        if environment:
            # update only the keys defined in environment
            env_keys = environment.keys
            local_environment = {key : update['concentrations'][key] for key in env_keys}
            environment.assign_values(local_environment)
            environment.assign_values({key: 0 for key in self.exchange_ids})  # reset exchange

    def generate_daughters(self):
        states = self.divide_state(self)
        volume = states[0][self.compartment]['volume']

        return [
            dict(
                id=str(uuid.uuid1()),
                volume=volume,
                boot_config=dict(
                    initial_time=self.time(),
                    initial_state=daughter_state,
                    volume=volume))
            for daughter_state in states]

    def generate_inner_update(self):
        environment = self.states.get(self.environment)
        if environment:
            changes = environment.state_for(self.exchange_ids)
            environment_change = {
                mol_id.replace(self.exchange_key, ''): value
                for mol_id, value in changes.items()}
        else:
            environment_change = {}

        state = self.states[self.compartment]
        values = state.state_for(['volume'])

        # check if state has motile_force and motile_torque
        if 'motile_force' in state.keys and 'motile_torque' in state.keys:
            forces = state.state_for(['motile_force', 'motile_torque'])
            motile_force = [
                forces['motile_force'],
                forces['motile_torque']]
        else:
            motile_force = [0.0, 0.0]

        if self.divide_condition(self):
            values['division'] = self.generate_daughters()

            # emit phylogeny info
            daughters = {'daughters': [daughter['id'] for daughter in values['division']]}
            self.emitter.emit({
                'table': 'phylogeny',
                'data': daughters})

        values.update({
            'motile_force': motile_force, # TODO -- get motile_force from compartment state
            'color': self.color,
            'environment_change': environment_change})

        return values

def generate_lattice_compartment(process, config):
    # declare the processes
    processes = {
        'process': process}

    # initialize states
    default_states = process.default_state()
    default_updaters = process.default_updaters()
    initial_state = config.get('initial_state', {})
    initial_time = config.get('initial_time', 0.0)

    # initialize keys for accumulate_delta updater
    # assumes all environmental updates have been set as `accumulate` updaters
    environment_ids = []
    initial_exchanges = {role: {} for role in process.roles.keys()}
    for role, state_ids in default_updaters.iteritems():
        for state_id, updater in state_ids.iteritems():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})

    default_states = deep_merge(default_states, initial_exchanges)

    states = {
        role: State(
            initial_state=deep_merge(
                default_states.get(role, {}),
                dict(initial_state.get(role, {}))),
            updaters=default_updaters.get(role, {}))
            for role in process.roles.keys()}

    # configure the states to the roles for each process
    topology = {
        'process': {
            role: role
            for role in process.roles.keys()}}

    # configure emitter
    emitter_config = config.get('emitter', {})
    emitter_config['keys'] = process.default_emitter_keys()
    emitter_config['experiment_id'] = config.get('experiment_id')
    emitter_config['simulation_id'] = config.get('simulation_id')
    emitter = get_emitter(emitter_config)

    options = {
        'topology': topology,
        'emitter': emitter,
        'initial_time': initial_time,
        'environment_ids': environment_ids,
        'exchange_key': exchange_key,
        'environment': config.get('environment', 'external'),
        'compartment': config.get('compartment', 'internal')}

    options.update(config['compartment_options'])

    # create the compartment
    return LatticeCompartment(processes, states, options)
