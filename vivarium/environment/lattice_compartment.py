from __future__ import absolute_import, division, print_function

import uuid

from vivarium.actor.process import Compartment, initialize_state, get_compartment_timestep
from vivarium.actor.emitter import get_emitter
from vivarium.actor.inner import Simulation

DEFAULT_COLOR = [color/255 for color in [153, 204, 255]]


# TODO -- remove these functions once all key manipulation is gone
def add_str_in_list(molecule_ids, key_str):
   return [mol_id + key_str for mol_id in molecule_ids]

def remove_str_in_list(molecule_ids, key_str):
   return [mol_id.replace(key_str, '') for mol_id in molecule_ids]

def add_str_to_keys(dct, key_str):
   ''' convert dictionary keys by adding key_str'''
   new_dct = {}
   for key, value in dct.items():
       if key_str in key:
           new_dct[key] = value
       else:
           new_dct[key + key_str] = value
   return new_dct



class LatticeCompartment(Compartment, Simulation):
    '''
    - environment_role holds the local concentrations from the external environment,
    and is updated at each exchange timestep
    - exchange_role holds accumulated molecules counts over the exchange timestep,
    and passes them to the environment upon exchange.
    '''
    def __init__(self, processes, states, configuration):
        self.color = DEFAULT_COLOR
        self.exchange_role = configuration.get('exchange_role', '')
        self.environment_role = configuration.get('environment_role', '')

        # set up exchange with lattice
        self.exchange_ids = []
        if self.exchange_role in states.keys():
            exchange_state = states[self.exchange_role]
            self.exchange_ids = exchange_state.keys()

        # find roles that contain volume, motile_force, division
        self.volume_role = False
        self.motile_role = False
        self.division_role = False
        for role, state in states.items():
            state_ids = state.keys()
            if 'volume' in state_ids:
                self.volume_role = role
            if all(item in state_ids for item in ['motile_force', 'motile_torque']):
                self.motile_role = role
            if 'division' in state_ids:
                self.division_role = role

        super(LatticeCompartment, self).__init__(processes, states, configuration)

    def run_incremental(self, run_until):
        while self.time() < run_until:
            self.update(self.time_step)

    def apply_outer_update(self, update):
        self.last_update = update
        env_keys = update['concentrations'].keys()
        environment = self.states.get(self.environment_role)

        if len(self.exchange_ids) > 0:
            # update only the states defined in both exchange and the external environment
            exchange = self.states.get(self.exchange_role)
            local_env_keys = self.exchange_ids
            exchange.assign_values({key: 0 for key in self.exchange_ids})  # reset exchange
        elif environment:
            local_env_keys = environment.keys()
        else:
            local_env_keys = []

        if environment:
            local_environment = {key : update['concentrations'][key]
                for key in local_env_keys if key in env_keys}
            environment.assign_values(local_environment)

    def generate_daughters(self):
        states = self.divide_state(self)
        volume = states[0][self.volume_role]['volume']  # TODO -- same volume for both daughters?

        return [
            dict(
                id=str(uuid.uuid1()),
                volume=volume,  # daughter_state[self.volume_role]['volume'],
                boot_config=dict(
                    initial_time=self.time(),
                    initial_state=daughter_state,
                    volume=volume))  # daughter_state[self.volume_role]['volume']))
            for daughter_state in states]

    def generate_inner_update(self):
        values = {}
        volume = 1.0
        motile_force = [0.0, 0.0]

        exchange = self.states.get(self.exchange_role)
        if exchange:
            environment_change = exchange.state_for(self.exchange_ids)
        else:
            environment_change = {}

        if self.volume_role:
            volume_state = self.states[self.volume_role].state_for(['volume'])
            volume = volume_state['volume']

        if self.motile_role:
            motile_state = self.states.get(self.motile_role)
            forces = motile_state.state_for(['motile_force', 'motile_torque'])
            motile_force = [
                forces['motile_force'],
                forces['motile_torque']]

        if self.divide_condition(self):
            values['division'] = self.generate_daughters()

            # emit phylogeny info
            daughters = {'daughters': [daughter['id'] for daughter in values['division']]}
            self.emitter.emit({
                'table': 'phylogeny',
                'data': daughters})

        values.update({
            'volume': volume,
            'motile_force': motile_force,
            'color': self.color,
            'environment_change': environment_change})

        return values


def generate_lattice_compartment(process, config):

    # get the process' default settings
    default_settings = process.default_settings()
    process_id = default_settings.get('process_id', 'process')

    # declare the processes layers (with a single layer)
    processes_layers = [{process_id: process}]

    # make a simple topology mapping 'role' to 'role'
    process_roles = process.roles.keys()
    topology = {process_id: {role: role for role in process_roles}}

    # initialize the states for each role
    states = initialize_state(processes_layers, topology, config.get('initial_state', {}))

    # get the time step
    time_step = get_compartment_timestep(processes_layers)

    # configure the emitter
    emitter_config = config.get('emitter', {})
    emitter_config['keys'] = default_settings['emitter_keys']
    emitter_config['experiment_id'] = config.get('experiment_id')
    emitter_config['simulation_id'] = config.get('simulation_id')
    emitter = get_emitter(emitter_config)

    options = {
        'emitter': emitter,
        'initial_time': config.get('initial_time', 0.0),
        'time_step': time_step,
        'exchange_role': 'exchange',  # TODO -- get this state id from a default_config() function in the process
        'environment_role': 'external',  # TODO -- get this state id from a default_config() function in the process
        'topology': topology,
    }

    # create the lattice compartment
    return LatticeCompartment(processes_layers, states, options)
