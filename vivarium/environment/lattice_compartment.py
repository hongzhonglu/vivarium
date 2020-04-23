from __future__ import absolute_import, division, print_function

import uuid

from vivarium.compartment.process import Compartment, initialize_state, get_minimum_timestep, Process, Store
from vivarium.compartment.emitter import get_emitter
from vivarium.actor.inner import Simulation
from vivarium.compartment.composition import get_derivers


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
    - environment_port holds the local concentrations from the external environment,
    and is updated at each exchange timestep
    - exchange_port holds accumulated molecules counts over the exchange timestep,
    and passes them to the environment upon exchange.
    '''
    def __init__(self, processes, derivers, states, configuration):
        self.exchange_port = configuration.get('exchange_port', '')
        self.environment_port = configuration.get('environment_port', '')

        # set up exchange with lattice
        self.exchange_ids = []
        if self.exchange_port in states.keys():
            exchange_state = states[self.exchange_port]
            self.exchange_ids = exchange_state.keys()

        # find ports that contain volume, motile_force, division
        self.volume_port = False
        self.motile_port = False
        self.division_port = False
        for port, state in states.items():
            state_ids = state.keys()
            if 'volume' in state_ids:
                self.volume_port = port
            if all(item in state_ids for item in ['motile_force', 'motile_torque']):
                self.motile_port = port
            if 'division' in state_ids:
                self.division_port = port

        super(LatticeCompartment, self).__init__(processes, derivers, states, configuration)

    def run_incremental(self, run_until):
        while self.time() < run_until:
            self.update(self.time_step)

    def apply_outer_update(self, update):
        self.last_update = update
        env_keys = update['concentrations'].keys()
        environment = self.states.get(self.environment_port)

        if len(self.exchange_ids) > 0:
            # update only the states defined in both exchange and the external environment
            exchange = self.states.get(self.exchange_port)
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
        states = self.divide_state()

        return [
            dict(
                id=str(uuid.uuid1()),
                volume=daughter_state[self.volume_port]['volume'],
                boot_config=dict(
                    initial_time=self.time(),
                    initial_state=daughter_state,
                    volume=daughter_state[self.volume_port]['volume']))
            for daughter_state in states]

    def generate_inner_update(self):
        values = {}
        volume = 1.0
        motile_force = [0.0, 0.0]

        exchange = self.states.get(self.exchange_port)
        if exchange:
            environment_change = exchange.state_for(self.exchange_ids)
        else:
            environment_change = {}

        if self.volume_port:
            volume_state = self.states[self.volume_port].state_for(['volume'])
            volume = volume_state['volume']

        if self.motile_port:
            motile_state = self.states.get(self.motile_port)
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
            'environment_change': environment_change})

        return values


def generate_lattice_compartment(process, config):

    # get the process' default settings
    default_settings = process.default_settings()
    process_id = default_settings.get('process_id', 'process')

    # declare the processes layers (with a single layer)
    processes_layers = [{process_id: process}]

    # make a simple topology mapping 'port' to 'port'
    process_ports = process.ports.keys()
    topology = {process_id: {port: port for port in process_ports}}

    # add derivers
    derivers = get_derivers(processes_layers, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes_layers + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])

    # initialize the states for each port
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

    # get the time step
    time_step = get_minimum_timestep(processes_layers)

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
        'exchange_port': 'exchange',  # TODO -- get this state id from a default_config() function in the process
        'environment_port': 'external',  # TODO -- get this state id from a default_config() function in the process
        'topology': topology,
    }

    # create the lattice compartment
    return LatticeCompartment(processes_layers, deriver_processes, states, options)



## functions for testing
def simulate_lattice_compartment(compartment, settings={}):
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

        values = compartment.generate_inner_update()
        if 'division' in values:
            break

    return saved_state


def divide_composite(config):

    def divide_condition(compartment):
        division_port = compartment.division_port
        division = compartment.states[division_port].state_for(['division'])
        if division.get('division', 0) == 0:  # 0 means false
            divide = False
        else:
            divide = True
        return divide

    # toy processes
    class ToyGrowth(Process):
        def __init__(self, initial_parameters={}):
            ports = {'pool': ['MASS', 'volume']}
            parameters = {
                'growth_rate': 0.1}
            parameters.update(initial_parameters)

            super(ToyGrowth, self).__init__(ports, parameters)

        def next_update(self, timestep, states):
            mass = states['pool']['MASS']
            volume = states['pool']['volume']
            new_mass = mass * self.parameters['growth_rate'] * timestep
            new_volume = volume * self.parameters['growth_rate'] * timestep
            return {
                'pool':
                    {'MASS': new_mass,
                     'volume': new_volume}}

    class ToyDivide(Process):
        def __init__(self, initial_parameters={}):
            self.division = 0
            ports = {'pool': ['MASS', 'division']}
            parameters = {
                'division_mass': 20}
            parameters.update(initial_parameters)

            super(ToyDivide, self).__init__(ports, parameters)

        def next_update(self, timestep, states):
            mass = states['pool']['MASS']

            if mass >= self.parameters['division_mass']:
                self.division = 1

            return {
                'pool': {
                    'division': self.division}}

    # declare processes in list
    processes = [
        {'growth': ToyGrowth()},
        {'divide': ToyDivide()}]

    # declare derivers
    derivers = []

    # declare the states
    states = {
        'cell': Store(
            initial_state={'MASS': 10, 'volume': 1, 'division': 0},
            schema={
                'division': {
                    'updater': 'set',
                    'divide': 'zero'}})}

    # hook up the ports in each process to compartment states
    topology = {
        'growth': {
            'pool': 'cell'},
        'divide': {
            'pool': 'cell'}}

    # emitter that prints to the terminal
    emitter = get_emitter({
        'type': 'print',
        'keys': {
            'cell': ['MASS']}})

    options = {
        'emitter': emitter,
        'topology': topology,
        'initial_time': 0.0,
        'divide_condition': divide_condition}

    return {
        'processes': processes,
        'derivers': derivers,
        'states': states,
        'options': options}

def test_divide(composite=divide_composite):
    # set up the the composite
    composite_config = composite({})
    processes = composite_config['processes']
    derivers = composite_config['derivers']
    states = composite_config['states']
    options = composite_config['options']
    # topology = options['topology']

    lattice_compartment = LatticeCompartment(processes, derivers, states, options)

    settings = {
        'timestep': 1,
        'total_time': 20}

    saved_state = simulate_lattice_compartment(lattice_compartment, settings)

    # assert division
    times = list(saved_state.keys())
    assert saved_state[times[-1]]['cell']['division'] == 1
    return saved_state

if __name__ == '__main__':
    saved_state = test_divide()
    for time, state in saved_state.items():
        print('{}: {}'.format(time,state))
