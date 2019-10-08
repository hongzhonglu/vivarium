from __future__ import absolute_import, division, print_function

from lens.actor.process import Compartment, State, dict_merge
from lens.actor.emitter import get_emitter
from lens.actor.inner import Simulation

DEFAULT_COLOR = [color/255 for color in [153, 204, 255]]

exchange_key = '__exchange'  # TODO -- this is declared in multiple locations

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
            environment.assign_values(update['concentrations'])
            environment.assign_values({key: 0 for key in self.exchange_ids})  # reset exchange

    def generate_inner_update(self):
        environment = self.states.get(self.environment)
        if environment:
            changes = environment.state_for(self.exchange_ids)
            environment_change = {mol_id.replace(self.exchange_key, ''): value for mol_id, value in changes.iteritems()}
        else:
            environment_change = {}

        compartment = self.states[self.compartment]
        values = compartment.state_for(['volume'])
        values.update({
            'motile_force': [0,0], # TODO -- get motile_force from compartment state
            'color': self.color,
            'environment_change': environment_change,
        })

        return values

def generate_lattice_compartment(process, config):
    # declare the processes
    processes = {
        'process': process}

    # initialize states
    default_state = process.default_state()
    default_updaters = process.default_updaters()

    # initialize keys for accumulate_delta updater
    # assumes all environmental updates have been set as `accumulate` updaters
    environment_ids = []
    initial_exchanges = {role: {} for role in process.roles.keys()}
    for role, state_ids in default_updaters.iteritems():
        for state_id, updater in state_ids.iteritems():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})

    states = {
        role: State(
            initial_state=dict_merge(default_state.get(role, {}), initial_exchanges.get(role, {})),
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
        'environment_ids': environment_ids,
        'exchange_key': exchange_key,
        'environment': config.get('environment', 'external'),
        'compartment': config.get('compartment', 'internal')}

    options.update(config['compartment_options'])

    # create the compartment
    return LatticeCompartment(processes, states, options)
