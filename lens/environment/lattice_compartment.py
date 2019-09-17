from __future__ import absolute_import, division, print_function

from lens.actor.process import Compartment, State
from lens.actor.emitter import get_emitter
from lens.actor.inner import Simulation

DEFAULT_COLOR = [color/255 for color in [153, 204, 255]]

class LatticeCompartment(Compartment, Simulation):
    def __init__(self, processes, states, configuration):
        self.environment = configuration['environment']
        self.compartment = configuration['compartment']
        self.exchange_key = configuration['exchange_key']
        self.environment_deltas = configuration['environment_deltas']
        self.configuration = configuration
        self.color = DEFAULT_COLOR

        super(LatticeCompartment, self).__init__(processes, states, configuration)

    def run_incremental(self, run_until):
        while self.time() < run_until:
            self.update(self.time_step)

    def apply_outer_update(self, update):
        self.last_update = update
        environment = self.states[self.environment]
        environment.assign_values(update['concentrations'])
        environment.assign_values({key: 0 for key in self.environment_deltas})  # reset delta counts to 0

    def generate_inner_update(self):
        environment = self.states[self.environment]
        changes = environment.state_for(self.environment_deltas)
        compartment = self.states[self.compartment]
        environment_change = {mol_id.replace(self.exchange_key, ''): value
            for mol_id, value in changes.iteritems()}

        values = compartment.state_for(['volume'])
        values.update({
            'motile_force': [0,0], # TODO -- get motile_force from compartment state
            'color': self.color,
            'environment_change': environment_change,
            'transport_fluxes': {},  # TODO -- remove this
        })

        return values


def generate_lattice_compartment(process, config):
    # declare the processes
    processes = {
        'process': process}

    # initialize states
    defaults = process.default_state()
    states = {
        role: State(defaults.get(role, {}))
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
        'environment': config.get('environment', 'external'),
        'compartment': config.get('compartment', 'internal'),
        'exchange_key': config['exchange_key'],
        'environment_deltas': defaults['environment_deltas']}

    options.update(config['compartment_options'])

    # create the compartment
    return LatticeCompartment(processes, states, options)
