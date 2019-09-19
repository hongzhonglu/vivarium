from __future__ import absolute_import, division, print_function

import collections

from lens.actor.inner import Inner
from lens.actor.boot import BootAgent
from lens.actor.process import State
from lens.actor.emitter import get_emitter, configure_emitter
from lens.environment.lattice_compartment import LatticeCompartment, generate_lattice_compartment

# processes
from lens.processes.transport_lookup import TransportLookup
from lens.processes.CovertPalsson2002_metabolism import Metabolism
from lens.processes.Kremling2007 import Transport
from lens.processes.derive_volume import DeriveVolume
from lens.processes.division import Division
from lens.processes.growth import Growth


def merge_initial_states(processes):
    initial_state = {}
    for process_id, process in processes.iteritems():
        default = process.default_state()
        dict_merge(initial_state, default)
    return initial_state

def dict_merge(dct, merge_dct):
    ''' Recursive dict merge '''
    for k, v in merge_dct.iteritems():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

# composites
def initialize_covert2008(config):
    config.update({
        'exchange_key': '__exchange',  # key for counting exchange with lattice
        'emitter': {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
        }
    })

    # declare the processes
    transport = Transport(config)
    deriver = DeriveVolume(config)
    processes = {
        'transport': transport,
        'deriver': deriver}

    # initialize the states
    initial_state = merge_initial_states(processes)
    states = {
        'environment': State(initial_state['external']),
        'cell': State(initial_state['internal'])}

    # configure the states to the roles for each process
    topology = {
        'transport': {
            'external': 'environment',
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # configure emitter
    emitter = configure_emitter(config, processes, topology)

    options = {
        'topology': topology,
        'emitter': emitter,
        'environment': 'environment',
        'compartment': 'cell',
        'exchange_key': config['exchange_key'],
        'environment_ids': initial_state['environment_ids'],
        'environment_deltas': initial_state['environment_deltas']}

    # create the compartment
    return LatticeCompartment(processes, states, options)

def initialize_growth_division(config):
    config.update({
        'exchange_key': '__exchange',  # key for counting exchange with lattice
        'emitter': {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
        }
    })

    # declare the processes
    growth = Growth(config)
    deriver = DeriveVolume(config)
    # division = Division(config)
    processes = {
        'growth': growth,
        'deriver': deriver,
        # 'division': division
    }

    # configure the states to the roles for each process
    topology = {
        'growth': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        # 'division': {
        #     'internal': 'cell'},
        }

    # initialize the states
    # TODO (Eran) -- this can be refactored
    initial_state = merge_initial_states(processes)
    states = {
        # 'environment': State(initial_state['external']),
        'cell': State(initial_state['internal'])}

    # configure emitter
    emitter = configure_emitter(config, processes, topology)

    options = {
        'topology': topology,
        'emitter': emitter,
        'environment': 'environment',
        'compartment': 'cell',
        'exchange_key': config['exchange_key'],
        'environment_ids': initial_state['environment_ids'],
        'environment_deltas': initial_state['environment_deltas']}

    # create the compartment
    return LatticeCompartment(processes, states, options)

def wrap_boot(initialize, initial_state):
    def boot(agent_id, agent_type, agent_config):
        initial_state.update(agent_config.get('declare', {}))
        agent_config['declare'] = initial_state  # 'declare' is for the environment

        return Inner(
            agent_id,
            agent_type,
            agent_config,
            initialize)

    return boot

def wrap_initialize(make_process):
    def initialize(boot_config):
        boot_config.update({
            'exchange_key': '__exchange',
            'emitter': {
                'type': 'database',
                'url': 'localhost:27017',
                'database': 'simulations',
                },
            'compartment_options':{
                'time_step': 10.0,
                },
            })
        process = make_process(boot_config)  # 'boot_config', set in environment.control is the process' initial_parameters
        return generate_lattice_compartment(process, boot_config)

    return initialize


class BootCompartment(BootAgent):
    def __init__(self):
        super(BootCompartment, self).__init__()
        self.agent_types = {
            'lookup': wrap_boot(wrap_initialize(TransportLookup), {'volume': 1.0}),
            'metabolism': wrap_boot(wrap_initialize(Metabolism), {'volume': 1.0}),
            'transport': wrap_boot(wrap_initialize(Transport), {'volume': 1.0}),
            'growth': wrap_boot(wrap_initialize(Growth), {'volume': 1.0}),
            'growth_division': wrap_boot(initialize_growth_division, {'volume': 1.0}),
            'covert2008': wrap_boot(initialize_covert2008, {'volume': 1.0})
        }

def run():
    boot = BootCompartment()
    boot.execute()

if __name__ == '__main__':
    run()
