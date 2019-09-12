from __future__ import absolute_import, division, print_function

from lens.actor.inner import Inner
from lens.actor.boot import BootAgent
from lens.actor.process import State
from lens.actor.emitter import get_emitter
from lens.environment.lattice_compartment import LatticeCompartment, generate_lattice_compartment

# processes
from lens.processes.transport_lookup import TransportLookup
from lens.processes.CovertPalsson2002_metabolism import Metabolism
from lens.processes.Kremling2007 import Transport

## example composition for one process:
# def initialize_lookup_transport(config):
#     # declare the processes
#     transport = TransportLookup()
#     processes = {
#         'transport': transport}
#
#     # initialize states
#     defaults = transport.default_state()
#     states = {
#         'environment': State(defaults['external']),
#         'cell': State(defaults['internal'])}
#
#     # configure the states to the roles for each process
#     topology = {
#         'transport': {
#             'external': 'environment',
#             'internal': 'cell'}}
#
#     # configure emitter
#     emitter_config = config.get('emitter', {})
#     emitter_config['keys'] = {'environment': [], 'cell': []}
#     emitter = get_emitter(emitter_config)
#
#     options = {
#         'topology': topology,
#         'emitter': emitter,
#         'environment': 'environment',
#         'compartment': 'cell',
#         'external_molecules': defaults['external_molecules']}
#
#     # create the compartment
#     return LatticeCompartment(processes, states, options)


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
                }
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
            'transport': wrap_boot(wrap_initialize(Transport), {'volume': 1.0})
        }

def run():
    boot = BootCompartment()
    boot.execute()

if __name__ == '__main__':
    run()
