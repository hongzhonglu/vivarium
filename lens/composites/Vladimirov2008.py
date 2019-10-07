from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_initial_states
from lens.actor.emitter import get_emitter, configure_emitter
from lens.environment.lattice_compartment import LatticeCompartment

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_flagella import FlagellaActivity
from lens.processes.derive_volume import DeriveVolume


def initialize_vladimirov2008(config):
    config.update({
        'emitter': {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
        }
    })

    # declare the processes
    receptor = ReceptorCluster(config)
    flagella = FlagellaActivity(config)
    deriver = DeriveVolume(config)
    processes = {
        'receptor': receptor,
        'flagella': flagella,
        'deriver': deriver}

    # initialize the states
    initial_state = merge_initial_states(processes)
    states = {
        'environment': State(initial_state['external']),
        'cell': State(initial_state['internal'])}

    # configure the states to the roles for each process
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cell'},
        'flagella': {
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