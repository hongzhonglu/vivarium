from __future__ import absolute_import, division, print_function

from lens.actor.process import merge_default_updaters, deep_merge, initialize_state
from lens.utils.dict_utils import merge_dicts

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity
from lens.processes.membrane_potential import MembranePotential
from lens.processes.Kremling2007_transport import Transport
from lens.processes.derive_volume import DeriveVolume
from lens.processes.division import Division, divide_condition, divide_state  # TODO -- division process can house all condition, state functions


def compose_pmf_chemotaxis(config):
    exchange_key = config.get('exchange_key')
    receptor_parameters = {'ligand': 'GLC'}
    receptor_parameters.update(config)

    # declare the processes
    receptor = ReceptorCluster(receptor_parameters)
    motor = MotorActivity(config)
    PMF = MembranePotential(config)
    transport = Transport(config)
    deriver = DeriveVolume(config)
    division = Division(config)

    # place processes in layers
    processes = [
        {'PMF': PMF},
        {'receptor': receptor, 'transport': transport},
        {'motor': motor},
        {'deriver': deriver,
         # 'division': division
         },
    ]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cell'},
        'transport': {
            'external': 'environment',
            'internal': 'cell'},
        'motor': {
            'external': 'environment',
            'internal': 'cell'},
        'PMF': {
            'external': 'environment',
            'membrane': 'membrane',
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        # 'division': {
        #     'internal': 'cell'},
        }


    # TODO -- remove this
    default_updaters = merge_default_updaters(processes)
    environment_ids = []
    initial_exchanges = {'cell': {}, 'environment': {}} #role: {} for role in process.roles.keys()}
    for roleP, state_ids in default_updaters.iteritems():
        if roleP is 'external':
            role = 'environment'
        elif roleP is 'internal':
            role = 'cell'
        for state_id, updater in state_ids.iteritems():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})



    # initialize the states
    states = initialize_state(processes, topology, initial_exchanges)  # TODO -- dont use initial_exchanges
    # states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
        # 'divide_condition': divide_condition,
        # 'divide_state': divide_state
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}


def test_PMF_chemotaxis():
    import numpy as np
    from lens.actor.process import Compartment
    from lens.environment.lattice_compartment import LatticeCompartment

    exchange_key = '__exchange'

    boot_config = {'exchange_key': exchange_key}
    composite_config = compose_pmf_chemotaxis(boot_config)
    processes = composite_config['processes']
    states = composite_config['states']
    options = composite_config['options']
    options.update({'exchange_key': exchange_key})

    # make compartment
    compartment = LatticeCompartment(processes, states, options)

    print(compartment.current_parameters())
    print(compartment.current_state())

    # # test compartment
    # compartment = Compartment(processes, states, options)
    #
    # print(compartment.current_parameters())
    # print(compartment.current_state())
    #
    # # evaluate compartment
    # timestep = 1
    # for steps in np.arange(13):
    #     compartment.update(timestep)
    #     print(compartment.current_state())


    # make lattice_compartment
    lattice_compartment = LatticeCompartment(processes, states, options)

    print(lattice_compartment.current_parameters())
    print(lattice_compartment.current_state())

    # evaluate compartment
    timestep = 1
    for steps in np.arange(13):
        lattice_compartment.update(timestep)
        print(lattice_compartment.current_state())


if __name__ == '__main__':
    saved_state = test_PMF_chemotaxis()
