from __future__ import absolute_import, division, print_function

import random

from lens.actor.process import State, merge_default_states, merge_default_updaters, deep_merge
from lens.utils.dict_utils import merge_dicts

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity
# from lens.processes.membrane_potential import MembranePotential
from lens.processes.Kremling2007_transport import Transport
from lens.processes.derive_volume import DeriveVolume
from lens.processes.division import Division


def divide_condition(compartment):
    division = compartment.states['cell'].state_for(['division'])
    if division.get('division', 0) == 0:  # 0 is false
        divide = False
    else:
        divide = True
    return divide

def divide_state(compartment):
    divided = [{}, {}]
    for state_key, state in compartment.states.items():
        left = random.randint(0, 1)
        for index in range(2):
            divided[index][state_key] = {}
            for key, value in state.to_dict().items():
                if key == 'division':
                    divided[index][state_key][key] = 0
                else:
                    divided[index][state_key][key] = value // 2 + (value % 2 if index == left else 0)

    print('divided {}'.format(divided))
    return divided


def compose_pmf_chemotaxis(config):
    exchange_key = config.get('exchange_key')
    receptor_parameters = {'ligand': 'GLC'}
    receptor_parameters.update(config)

    # declare the processes
    # TODO -- add flagella process.  for now assume fixed flagella.
    receptor = ReceptorCluster(receptor_parameters)
    motor = MotorActivity(config)
    # PMF = MembranePotential(config)
    transport = Transport(config)
    deriver = DeriveVolume(config)
    division = Division(config)

    processes = [
        # {'PMF': PMF},
        {'receptor': receptor,
         'transport': transport},
        {'motor': motor},
        {'deriver': deriver,
         'division': division},
    ]

    # initialize the states
    default_states = merge_default_states(processes)
    default_updaters = merge_default_updaters(processes)
    initial_state = config.get('initial_state', {})
    initial_time = config.get('initial_time', 0.0)

    # get environment ids, and make exchange_ids for external state
    environment_ids = []
    initial_exchanges = {}
    for process_id, process in merge_dicts(processes).items():
        roles = {role: {} for role in process.roles.keys()}
        initial_exchanges.update(roles)

    for role, state_ids in default_updaters.items():
        for state_id, updater in state_ids.items():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})

    default_states = deep_merge(default_states, initial_exchanges)

    # set states according to the compartment_roles mapping.
    # This will not generalize to composites with processes that have different roles
    compartment_roles = {
        'external': 'environment',
        'internal': 'cell'}

    states = {
        compartment_roles[role]: State(
            initial_state=deep_merge(
                default_states.get(role, {}),
                dict(initial_state.get(compartment_roles[role], {}))),
            updaters=default_updaters.get(role, {}))
        for role in default_states.keys()}

    # configure the states to the roles for each process
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
        'deriver': {
            'internal': 'cell'},
        'division': {
            'internal': 'cell'},
        }

    options = {
        'topology': topology,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
        'divide_condition': divide_condition,
        'divide_state': divide_state
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
