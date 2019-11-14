from __future__ import absolute_import, division, print_function

import random

from lens.actor.process import State, merge_default_states, merge_default_updaters, deep_merge
from lens.utils.dict_utils import merge_dicts

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.growth import Growth
from lens.processes.division import Division
from lens.processes.protein_expression import ProteinExpression


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

def compose_growth_division(config):
    exchange_key = config.get('exchange_key')

    # declare the processes
    growth = Growth(config)
    division = Division(config)
    expression = ProteinExpression(config)
    deriver = DeriveVolume(config)
    processes = [
        {'growth': growth,
         'division': division,
         'expression': expression},
        {'deriver': deriver}]

    # configure the states to the roles for each process
    topology = {
        'growth': {
            'internal': 'cell'},
        'division': {
            'internal': 'cell'},
        'expression': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # initialize the states
    default_states = merge_default_states(processes)
    default_updaters = merge_default_updaters(processes)
    initial_state = config.get('initial_state', {})
    initial_time = config.get('initial_time', 0.0)

    # get environment ids, and make exchange_ids for external state
    environment_ids = []
    initial_exchanges = {}
    for process_id, process in merge_dicts(processes).iteritems():
        roles = {role: {} for role in process.roles.keys()}
        initial_exchanges.update(roles)

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

    options = {
        'topology': topology,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}
