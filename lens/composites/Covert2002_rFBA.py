from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states, merge_default_updaters, deep_merge
from lens.utils.dict_utils import merge_dicts

# processes
from lens.processes.CovertPalsson2002_metabolism import Covert2002Metabolism
from lens.processes.CovertPalsson2002_regulation import Regulation
from lens.processes.derive_volume import DeriveVolume


def compose_covert2002(config):
    exchange_key = config.get('exchange_key')

    # declare the processes
    metabolism = Covert2002Metabolism(config)
    regulation = Regulation(config)
    deriver = DeriveVolume(config)
    processes = [
        {'regulation': regulation,
         'metabolism': metabolism},
        {'deriver': deriver}]

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

    for role, state_ids in default_updaters.iteritems():
        for state_id, updater in state_ids.iteritems():
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
                default_states[role],
                dict(initial_state.get(compartment_roles[role], {}))),
            updaters=default_updaters.get(role, {}))
            for role in default_states.keys()}

    # configure the states to the roles for each process
    topology = {
        'regulation': {
            'external': 'environment',
            'internal': 'cell'},
        'metabolism': {
            'external': 'environment',
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    options = {
        'topology': topology,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
