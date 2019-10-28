from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states, merge_default_updaters, dict_merge

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.growth import Growth
from lens.processes.division import Division


def compose_growth_division(config):
    exchange_key = config.get('exchange_key')

    # declare the processes
    growth = Growth(config)
    division = Division(config)
    deriver = DeriveVolume(config)
    processes = {
        'growth': growth,
        'division': division,
        'deriver': deriver,
    }

    # configure the states to the roles for each process
    topology = {
        'growth': {
            'internal': 'cell'},
        'division': {
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
    for process_id, process in processes.iteritems():
        roles = {role: {} for role in process.roles.keys()}
        initial_exchanges.update(roles)

    # set states according to the compartment_roles mapping.
    # This will not generalize to composites with processes that have different roles
    compartment_roles = {
        'external': 'environment',
        'internal': 'cell'}

    states = {
        compartment_roles[role]: State(
            initial_state=dict_merge(
                default_states.get(role, {}),
                dict(initial_state.get(compartment_roles[role], {}))),
            updaters=default_updaters.get(role, {}))
        for role in default_states.keys()}

    options = {
        'topology': topology,
        'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids}

    return {
        'processes': processes,
        'states': states,
        'options': options}
