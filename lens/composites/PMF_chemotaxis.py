from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states, merge_default_updaters, dict_merge

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity
from lens.processes.membrane_potential import MembranePotential


def compose_pmf_chemotaxis(config):
    exchange_key = config.get('exchange_key')

    # declare the processes
    receptor = ReceptorCluster(config)
    motor = MotorActivity(config)
    PMF = MembranePotential(config)
    # TODO -- add flagella process.  for now assume fixed flagella.

    processes = {
        'receptor': receptor,
        'motor': motor,
        'PMF': PMF,
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

    for role, state_ids in default_updaters.iteritems():
        for state_id, updater in state_ids.iteritems():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})

    default_states = dict_merge(default_states, initial_exchanges)

    # set states according to the compartment_roles mapping.
    # This will not generalize to composites with processes that have different roles
    compartment_roles = {
        'external': 'environment',
        'internal': 'cytoplasm'}

    states = {
        compartment_roles[role]: State(
            initial_state=dict_merge(
                default_states.get(role, {}),
                dict(initial_state.get(compartment_roles[role], {}))),
            updaters=default_updaters.get(role, {}))
        for role in default_states.keys()}

    # configure the states to the roles for each process
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cytoplasm'},
        'motor': {
            'external': 'environment',
            'internal': 'cytoplasm'},
        'PMF': {
            'external': 'environment',
            'membrane': 'membrane',
            'internal': 'cytoplasm'},
        }

    options = {
        'topology': topology,
		'initial_time': initial_time,
        'environment': 'environment',
        'compartment': 'cytoplasm',   # TODO -- should this include both cytoplasm and membrane?
        'environment_ids': environment_ids,
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}
