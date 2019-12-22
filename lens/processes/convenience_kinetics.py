from __future__ import absolute_import, division, print_function

from lens.actor.process import Process
from lens.utils.kinetic_rate_laws import KineticFluxModel
from lens.utils.dict_utils import merge_dicts


class ConvenienceKinetics(Process):

    def __init__(self, initial_parameters={}):
        self.reactions = initial_parameters.get('reactions')
        kinetic_parameters = initial_parameters.get('kinetic_parameters')
        roles = initial_parameters.get('roles')
        self.initial_state = initial_parameters.get('initial_state')

        # Make the kinetic model
        self.kinetic_rate_laws = KineticFluxModel(self.reactions, kinetic_parameters)

        roles.update({'fluxes': self.kinetic_rate_laws.reaction_ids})

        parameters = {}
        super(ConvenienceKinetics, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state

        # default emitter keys
        default_emitter_keys = {}

        # default updaters
        default_updaters = {}

        default_settings = {
            'process_id': 'template',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'time_step': 1.0}

        return default_settings

    def next_update(self, timestep, states):

        # kinetic rate law requires a flat dict with states.
        # TODO -- might need to distinguish different roles with suffixes (_external, _internal, etc)
        flattened_states = merge_dicts([role_states for role, role_states in states.items()])

        # get flux
        fluxes = self.kinetic_rate_laws.get_fluxes(flattened_states)

        # apply fluxes to state
        update = {role: {} for role in self.roles.keys()}
        update.update({'fluxes': fluxes})
        for reaction_id, flux in fluxes.items():
            stoichiometry = self.reactions[reaction_id]['stoichiometry']
            for state_id, coeff in stoichiometry.items():
                # which role is state_id in? TODO -- there must be a better way to do this.
                for role_id, state_list in self.roles.items():
                    if state_id in state_list:
                        update[role_id][state_id] = coeff * flux

        return update



# testing functions
toy_reactions = {
    'reaction1': {
        'stoichiometry': {
            'A': 1,
            'B': -1},
        'is reversible': False,
        'catalyzed by': ['enzyme1']}
    }

toy_kinetics = {
    'reaction1': {
        'enzyme1': {
            'B': 0.1,
            'kcat_f': 0.1}
        }
    }

toy_roles = {
    'internal': ['A', 'enzyme1'],
    'external': ['B'],
    }

toy_initial_state = {
    'internal': {
        'A': 1.0,
        'enzyme1': 1.0},
    'external': {
        'B': 1.0},
    'fluxes': {
        'reaction1': 0.0}
    }

# test
def test_convenience_kinetics():
    toy_config = {
        'reactions': toy_reactions,
        'kinetic_parameters': toy_kinetics,
        'initial_state': toy_initial_state,
        'roles': toy_roles}

    kinetic_process = ConvenienceKinetics(toy_config)

    # get initial state and parameters
    settings = kinetic_process.default_settings()
    state = settings['state']

    # run the simulation
    timestep = 1
    for step in range(10):
        # get update
        update = kinetic_process.next_update(timestep, state)

        # apply update
        for role_id, states_update in update.items():
            for state_id, change in states_update.items():
                state[role_id][state_id] += change
        print(state)


if __name__ == '__main__':
    test_convenience_kinetics()
