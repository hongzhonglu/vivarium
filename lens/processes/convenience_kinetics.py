from __future__ import absolute_import, division, print_function

from lens.actor.process import Process
from lens.utils.kinetic_rate_laws import KineticFluxModel
from lens.utils.units import units
from lens.utils.dict_utils import merge_dicts


class ConvenienceKinetics(Process):

    def __init__(self, initial_parameters={}):
        reactions = initial_parameters.get('reactions')
        kinetic_parameters = initial_parameters.get('kinetic_parameters')
        roles = initial_parameters.get('roles')
        self.initial_state = initial_parameters.get('initial_state')

        # Make the kinetic model
        self.kinetic_rate_laws = KineticFluxModel(reactions, kinetic_parameters)

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
        # TODO -- might need to distinguish different roles with suffixes,
        flattened_states = merge_dicts([role_states for role, role_states in states.items()])

        fluxes = self.kinetic_rate_laws.get_fluxes(flattened_states)

        return {'fluxes': fluxes}



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
            'kcat_f': 1.0}
        }
    }

toy_roles = {
    'internal': [
        'A',
        'enzyme1'],
    'external': [
        'B']
    }

toy_initial_state = {
    'internal': {
        'A': 1.0,
        'enzyme1': 1.0},
    'external': {
        'B': 1.0}
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

    timestep = 1
    update = kinetic_process.next_update(timestep, state)


if __name__ == '__main__':
    test_convenience_kinetics()