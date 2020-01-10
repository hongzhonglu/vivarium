from __future__ import absolute_import, division, print_function

import os

from scipy import constants

from vivarium.actor.process import Process, convert_to_timeseries, plot_simulation_output
from vivarium.utils.kinetic_rate_laws import KineticFluxModel
from vivarium.utils.dict_utils import flatten_role_dicts
from vivarium.utils.units import units

EMPTY_ROLES = {
    'internal': [],
    'external': []}

EMPTY_STATES = {
    'internal': {},
    'external': {}}


class ConvenienceKinetics(Process):

    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1 / units.mol

        self.reactions = initial_parameters.get('reactions')
        kinetic_parameters = initial_parameters.get('kinetic_parameters')
        roles = initial_parameters.get('roles', EMPTY_ROLES)
        self.initial_state = initial_parameters.get('initial_state', EMPTY_STATES)

        # Make the kinetic model
        self.kinetic_rate_laws = KineticFluxModel(self.reactions, kinetic_parameters)

        # add volume to internal role
        if 'volume' not in roles.get('internal'):
            roles['internal'].append('volume')

        roles.update({
            'fluxes': self.kinetic_rate_laws.reaction_ids,
            'exchange': roles['external']
        })

        parameters = {}
        super(ConvenienceKinetics, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state
        default_state['internal'].update({'volume': 1.2})  # (fL)

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

        # get mmol_to_count for converting flux to exchange counts
        volume = states['internal']['volume'] * units.L
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume

        # kinetic rate law requires a flat dict with 'state_role' keys.
        flattened_states = flatten_role_dicts(states)

        # get flux
        fluxes = self.kinetic_rate_laws.get_fluxes(flattened_states)

        # apply fluxes to state
        update = {role: {} for role in self.roles.keys()}
        update.update({'fluxes': fluxes})

        # get exchange
        for reaction_id, flux in fluxes.items():
            stoichiometry = self.reactions[reaction_id]['stoichiometry']
            for state_role_id, coeff in stoichiometry.items():
                for role_id, state_list in self.roles.items():
                    # separate the state_id and role_id
                    if role_id in state_role_id:
                        role_string = '_{}'.format(role_id)
                        state_id = state_role_id.replace(role_string, '')
                        state_flux = coeff * flux

                        if role_id == 'external':
                            # convert exchange fluxes to counts with mmol_to_count
                            delta_counts = int((state_flux * mmol_to_count).magnitude)
                            update['exchange'][state_id] = delta_counts
                        else:
                            update[role_id][state_id] = state_flux

        # note: external and internal roles get update in change in mmol.
        return update



# testing functions
toy_reactions = {
    'reaction1': {
        'stoichiometry': {
            'A_internal': 1,
            'B_external': -1},
        'is reversible': False,
        'catalyzed by': ['enzyme1_internal']}
    }

toy_kinetics = {
    'reaction1': {
        'enzyme1_internal': {
            'B_external': 0.1,
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
    skip_roles = ['exchange']

    # initialize saved data
    saved_state = {}

    # run the simulation
    time = 0
    timestep = 1
    end_time = 10
    saved_state[time] = state
    while time < end_time:
        time += timestep
        # get update
        update = kinetic_process.next_update(timestep, state)

        # apply update
        for role_id, states_update in update.items():
            if role_id not in skip_roles:
                for state_id, change in states_update.items():
                    state[role_id][state_id] += change
        saved_state[time] = state

    return saved_state


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'convenience_kinetics')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}

    saved_data = test_convenience_kinetics()
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
