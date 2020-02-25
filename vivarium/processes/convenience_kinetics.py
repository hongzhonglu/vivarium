from __future__ import absolute_import, division, print_function

import os

from scipy import constants

from vivarium.actor.process import Process
from vivarium.actor.composition import simulate_process_with_environment, convert_to_timeseries, plot_simulation_output
from vivarium.utils.kinetic_rate_laws import KineticFluxModel
from vivarium.utils.dict_utils import tuplify_role_dicts
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

        # retrieve initial parameters
        self.reactions = initial_parameters.get('reactions', {})
        self.initial_state = initial_parameters.get('initial_state', EMPTY_STATES)
        kinetic_parameters = initial_parameters.get('kinetic_parameters', {})
        roles = initial_parameters.get('roles', EMPTY_ROLES)

        # make the kinetic model
        self.kinetic_rate_laws = KineticFluxModel(self.reactions, kinetic_parameters)

        # roles
        # fluxes role is used to pass constraints
        # exchange is equivalent to external, for lattice_compartment
        roles.update({
            'fluxes': self.kinetic_rate_laws.reaction_ids,
            'exchange': roles['external'],
            'global': ['mmol_to_counts']})

        # parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(ConvenienceKinetics, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'fluxes': {
                flux_id : {
                    'updater': 'set'}
                for flux_id in self.kinetic_rate_laws.reaction_ids}}

        default_settings = {
            'process_id': 'convenience_kinetics',
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'time_step': 1.0}

        return default_settings

    def next_update(self, timestep, states):

        # get mmol_to_counts for converting flux to exchange counts
        mmol_to_counts = states['global']['mmol_to_counts'] * units.L / units.mmol

        # kinetic rate law requires a flat dict with ('role', 'state') keys.
        flattened_states = tuplify_role_dicts(states)

        # get flux
        fluxes = self.kinetic_rate_laws.get_fluxes(flattened_states)

        # make the update
        # add fluxes to update
        update = {role: {} for role in self.roles.keys()}
        update.update({'fluxes': fluxes})

        # get exchange
        for reaction_id, flux in fluxes.items():
            stoichiometry = self.reactions[reaction_id]['stoichiometry']
            for role_state_id, coeff in stoichiometry.items():
                for role_id, state_list in self.roles.items():
                    # separate the state_id and role_id
                    if role_id in role_state_id:
                        state_id = role_state_id[1]
                        state_flux = (coeff * flux * timestep *
                            units.mmol / units.L)
                        if role_id == 'external':
                            # convert exchange fluxes to counts with mmol_to_counts
                            # TODO -- use deriver to get exchanges
                            delta_counts = int((state_flux * mmol_to_counts).magnitude)
                            update['exchange'][state_id] = delta_counts
                        else:
                            update[role_id][state_id] = (
                                update[role_id].get(state_id, 0)
                                + state_flux.magnitude
                            )

        # note: external and internal roles update change in mmol.
        return update



# functions
def get_glc_lct_config():
    """
    Convenience kinetics configuration for simplified glucose and lactose transport.
    Glucose updake simplifies the PTS/GalP system to a single uptake kinetic
    with glc__D_e_external as the only cofactor.
    """
    transport_reactions = {
        'EX_glc__D_e': {
            'stoichiometry': {
                ('internal', 'g6p_c'): 1.0,
                ('external', 'glc__D_e'): -1.0,
                ('internal', 'pep_c'): -1.0,  # TODO -- PEP requires homeostasis mechanism to avoid depletion
                ('internal', 'pyr_c'): 1.0},
            'is reversible': False,
            'catalyzed by': [('internal', 'PTSG')]},
        'EX_lac__D_e': {
            'stoichiometry': {
                ('external', 'lac__D_e'): -1.0,
                ('external', 'h_e'): -1.0,
                ('internal', 'lac__D_c'): 1.0,
                ('internal', 'h_c'): 1.0},
            'is reversible': False,
            'catalyzed by': [('internal', 'LacY')]}}

    transport_kinetics = {
        'EX_glc__D_e': {
            ('internal', 'PTSG'): {
                # k_m for external [glc__D_e]
                ('external', 'glc__D_e'): 1e-1,
                # Set k_m = None to make a reactant non-limiting
                ('internal', 'pep_c'): None,
                'kcat_f': 3e5}},  # kcat for the forward direction
        'EX_lac__D_e': {
            ('internal', 'LacY'): {
                ('external', 'lac__D_e'): 1e-1,
                ('external', 'h_e'): None,
                'kcat_f': 5e4}}}

    transport_initial_state = {
        'internal': {
            'PTSG': 1.8e-6,  # concentration (mmol/L)
            'g6p_c': 0.0,
            'pep_c': 1.8e-1,
            'pyr_c': 0.0,
            'LacY': 0.0,
            'lac__D_c': 0.0,
            'h_c': 100.0},
        'external': {
            'glc__D_e': 12.0,
            'lac__D_e': 10.0,
            'h_e': 100.0},
        'fluxes': {  # TODO -- is this needed?
            'EX_glc__D_e': 0.0,
            'EX_lac__D_e': 0.0}}

    transport_roles = {
        'internal': [
            'g6p_c', 'pep_c', 'pyr_c', 'h_c', 'PTSG', 'LacY'],
        'external': [
            'glc__D_e', 'lac__D_e', 'h_e']}

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}

def get_toy_config():
    toy_reactions = {
        'reaction1': {
            'stoichiometry': {
                ('internal', 'A'): 1,
                ('external', 'B'): -1},
            'is reversible': False,
            'catalyzed by': [('internal', 'enzyme1')]}}

    toy_kinetics = {
        'reaction1': {
            ('internal', 'enzyme1'): {
                ('external', 'B'): 0.2,
                'kcat_f': 5e1}}}

    toy_roles = {
        'internal': ['A', 'enzyme1'],
        'external': ['B']}

    toy_initial_state = {
        'internal': {
            'A': 1.0,
            'enzyme1': 1e-1},
        'external': {
            'B': 10.0},
        'fluxes': {
            'reaction1': 0.0}}

    return {
        'reactions': toy_reactions,
        'kinetic_parameters': toy_kinetics,
        'initial_state': toy_initial_state,
        'roles': toy_roles}

# test
def test_convenience_kinetics(end_time=10):
    # config = get_toy_config()
    config = get_glc_lct_config()
    kinetic_process = ConvenienceKinetics(config)

    settings = {
        'environment_role': 'external',
        'exchange_role': 'exchange',
        'environment_volume': 1e-13,  # L
        'timestep': 1,
        'total_time': end_time}

    saved_state = simulate_process_with_environment(kinetic_process, settings)
    return saved_state


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'convenience_kinetics')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}

    saved_data = test_convenience_kinetics(1000)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
