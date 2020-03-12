from __future__ import absolute_import, division, print_function

import os

from scipy import constants

from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process_with_environment,
    plot_simulation_output,
    flatten_timeseries,
    save_timeseries,
    load_timeseries,
    REFERENCE_DATA_DIR,
    TEST_OUT_DIR,
    assert_timeseries_close,
)
from vivarium.utils.kinetic_rate_laws import KineticFluxModel
from vivarium.utils.dict_utils import tuplify_port_dicts
from vivarium.utils.units import units

EMPTY_ROLES = {
    'internal': [],
    'external': []}

EMPTY_STATES = {
    'internal': {},
    'external': {}}
NAME = 'convenience_kinetics'


class ConvenienceKinetics(Process):

    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1 / units.mol

        # retrieve initial parameters
        self.reactions = initial_parameters.get('reactions', {})
        self.initial_state = initial_parameters.get('initial_state', EMPTY_STATES)
        kinetic_parameters = initial_parameters.get('kinetic_parameters', {})
        ports = initial_parameters.get('ports', EMPTY_ROLES)

        # make the kinetic model
        self.kinetic_rate_laws = KineticFluxModel(self.reactions, kinetic_parameters)

        # ports
        # fluxes port is used to pass constraints
        # exchange is equivalent to external, for lattice_compartment
        ports.update({
            'fluxes': self.kinetic_rate_laws.reaction_ids,
            'exchange': ports['external'],
            'global': ['mmol_to_counts']})

        # parameters
        parameters = {}
        parameters.update(initial_parameters)

        super(ConvenienceKinetics, self).__init__(ports, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state

        # default emitter keys
        emit_ports = ['internal', 'external']
        default_emitter_keys = {
            port: state_list
            for port, state_list in self.ports.items() if port in emit_ports}

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

        # kinetic rate law requires a flat dict with ('port', 'state') keys.
        flattened_states = tuplify_port_dicts(states)

        # get flux
        fluxes = self.kinetic_rate_laws.get_fluxes(flattened_states)

        # make the update
        # add fluxes to update
        update = {port: {} for port in self.ports.keys()}
        update.update({'fluxes': fluxes})

        # get exchange
        for reaction_id, flux in fluxes.items():
            stoichiometry = self.reactions[reaction_id]['stoichiometry']
            for port_state_id, coeff in stoichiometry.items():
                for port_id, state_list in self.ports.items():
                    # separate the state_id and port_id
                    if port_id in port_state_id:
                        state_id = port_state_id[1]
                        state_flux = coeff * flux * timestep

                        if port_id == 'external':
                            # convert exchange fluxes to counts with mmol_to_counts
                            # TODO -- use deriver to get exchanges
                            delta_counts = int((state_flux * mmol_to_counts).magnitude)
                            update['exchange'][state_id] = (
                                update['exchange'].get(state_id, 0)
                                + delta_counts
                            )
                        else:
                            update[port_id][state_id] = (
                                update[port_id].get(state_id, 0)
                                + state_flux
                            )

        # note: external and internal ports update change in mmol.
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

    transport_ports = {
        'internal': [
            'g6p_c', 'pep_c', 'pyr_c', 'h_c', 'PTSG', 'LacY'],
        'external': [
            'glc__D_e', 'lac__D_e', 'h_e']}

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'ports': transport_ports}

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

    toy_ports = {
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
        'ports': toy_ports}


def test_convenience_kinetics(end_time=1000):
    config = get_glc_lct_config()
    kinetic_process = ConvenienceKinetics(config)

    settings = {
        'environment_port': 'external',
        'exchange_port': 'exchange',
        'emit_timeseries': True,
        'environment_volume': 1e-13,  # L
        'timestep': 1,
        'total_time': end_time}

    saved_state = simulate_process_with_environment(kinetic_process, settings)
    return saved_state


def test_convenience_kinetics_correlated_to_reference():
    timeseries = test_convenience_kinetics()
    flattened = flatten_timeseries(timeseries)
    reference_timeseries = load_timeseries(
        os.path.join(REFERENCE_DATA_DIR, NAME + '.csv'))
    assert_timeseries_close(flattened, reference_timeseries)


if __name__ == '__main__':
    out_dir = os.path.join(TEST_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}

    timeseries = test_convenience_kinetics()
    plot_simulation_output(timeseries, plot_settings, out_dir)
    save_timeseries(timeseries, out_dir)
