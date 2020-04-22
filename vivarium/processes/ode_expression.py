from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.process import Process
from vivarium.utils.dict_utils import deep_merge, tuplify_port_dicts
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_with_environment,
    plot_simulation_output
)
from vivarium.utils.regulation_logic import build_rule
from vivarium.utils.units import units
from vivarium.processes.derive_globals import AVOGADRO


class ODE_expression(Process):
    '''
    a ode-based mRNA and protein expression process, with boolean logic for regulation
    TODO -- kinetic regulation, cooperativity, autoinhibition, autactivation
    '''

    defaults = {
        'transcription_rates': {},
        'translation_rates': {},
        'degradation_rates': {},
        'protein_map': {},
        'regulation': {},
        'regulators': [],
        'initial_state': {},
    }

    def __init__(self, initial_parameters={}):

        # ode gene expression
        self.transcription = initial_parameters.get(
            'transcription_rates', self.defaults['transcription_rates'])
        self.translation = initial_parameters.get(
            'translation_rates', self.defaults['translation_rates'])
        self.degradation = initial_parameters.get(
            'degradation_rates', self.defaults['degradation_rates'])
        self.protein_map = initial_parameters.get(
            'protein_map', self.defaults['protein_map'])

        # boolean regulation
        regulation_logic = initial_parameters.get(
            'regulation', self.defaults['regulation'])
        self.regulation = {
            gene_id: build_rule(logic) for gene_id, logic in regulation_logic.items()}
        regulators = initial_parameters.get('regulators', self.defaults['regulators'])
        internal_regulators = [state_id for port_id, state_id in regulators if port_id == 'internal']
        external_regulators = [state_id for port_id, state_id in regulators if port_id == 'external']

        # get initial state
        states = list(self.transcription.keys()) + list(self.translation.keys())
        null_states = {'internal': {
            state_id: 0 for state_id in states}}
        initialized_states = initial_parameters.get('initial_state', self.defaults['initial_state'])
        self.initial_state = deep_merge(null_states, initialized_states)
        internal = list(self.initial_state.get('internal', {}).keys())
        external = list(self.initial_state.get('external', {}).keys())

        ports = {
            'internal': internal + internal_regulators,
            'external': external + external_regulators,
            'counts': []}

        parameters = {}
        parameters.update(initial_parameters)

        super(ODE_expression, self).__init__(ports, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state

        # schema
        # don't include if it uses the default
        schema = {
            'internal': {
                state : {
                    'divide': 'set',
                    'units': 'mmol',
                    'updater': 'accumulate'}
                for state in self.ports['internal']}}

        # default emitter keys
        default_emitter_keys = {
            'internal': self.ports['internal'],
            'counts': self.ports['internal']}

        # derivers
        deriver_setting = [{
            'type': 'mmol_to_counts',
            'source_port': 'internal',
            'derived_port': 'counts',
            'keys': self.ports['internal']}]

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'deriver_setting': deriver_setting}

        return default_settings

    def next_update(self, timestep, states):
        internal_state = states['internal']

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
        regulation_state = {}
        for gene_id, reg_logic in self.regulation.items():
            regulation_state[gene_id] = reg_logic(flattened_states)

        internal_update = {}
        # transcription: dM/dt = k_M - d_M * M
        # M: conc of mRNA, k_M: transcription rate, d_M: degradation rate
        # transcription rate abstracts over a number of factors, including gene copy number,
        # RNA polymerase abundance, strength of gene promoters, and availability of nucleotides
        for transcript, rate in self.transcription.items():
            transcript_state = internal_state[transcript]
            # do not transcribe inhibited genes
            if transcript in regulation_state and not regulation_state[transcript]:
                rate = 0

            internal_update[transcript] = \
                (rate - self.degradation.get(transcript, 0) * transcript_state) * timestep

        # translation: dP/dt = k_P * m_P - d_P * P
        # P: conc of protein, m_P: conc of P's transcript, k_P: translation rate, d_P: degradation rate
        # translation rate abstracts over a number of factors, including ribosome availability,
        # strength of mRNA's ribosome binding, availability of tRNAs and free amino acids
        for protein, rate in self.translation.items():
            transcript = self.protein_map[protein]
            transcript_state = internal_state[transcript]
            protein_state = internal_state[protein]

            internal_update[protein] = \
                (rate * transcript_state - self.degradation.get(protein, 0) * protein_state) * timestep

        return {'internal': internal_update}


# test functions
# toy config
def get_lacy_config():
    toy_transcription_rates = {
        'lacy_RNA': 2e-6}

    toy_translation_rates = {
        'LacY': 2e-4}

    toy_protein_map = {
        'LacY': 'lacy_RNA'}

    toy_degradation_rates = {
        'lacy_RNA': 1e-1,
        'LacY': 1e-4}

    initial_state = {
        'internal': {
            'lacy_RNA': 0,
            'LacY': 0.0}}

    return {
        'transcription_rates': toy_transcription_rates,
        'translation_rates': toy_translation_rates,
        'degradation_rates': toy_degradation_rates,
        'protein_map': toy_protein_map,
        'initial_state': initial_state}

def get_flagella_expression():
    transcription = {
        'flag_RNA': 1e-6}

    translation = {
        'flagella': 8e-5}

    degradation = {
        'flag_RNA': 2e-2}  # 1e-23}

    protein_map = {
        'flagella': 'flag_RNA'}

    # get initial concentrations from counts
    volume = 1.2 * units.fL
    mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
    counts = {
        'flagella': 5,
        'flag_RNA': 30}
    concentrations = {}
    for state_id, count in counts.items():
        concentrations[state_id] = (count / mmol_to_counts).magnitude

    initial_state = {
        'counts': counts,
        'internal': concentrations}
        # 'global': {
        #     'volume': 1.2}}

    return  {
        'transcription_rates': transcription,
        'translation_rates': translation,
        'degradation_rates': degradation,
        'protein_map': protein_map,
        'initial_state': initial_state}

def test_expression(time=10):
    expression_config = get_lacy_config()

    # load process
    expression = ODE_expression(expression_config)

    settings = {
        'total_time': time,
        # 'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-12,
    }

    compartment = process_in_compartment(expression)
    return simulate_with_environment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'ode_expression_process')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    timeseries = test_expression(2520) # 2520 sec (42 min) is the expected doubling time in minimal media
    plot_simulation_output(timeseries, {}, out_dir)
