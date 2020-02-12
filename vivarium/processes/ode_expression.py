from __future__ import absolute_import, division, print_function

import os
from scipy import constants

from vivarium.actor.process import Process, convert_to_timeseries, \
    plot_simulation_output, simulate_process_with_environment
from vivarium.utils.dict_utils import deep_merge, tuplify_role_dicts
from vivarium.utils.units import units
from vivarium.utils.regulation_logic import build_rule



class ODE_expression(Process):
    '''
    a ode-based mRNA and protein expression process, with boolean logic for regulation

    TODO -- kinetic regulation, cooperativity, autoinhibition, autactivation
    '''
    def __init__(self, initial_parameters={}):

        self.nAvogadro = constants.N_A * 1 / units.mol

        # ode gene expression
        self.transcription = initial_parameters.get('transcription_rates', {})
        self.translation = initial_parameters.get('translation_rates', {})
        self.degradation = initial_parameters.get('degradation_rates', {})
        self.protein_map = initial_parameters.get('protein_map', {})

        # boolean regulation
        regulation_logic = initial_parameters.get('regulation', {})
        self.regulation = {
            gene_id: build_rule(logic) for gene_id, logic in regulation_logic.items()}
        regulators = initial_parameters.get('regulators', [])
        internal_regulators = [state_id for role_id, state_id in regulators if role_id == 'internal']
        external_regulators = [state_id for role_id, state_id in regulators if role_id == 'external']

        # get initial state
        states = list(self.transcription.keys()) + list(self.translation.keys())
        null_states = {'internal': {
            state_id: 0 for state_id in states}}
        initialized_states = initial_parameters.get('initial_state', {})
        self.initial_state = deep_merge(null_states, initialized_states)
        internal = list(self.initial_state.get('internal', {}).keys())
        external = list(self.initial_state.get('external', {}).keys())

        self.partial_expression = {
            mol_id: 0 for mol_id in internal}

        # convert initial counts to concs.
        # TODO -- concs to counts (use deriver?)
        counts = self.initial_state['counts']
        volume = self.initial_state['global'].get('volume', 1.2)
        # concentrations = self.initial_state['internal']
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume
        new_concentrations = {}
        for state_id, count in counts.items():
            conc = (count / mmol_to_count).magnitude
            new_concentrations[state_id] = conc
        self.initial_state['internal'].update(new_concentrations)


        roles = {
            'counts': states,
            'internal': internal + internal_regulators,
            'external': external + external_regulators,
            'global': ['volume']}

        parameters = {}
        parameters.update(initial_parameters)

        super(ODE_expression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state
        default_emitter_keys = self.roles
        default_updaters = {}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        internal_state = states['internal']
        volume = states['global']['volume'] * units.fL
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_role_dicts(states)
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

        # convert concentrations to counts.
        # keep partially expressed molecules for the next iteration
        expression_levels = {
            mol_id: (delta_conc * mmol_to_count).magnitude + self.partial_expression[mol_id] for
            mol_id, delta_conc in internal_update.items()}

        counts_update = {
            mol_id: int(level) for mol_id, level in expression_levels.items()}

        self.partial_expression = {
             mol_id: level - int(level)
             for mol_id, level in expression_levels.items()}

        return {
            'internal': internal_update,
            'counts': counts_update}



# functions
def get_toy_expression():
    # toy config
    transcription = {
        'lacy_RNA': 1e-20}

    translation = {
        'LacY': 1e-2}

    protein_map = {
        'LacY': 'lacy_RNA'}

    degradation = {
        'lacy_RNA': 0.2,
        'LacY': 0.001}

    initial_state = {
        'internal': {
            'lacy_RNA': 0,
            'LacY': 0.0},
        'global': {
            'volume': 1.2}}

    return {
        'transcription_rates': transcription,
        'translation_rates': translation,
        'degradation_rates': degradation,
        'protein_map': protein_map,
        'initial_state': initial_state}


def get_flagella_expression():
    transcription = {
        'flag_RNA': 1e-21}

    translation = {
        'flagella': 8e-5}

    degradation = {
        'flag_RNA': 2e-2}  # 1e-23}

    protein_map = {
        'flagella': 'flag_RNA'}

    initial_state = {
        'counts': {
            'flagella': 5,
            'flag_RNA': 30},
        'internal': {
            'flagella': 0,
            'flag_RNA': 0},
        'global': {
            'volume': 1.2}}

    return  {
        'transcription_rates': transcription,
        'translation_rates': translation,
        'degradation_rates': degradation,
        'protein_map': protein_map,
        'initial_state': initial_state}

default_settings = {
    'total_time': 10,
    # 'exchange_role': 'exchange',
    'environment_role': 'external',
    'environment_volume': 1e-12}

def test_expression(expression_config=get_toy_expression(), settings=default_settings):
    expression = ODE_expression(expression_config)
    saved_data = simulate_process_with_environment(expression, settings)
    return saved_data


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'ode_expression')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    settings = {
        'total_time': 2520,
        # 'exchange_role': 'exchange',
        'environment_role': 'external',
        'environment_volume': 1e-12}

    # saved_data = test_expression()
    saved_data = test_expression(get_flagella_expression(), settings)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)
    