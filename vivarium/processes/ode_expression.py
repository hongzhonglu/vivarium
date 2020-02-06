from __future__ import absolute_import, division, print_function

import os
from scipy import constants

from vivarium.actor.process import Process, convert_to_timeseries, \
    plot_simulation_output, simulate_process_with_environment
from vivarium.utils.dict_utils import deep_merge
from vivarium.utils.units import units

default_step_size = 1



class ODE_expression(Process):
    '''
    a ode-based mRNA and protein expression process

    TODO -- add regulation
    '''
    def __init__(self, initial_parameters={}):

        self.nAvogadro = constants.N_A * 1 / units.mol

        self.transcription = initial_parameters.get('transcription_rates', {})
        self.translation = initial_parameters.get('translation_rates', {})
        self.degradation = initial_parameters.get('degradation_rates', {})
        self.protein_map = initial_parameters.get('protein_map', {})

        # get initial state
        states = list(self.transcription.keys()) + list(self.translation.keys())
        null_states = {'internal': {
            state_id: 0 for state_id in states}}
        initialized_states = initial_parameters.get('initial_state', {})
        self.initial_state = deep_merge(null_states, initialized_states)

        # TODO -- get initial counts....

        roles = {
            'counts': states,
            'internal': list(self.initial_state.get('internal', {}).keys()) + ['volume'],
            'external': list(self.initial_state.get('external', {}).keys())}

        parameters = {}
        parameters.update(initial_parameters)


        super(ODE_expression, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        default_state = self.initial_state

        # default emitter keys
        default_emitter_keys = {}

        # default updaters
        default_updaters = {}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        internal = states['internal']
        volume = internal['volume'] * units.fL
        mmol_to_count = self.nAvogadro.to('1/mmol') * volume

        internal_update = {}
        ## transcription
        # transcription rate abstracts over a number of factors, including gene copy number,
        # RNA polymerase abundance, strength of gene protomoters, and availability of nucleotides
        for transcript, rate in self.transcription.items():
            internal_update[transcript] = \
                rate - self.degradation.get(transcript, 0) * internal[transcript]

        ## translation
        # translation rate abstracts over a number of factors, including ribosome availability,
        # strength of mRNA's ribosome binding, availability of tRNAs and free amino acids
        for protein, rate in self.translation.items():
            transcript = self.protein_map[protein]
            transcript_state = internal[transcript]
            internal_update[protein] = \
                rate * transcript_state - self.degradation.get(protein, 0) * internal[protein]

        # convert concentrations to counts
        counts_update = {}
        for mol_id, delta_conc in internal_update.items():
            delta_counts = int((delta_conc * mmol_to_count).magnitude)
            counts_update[mol_id] = delta_counts

        return {
            'internal': internal_update,
            'counts': counts_update}



# test functions
# toy config
toy_transcription_rates = {
    'transcript1': 1.0}

toy_translation_rates = {
    'protein1': 1.0}

toy_protein_map = {
    'protein1': 'transcript1'}

toy_degradation_rates = {
    'transcript1': 0.1,
    'protein1': 0.01}

initial_state = {
    'internal': {
        'volume': 1.2,
        'transcript1': 1.0,
        # 'protein1': 1.0
    }}

def test_expression():
    expression_config = {
        'transcription_rates': toy_transcription_rates,
        'translation_rates': toy_translation_rates,
        'degradation_rates': toy_degradation_rates,
        'protein_map': toy_protein_map,
        'initial_state': initial_state
    }

    # load process
    expression = ODE_expression(expression_config)

    settings = {
        'total_time': 20,
        # 'exchange_role': 'exchange',
        'environment_role': 'external',
        'environment_volume': 1e-12}

    saved_data = simulate_process_with_environment(expression, settings)

    return saved_data



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'ode_expression')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    saved_data = test_expression()
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, {}, out_dir)
    