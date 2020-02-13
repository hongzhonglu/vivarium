from __future__ import absolute_import, division, print_function

import os
import copy

from vivarium.actor.process import Process
from vivarium.actor.composition import simulate_process, convert_to_timeseries, plot_simulation_output
from vivarium.data.nucleotides import nucleotides
from vivarium.utils.units import units

def all_subkeys(d):
    subkeys = set([])
    for ase in d.keys():
        subkeys = subkeys.union(set(d[ase].keys()))
    return list(subkeys)

def kinetics(E, S, kcat, km):
    return kcat * E * S / (S + km)

def keys_list(d):
    return list(d.keys())

DEFAULT_TRANSCRIPT_DEGRADATION_KM = 1e-23

default_degradation_parameters = {
    'sequences': {
        'oA': 'GCC',
        'oAZ': 'GCCGUGCAC',
        'oB': 'AGUUGA',
        'oBY': 'AGUUGACGG'},

    'catalysis_rates': {
        'endoRNAse': 2.0},

    'degradation_rates': {
        'transcripts': {
            'endoRNAse': {
                'oA': DEFAULT_TRANSCRIPT_DEGRADATION_KM,
                'oAZ': DEFAULT_TRANSCRIPT_DEGRADATION_KM,
                'oB': DEFAULT_TRANSCRIPT_DEGRADATION_KM,
                'oBY': DEFAULT_TRANSCRIPT_DEGRADATION_KM}}}}

class RnaDegradation(Process):
    def __init__(self, initial_parameters={}):
        self.default_parameters = default_degradation_parameters

        self.derive_defaults(initial_parameters, 'sequences', 'transcript_order', keys_list)
        self.derive_defaults(initial_parameters, 'catalysis_rates', 'protein_order', keys_list)

        self.parameters = copy.deepcopy(self.default_parameters)
        self.parameters.update(initial_parameters)

        self.sequences = self.parameters['sequences']
        self.catalysis_rates = self.parameters['catalysis_rates']
        self.degradation_rates = self.parameters['degradation_rates']
        self.transcript_order = self.parameters['transcript_order']
        self.protein_order = self.parameters['protein_order']
        self.molecule_order = list(nucleotides.values())

        self.partial_transcripts = {
            transcript: 0
            for transcript in self.transcript_order}

        self.roles = {
            'transcripts': self.transcript_order,
            'proteins': self.protein_order,
            'molecules': self.molecule_order,
            'global': ['mmol_to_counts']}

        super(RnaDegradation, self).__init__(self.roles, self.parameters)

    def default_settings(self):
        default_state = {
            'transcripts': {
                transcript: 10
                for transcript in self.transcript_order},
            'proteins': {
                protein: 1
                for protein in self.protein_order},
            'molecules': {
                nucleotide: 100
                for nucleotide in self.molecule_order}}

        default_emitter_keys = {
            'transcripts': self.transcript_order,
            'proteins': self.protein_order,
            'molecules': self.molecule_order,
            'global': []}

        default_updaters = {
            'transcripts': {},
            'proteins': {},
            'molecules': {},
            'global': {}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'parameters': self.parameters}

    def next_update(self, timestep, states):
        transcripts = states['transcripts']
        proteins = states['proteins']
        molecules = states['molecules']
        mmol_to_counts = states['global']['mmol_to_counts'] * units.L / units.mmol

        delta_transcripts = {
            transcript: 0
            for transcript in self.transcript_order}

        for protein, kcat in self.catalysis_rates.items():
            for transcript, km in self.degradation_rates['transcripts'][protein].items():
                km *= units.mole / units.fL
                delta_transcripts[transcript] += kinetics(
                    proteins[protein] / mmol_to_counts,
                    transcripts[transcript] / mmol_to_counts,
                    kcat,
                    km)

        degradation_levels = {
            transcript: (
                level * mmol_to_counts * timestep).magnitude + self.partial_transcripts[transcript]
            for transcript, level in delta_transcripts.items()}

        transcript_counts = {
            transcript: -int(level)
            for transcript, level in degradation_levels.items()}

        self.partial_transcripts = {
            transcript: level - int(level)
            for transcript, level in degradation_levels.items()}

        delta_molecules = {
            molecule: 0
            for molecule in self.molecule_order}

        for transcript, count in transcript_counts.items():
            sequence = self.sequences[transcript]
            for base in sequence:
                delta_molecules[nucleotides[base]] -= count

        return {
            'transcripts': transcript_counts,
            'molecules': delta_molecules}


def test_rna_degradation(end_time=10):
    rna_degradation = RnaDegradation({})
    settings = {
        'timestep': 1,
        'total_time': end_time}
    return simulate_process(rna_degradation, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'degradation')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}

    saved_data = test_rna_degradation()
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
