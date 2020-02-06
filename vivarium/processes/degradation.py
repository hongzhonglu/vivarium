import copy

from vivarium.actor.process import Process
from vivarium.data.nucleotides import nucleotides
from vivarium.processes.deriver import AVOGADRO
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

class RnaDegradation(Process):
    def __init__(self, initial_parameters={}):
        self.default_km = 1e-23
        self.default_parameters = {

            'sequences': {
                'oA': 'GCC',
                'oAZ': 'GCCGUGCAC',
                'oB': 'AGUUGA',
                'oBY': 'AGUUGACGG'},

            'catalysis_rates': {
                'endoRNAse': 10.0},

            'degradation_rates': {
                'transcripts': {
                    'endoRNAse': {
                        'oA': self.default_km,
                        'oAZ': self.default_km,
                        'oB': self.default_km,
                        'oBY': self.default_km}}}}

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

        self.roles = {
            'transcripts': self.transcript_order,
            'proteins': self.protein_order,
            'molecules': self.molecule_order,
            'global': ['volume']}

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
                for nucleotide in self.molecule_order},
            'global': {
                'volume': 1.2}}

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
        volume = states['global']['volume']

        delta_transcripts = {
            transcript: 0
            for transcript in self.transcript_order}

        mmol_to_count = (AVOGADRO * volume * units.fL).to(units.L / units.mmol)

        for protein, kcat in self.catalysis_rates.items():
            for transcript, km in self.degradation_rates['transcripts'][protein].items():
                km *= units.mole / units.fL
                delta_transcripts[transcript] += kinetics(
                    proteins[protein] / mmol_to_count,
                    transcripts[transcript] / mmol_to_count,
                    kcat,
                    km)

        delta_transcripts = {
            transcript: -int((level * mmol_to_count * timestep).magnitude)
            for transcript, level in delta_transcripts.items()}

        delta_molecules = {
            molecule: 0
            for molecule in self.molecule_order}

        for transcript, count in delta_transcripts.items():
            sequence = self.sequences[transcript]
            for base in sequence:
                delta_molecules[nucleotides[base]] += count

        return {
            'transcripts': delta_transcripts,
            'molecules': delta_molecules}


def test_rna_degradation():
    rna_degradation = RnaDegradation({})
    states = rna_degradation.default_settings()['state']
    print('states: {}'.format(states))
    update = rna_degradation.next_update(
        1.0, states)
    print('update: {}'.format(update))


if __name__ == '__main__':
    test_rna_degradation()
