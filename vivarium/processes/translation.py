import copy
import numpy as np
import random
from arrow import StochasticSystem

from vivarium.actor.process import Process
from vivarium.data.amino_acids import amino_acids
from vivarium.utils.datum import Datum
from vivarium.utils.polymerize import Elongation, Polymerase, Template, build_stoichiometry, build_rates, all_products

class Ribosome(Polymerase):
    pass

class Transcript(Template):
    pass

def generate_template(id, length, products):
    return {
        'id': id,
        'position': 0,
        'direction': 1,
        'sites': [],
        'terminators': [
            {'position': length,
             'strength': 1.0,
             'products': products}]}

def shuffle(l):
    l = [item for item in l]
    np.random.shuffle(l)
    return l

def random_string(alphabet, length):
    string = ''
    for step in range(length):
        string += random.choice(alphabet)
    return string

UNBOUND_RIBOSOME_KEY = 'Ribosome'

monomer_symbols = []
monomer_ids = []

for symbol, id in amino_acids.items():
    monomer_symbols.append(symbol)
    monomer_ids.append(id)

default_translation_parameters = {
    'sequences': {
        'oA': random_string(monomer_symbols, 20),
        'oAZ': random_string(monomer_symbols, 50),
        'oB': random_string(monomer_symbols, 30),
        'oBY': random_string(monomer_symbols, 70)},
    'templates': {
        'oA': generate_template('oA', 20, ['eA']),
        'oAZ': generate_template('oAZ', 50, ['eA', 'eZ']),
        'oB': generate_template('oB', 30, ['eB']),
        'oBY': generate_template('oBY', 70, ['eB', 'eY'])},
    'transcript_affinities': {
        'oA': 1.0,
        'oAZ': 1.0,
        'oB': 1.0,
        'oBY': 1.0},
    'elongation_rate': 5.0,
    'advancement_rate': 10.0,
    'symbol_to_monomer': amino_acids,
    'monomer_ids': monomer_ids}

class Translation(Process):
    def __init__(self, initial_parameters={}):
        self.monomer_symbols = list(amino_acids.keys())
        self.monomer_ids = list(amino_acids.values())

        self.default_parameters = default_translation_parameters

        templates = initial_parameters.get(
            'templates',
            self.default_parameters['templates'])

        self.default_parameters['protein_ids'] = all_products({
            key: Template(config)
            for key, config in templates.items()})

        self.default_parameters['transcript_order'] = list(
            initial_parameters.get(
                'transcript_affinities',
                self.default_parameters['transcript_affinities']).keys())
        self.default_parameters['molecule_ids'] = self.monomer_ids + [
            UNBOUND_RIBOSOME_KEY]

        self.parameters = copy.deepcopy(self.default_parameters)
        self.parameters.update(initial_parameters)

        self.sequences = self.parameters['sequences']
        self.templates = self.parameters['templates']

        self.transcript_affinities = self.parameters['transcript_affinities']
        self.transcript_order = self.parameters['transcript_order']
        self.transcript_count = len(self.transcript_order)

        self.monomer_ids = self.parameters['monomer_ids']
        self.molecule_ids = self.parameters['molecule_ids']
        self.protein_ids = self.parameters['protein_ids']
        self.symbol_to_monomer = self.parameters['symbol_to_monomer']
        self.elongation = 0
        self.elongation_rate = self.parameters['elongation_rate']
        self.advancement_rate = self.parameters['advancement_rate']

        self.affinity_vector = np.array([
            self.transcript_affinities[transcript_key]
            for transcript_key in self.transcript_order], dtype=np.float64)

        self.stoichiometry = build_stoichiometry(self.transcript_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        self.initiation = StochasticSystem(self.stoichiometry, self.rates)

        self.ribosome_id = 0

        self.roles = {
            'ribosomes': ['ribosomes'],
            'molecules': self.molecule_ids,
            'transcripts': self.transcript_order,
            'proteins': self.protein_ids}

        print('translation parameters: {}'.format(self.parameters))

        super(Translation, self).__init__(self.roles, self.parameters)

    def default_settings(self):
        default_state = {
            'ribosomes': {
                'ribosomes': []},
            'molecules': dict({
                UNBOUND_RIBOSOME_KEY: 10}),
            'transcripts': {
                transcript_id: 1
                for transcript_id in self.transcript_order},
            'proteins': {
                protein_id: 0
                for protein_id in self.protein_ids}}

        default_state['molecules'].update({
            monomer_id: 200
            for monomer_id in self.monomer_ids})

        operons = list(default_state['transcripts'].keys())
        default_emitter_keys = {
            'ribosomes': ['ribosomes'],
            'molecules': self.monomer_ids + [UNBOUND_RIBOSOME_KEY],
            'transcripts': operons,
            'proteins': self.protein_ids}

        default_updaters = {
            'ribosomes': {'ribosomes': 'set'},
            'molecules': {},
            'transcripts': {},
            'proteins': {}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'parameters': self.parameters}

    def next_update(self, timestep, states):
        ribosomes = list(map(Ribosome, states['ribosomes']['ribosomes']))
        molecules = states['molecules']
        transcripts = states['transcripts']
        transcript_counts = np.array([
            transcripts.get(transcript_key, 0)
            for transcript_key in self.transcript_order], dtype=np.int64)

        # Find out how many transcripts are currently blocked by a
        # newly initiated ribosome
        bound_transcripts = np.zeros(self.transcript_count, dtype=np.int64)
        ribosomes_by_transcript = {
            transcript_key: []
            for transcript_key in self.transcript_order}
        for ribosome in ribosomes:
            ribosomes_by_transcript[ribosome.template].append(ribosome)
        for index, transcript in enumerate(self.transcript_order):
            bound_transcripts[index] = len([
                ribosome
                for ribosome in ribosomes_by_transcript[transcript]
                if ribosome.is_bound()])

        # Make the state for a gillespie simulation out of total number of each
        # transcript not blocked by a bound ribosome, concatenated with the number
        # of each transcript that is bound by a ribosome.
        # These are the two states for each transcript the simulation
        # will operate on, essentially going back and forth between
        # bound and unbound states.

        original_unbound_ribosomes = molecules[UNBOUND_RIBOSOME_KEY]
        monomer_limits = {
            monomer: molecules[monomer]
            for monomer in self.monomer_ids}
        unbound_ribosomes = original_unbound_ribosomes

        templates = {
            key: Template(template)
            for key, template in self.templates.items()}

        time = 0
        now = 0
        elongation = Elongation(
            self.sequences,
            templates,
            monomer_limits,
            self.symbol_to_monomer,
            self.elongation)

        while time < timestep:
            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                transcript_counts - bound_transcripts,
                bound_transcripts,
                [unbound_ribosomes]])

            # find number of monomers until next terminator
            # distance = chromosome.terminator_distance()
            distance = 1

            # find interval of time that elongates to the point of the next terminator
            interval = min(distance / self.elongation_rate, timestep - time)

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(interval, substrate)

            # go through each event in the simulation and update the state
            ribosome_bindings = 0
            for now, event in zip(result['time'], result['events']):

                # perform the elongation until the next event
                terminations, monomer_limits, ribosomes = elongation.elongate(
                    time + now,
                    self.elongation_rate,
                    monomer_limits,
                    ribosomes)
                unbound_ribosomes += terminations

                # ribosome has bound the transcript
                if event < self.transcript_count:
                    transcript_key = self.transcript_order[event]
                    # transcript = self.templates[transcript_key]
                    bound_transcripts[event] += 1

                    self.ribosome_id += 1
                    new_ribosome = Ribosome({
                        'id': self.ribosome_id,
                        'template': transcript_key,
                        'position': 0})
                    new_ribosome.bind()
                    ribosomes.append(new_ribosome)
                    ribosomes_by_transcript[transcript_key].append(new_ribosome)

                    ribosome_bindings += 1
                    unbound_ribosomes -= 1
                # ribosome has begun polymerizing its protein
                else:
                    transcript_index = event - self.transcript_count
                    transcript_key = self.transcript_order[transcript_index]

                    bound_transcripts[transcript_index] -= 1

                    ribosome = ribosomes_by_transcript[transcript_key].pop()
                    ribosome.start_transcribing()

            # now that all events have been accounted for, elongate
            # until the end of this interval.
            terminations, monomer_limits, ribosomes = elongation.elongate(
                time + interval,
                self.elongation_rate,
                monomer_limits,
                ribosomes)
            unbound_ribosomes += terminations

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        molecules = {
            UNBOUND_RIBOSOME_KEY: unbound_ribosomes - original_unbound_ribosomes}

        molecules.update({
            key: count * -1
            for key, count in elongation.monomers.items()})

        update = {
            'ribosomes': {
                'ribosomes': [ribosome.to_dict() for ribosome in ribosomes]},
            'molecules': molecules,
            'proteins': elongation.complete_polymers}

        return update


def test_translation():
    parameters = {
        'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
        'elongation_rate': 10.0,
        'advancement_rate': 10.0}

    parameters = {}
    translation = Translation(parameters)

    states = {
        'ribosomes': {'ribosomes': []},
        'molecules': {UNBOUND_RIBOSOME_KEY: 10},
        'transcripts': {
            'oA': 10,
            'oAZ': 10,
            'oB': 10,
            'oBY': 10}}
    states['molecules'].update({
        molecule_id: 100
        for molecule_id in translation.monomer_ids})

    update = translation.next_update(10.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_translation()
        
