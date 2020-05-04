import copy
import numpy as np
import random
from arrow import StochasticSystem

from vivarium.compartment.process import Process
from vivarium.data.amino_acids import amino_acids
from vivarium.data.molecular_weight import molecular_weight
from vivarium.utils.datum import Datum
from vivarium.utils.polymerize import Elongation, Polymerase, Template, build_stoichiometry, all_products, generate_template
from vivarium.utils.dict_utils import deep_merge

class Ribosome(Polymerase):
    pass

class Transcript(Template):
    pass

def shuffle(l):
    l = [item for item in l]
    np.random.shuffle(l)
    return l

def random_string(alphabet, length):
    string = ''
    for step in range(length):
        string += random.choice(alphabet)
    return string

VERBOSE = True
UNBOUND_RIBOSOME_KEY = 'Ribosome'

monomer_symbols = []
monomer_ids = []

for symbol, id in amino_acids.items():
    monomer_symbols.append(symbol)
    monomer_ids.append(id)

A = random_string(monomer_symbols, 20)
Z = random_string(monomer_symbols, 60)
B = random_string(monomer_symbols, 30)
Y = random_string(monomer_symbols, 40)

default_translation_parameters = {

    'sequences': {
        ('oA', 'eA'): A,
        ('oAZ', 'eA'): A,
        ('oAZ', 'eZ'): Z,
        ('oB', 'eB'): B,
        ('oBY', 'eB'): B,
        ('oBY', 'eY'): Y},

    'templates': {
        ('oA', 'eA'): generate_template(('oA', 'eA'), 20, ['eA']),
        ('oAZ', 'eA'): generate_template(('oAZ', 'eA'), 20, ['eA']),
        ('oAZ', 'eZ'): generate_template(('oAZ', 'eZ'), 60, ['eZ']),
        ('oB', 'eB'): generate_template(('oB', 'eB'), 30, ['eB']),
        ('oBY', 'eB'): generate_template(('oBY', 'eB'), 30, ['eB']),
        ('oBY', 'eY'): generate_template(('oBY', 'eY'), 40, ['eY'])},

    'transcript_affinities': {
        ('oA', 'eA'): 1.0,
        ('oAZ', 'eA'): 2.0,
        ('oAZ', 'eZ'): 5.0,
        ('oB', 'eB'): 1.0,
        ('oBY', 'eB'): 2.0,
        ('oBY', 'eY'): 5.0},

    'elongation_rate': 5.0,
    'polymerase_occlusion': 10,
    'symbol_to_monomer': amino_acids,
    'monomer_ids': monomer_ids,
    'concentration_keys': []}

def gather_genes(affinities):
    genes = {}
    for operon, product in affinities.keys():
        if not operon in genes:
            genes[operon] = []
        genes[operon].append(product)
    return genes

def transcripts_to_gene_counts(transcripts, operons):
    counts = {}
    for transcript, genes in operons.items():
        for gene in genes:
            counts[(transcript, gene)] = transcripts.get(transcript, 0)
    return counts

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
        self.default_parameters['molecule_ids'] = self.monomer_ids

        self.parameters = copy.deepcopy(self.default_parameters)
        self.parameters.update(initial_parameters)

        self.sequences = self.parameters['sequences']
        self.templates = self.parameters['templates']

        self.transcript_affinities = self.parameters['transcript_affinities']
        self.operons = gather_genes(self.transcript_affinities)
        self.operon_order = list(self.operons.keys())
        self.transcript_order = self.parameters['transcript_order']
        self.transcript_count = len(self.transcript_order)

        self.monomer_ids = self.parameters['monomer_ids']
        self.molecule_ids = self.parameters['molecule_ids']
        self.protein_ids = self.parameters['protein_ids']
        self.symbol_to_monomer = self.parameters['symbol_to_monomer']
        self.elongation = 0
        self.elongation_rate = self.parameters['elongation_rate']
        self.polymerase_occlusion = self.parameters['polymerase_occlusion']
        self.concentration_keys = self.parameters['concentration_keys']

        self.affinity_vector = np.array([
            self.transcript_affinities[transcript_key]
            for transcript_key in self.transcript_order], dtype=np.float64)

        self.stoichiometry = build_stoichiometry(self.transcript_count)

        self.initiation = StochasticSystem(self.stoichiometry)

        self.ribosome_id = 0

        concentration_keys = self.concentration_keys + self.protein_ids
        self.ports = {
            'ribosomes': ['ribosomes'],
            'molecules': self.molecule_ids,
            'transcripts': list(self.operons.keys()),
            'proteins': concentration_keys + [UNBOUND_RIBOSOME_KEY],
            'concentrations': concentration_keys,
            'global': []}

        if VERBOSE:
            print('translation parameters: {}'.format(self.parameters))

        super(Translation, self).__init__(self.ports, self.parameters)

    def default_settings(self):
        default_state = {
            'ribosomes': {
                'ribosomes': []},
            'molecules': {},
            'transcripts': {
                transcript_id: 1
                for transcript_id in self.operons.keys()},
            'proteins': dict({
                UNBOUND_RIBOSOME_KEY: 10})}

        default_state['proteins'].update({
            protein_id: 0
            for protein_id in self.protein_ids})

        default_state['molecules'].update({
            monomer_id: 200
            for monomer_id in self.monomer_ids})

        default_state = deep_merge(
            default_state,
            self.parameters.get('initial_state', {}))

        operons = list(default_state['transcripts'].keys())
        default_emitter_keys = {
            'ribosomes': ['ribosomes'],
            'molecules': self.monomer_ids,
            'transcripts': operons,
            'proteins': self.protein_ids + [UNBOUND_RIBOSOME_KEY]}

        # schema
        mols_with_mass = [
            mol_id for mol_id in self.ports['proteins']
            if mol_id in molecular_weight]
        schema = {
            'ribosomes': {
                'ribosomes': {'updater': 'set'}},
            'proteins': {mol_id: {
                    'mass': molecular_weight.get(mol_id)}
                for mol_id in mols_with_mass}}

        # deriver_settings
        deriver_setting = [{
            'type': 'mass',
            'source_port': 'proteins',
            'derived_port': 'global',
            'keys': mols_with_mass},
            {
            'type': 'counts_to_mmol',
            'source_port': 'proteins',
            'derived_port': 'concentrations',
            'keys': self.protein_ids + self.concentration_keys
            }]

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'deriver_setting': deriver_setting,
            'parameters': self.parameters}

    def next_update(self, timestep, states):
        molecules = states['molecules']
        transcripts = states['transcripts']
        proteins = states['proteins']
        ribosomes = list(map(Ribosome, states['ribosomes']['ribosomes']))

        gene_counts = np.array(
            list(transcripts_to_gene_counts(transcripts, self.operons).values()),
            dtype=np.int64)

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

        original_unbound_ribosomes = proteins[UNBOUND_RIBOSOME_KEY]
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
                gene_counts - bound_transcripts,
                bound_transcripts,
                [unbound_ribosomes]])

            # find number of monomers until next terminator
            distance = 1 / self.elongation_rate

            # find interval of time that elongates to the point of the next terminator
            interval = min(distance, timestep - time)

            if interval == distance:
                # perform the elongation until the next event
                terminations, monomer_limits, ribosomes = elongation.step(
                    interval,
                    monomer_limits,
                    ribosomes)
                unbound_ribosomes += terminations
            else:
                elongation.store_partial(interval)
                terminations = 0

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(
                interval,
                substrate,
                self.affinity_vector)

            # go through each event in the simulation and update the state
            ribosome_bindings = 0
            for now, event in zip(result['time'], result['events']):
                # ribosome has bound the transcript
                transcript_key = self.transcript_order[event]
                bound_transcripts[event] += 1

                self.ribosome_id += 1
                new_ribosome = Ribosome({
                    'id': self.ribosome_id,
                    'template': transcript_key,
                    'position': 0})
                new_ribosome.bind()
                new_ribosome.start_polymerizing()
                ribosomes.append(new_ribosome)

                ribosome_bindings += 1
                unbound_ribosomes -= 1

            # deal with occluding rnap
            for ribosome in ribosomes:
                if ribosome.is_unoccluding(self.polymerase_occlusion):
                    bound_transcripts[ribosome.template_index] -= 1
                    ribosome.unocclude()

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        proteins = {
            UNBOUND_RIBOSOME_KEY: unbound_ribosomes - original_unbound_ribosomes}
        proteins.update(elongation.complete_polymers)

        molecules = {
            key: count * -1
            for key, count in elongation.monomers.items()}

        update = {
            'ribosomes': {
                'ribosomes': [ribosome.to_dict() for ribosome in ribosomes]},
            'molecules': molecules,
            'proteins': proteins}

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
        'molecules': {},
        'proteins': {UNBOUND_RIBOSOME_KEY: 10},
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
        
