import copy
import numpy as np
import logging as log
from arrow import StochasticSystem

from vivarium.utils.dict_utils import deep_merge
from vivarium.compartment.process import Process, keys_list
from vivarium.states.chromosome import Chromosome, Rnap, Promoter, frequencies, add_merge, toy_chromosome_config
from vivarium.utils.polymerize import Elongation, build_stoichiometry, template_products
from vivarium.data.nucleotides import nucleotides

# log.basicConfig(level=log.DEBUG)

def choose_element(elements):
    if elements:
        choice = np.random.choice(len(elements), 1)
        return list(elements)[int(choice)]

UNBOUND_RNAP_KEY = 'RNA Polymerase'

monomer_ids = list(nucleotides.values())

default_transcription_parameters = {
    'promoter_affinities': {
        ('pA', None): 1.0,
        ('pA', 'tfA'): 10.0,
        ('pB', None): 1.0,
        ('pB', 'tfB'): 10.0},
    'transcription_factors': ['tfA', 'tfB'],
    'sequence': toy_chromosome_config['sequence'],
    'templates': toy_chromosome_config['promoters'],
    'genes': toy_chromosome_config['genes'],
    'elongation_rate': 1.0,
    'polymerase_occlusion': 5,
    'symbol_to_monomer': nucleotides,
    'monomer_ids': monomer_ids,
    'molecule_ids': monomer_ids}

class Transcription(Process):
    def __init__(self, initial_parameters={}):
        '''
        The main parameters here are:
        * promoter_affinities - a dict where keys are tuples representing promoter
            states with respect to their transcription factor binding sites, and the
            values are the affinity RNAP have for that promoter state.
        * promoter_order - a list representing a canonical ordering of the promoters.
        * elongation_rate - elongation rate of polymerizing mRNA.
        '''

        log.debug('inital_parameters: {}'.format(initial_parameters))

        self.default_parameters = default_transcription_parameters
        self.derive_defaults(initial_parameters, 'templates', 'promoter_order', keys_list)
        self.derive_defaults(initial_parameters, 'templates', 'transcript_ids', template_products)

        self.parameters = copy.deepcopy(self.default_parameters)
        self.parameters.update(initial_parameters)

        self.sequence = self.parameters['sequence']
        self.templates = self.parameters['templates']
        self.genes = self.parameters['genes']
        empty_chromosome = Chromosome({
            'sequence': self.sequence,
            'promoters': self.templates,
            'genes': self.genes})
        self.sequences = empty_chromosome.sequences()
        self.symbol_to_monomer = self.parameters['symbol_to_monomer']

        log.debug('chromosome sequence: {}'.format(self.sequence))

        self.promoter_affinities = self.parameters['promoter_affinities']
        self.promoter_order = self.parameters['promoter_order']
        self.promoter_count = len(self.promoter_order)

        self.transcription_factors = self.parameters['transcription_factors']
        self.molecule_ids = self.parameters['molecule_ids']
        self.monomer_ids = self.parameters['monomer_ids']
        self.transcript_ids = self.parameters['transcript_ids']
        self.elongation = 0
        self.elongation_rate = self.parameters['elongation_rate']
        self.polymerase_occlusion = self.parameters['polymerase_occlusion']

        self.stoichiometry = build_stoichiometry(self.promoter_count)
        self.initiation = StochasticSystem(self.stoichiometry, random_seed=np.random.randint(2**31))

        self.ports = {
            'chromosome': ['rnaps', 'rnap_id', 'domains', 'root_domain'],
            'molecules': self.molecule_ids,
            'factors': self.transcription_factors,
            'transcripts': self.transcript_ids,
            'proteins': [UNBOUND_RNAP_KEY] + self.transcription_factors}

        log.debug('transcription parameters: {}'.format(self.parameters))

        super(Transcription, self).__init__(self.ports, self.parameters)

    def build_affinity_vector(self, promoters, factors):
        vector = np.zeros(len(self.promoter_order), dtype=np.float64)
        for index, promoter_key in enumerate(self.promoter_order):
            promoter = promoters[promoter_key]
            binding = promoter.binding_state(factors)
            affinity = self.promoter_affinities.get(binding, 0.0)
            # print('promoter state - {}: {}'.format(binding, affinity))
            vector[index] = affinity
        return vector

    def chromosome_config(self, chromosome_states):
        return dict(
            chromosome_states,
            sequence=self.sequence,
            promoters=self.templates,
            promoter_order=self.promoter_order,
            genes=self.genes)

    def default_settings(self):
        default_state = {
            'chromosome': {
                'rnaps': [],
                'rnap_id': 0,
                'root_domain': 0,
                'domains': {
                    0: {
                        'id': 0,
                        'lead': 0,
                        'lag': 0,
                        'children': []}}},
            'molecules': {},
            'proteins': {UNBOUND_RNAP_KEY: 10},
            'factors': {
                key: 0.0
                for key in self.transcription_factors}}

        default_state['molecules'].update({
            nucleotide: 100
            for nucleotide in self.monomer_ids})

        chromosome = Chromosome(
            self.chromosome_config(
                default_state['chromosome']))

        operons = [operon.id for operon in chromosome.operons()]

        default_state['transcripts'] = {
            operon: 0
            for operon in operons}

        default_state = deep_merge(
            default_state,
            self.parameters.get('initial_state', {}))

        default_emitter_keys = {
            'chromosome': ['rnaps'],
            'molecules': self.monomer_ids,
            'proteins': [UNBOUND_RNAP_KEY],
            'transcripts': operons}

        schema = {
            'chromosome': {
                state_id : {
                    'updater': 'set'}
                for state_id in self.ports['chromosome']}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema,
            'parameters': self.parameters}

    def next_update(self, timestep, states):
        chromosome = Chromosome(
            self.chromosome_config(
                states['chromosome']))
        molecules = states['molecules']
        proteins = states['proteins']
        factors = states['factors'] # as concentrations

        # print('concentrations: flhDC - {}, fliA - {}'.format(
        #     factors['flhDC'],
        #     factors['fliA']))

        promoter_rnaps = chromosome.promoter_rnaps()
        promoter_domains = chromosome.promoter_domains()

        # Find out how many promoters are currently blocked by a
        # newly initiated or occluding rnap
        promoter_count = len(chromosome.promoter_order)
        blocked_promoters = np.zeros(promoter_count, dtype=np.int64)
        open_domains = {}
        bound_domains = {}
        for promoter_index, promoter_key in enumerate(chromosome.promoter_order):
            domains = []
            for rnap in promoter_rnaps.get(promoter_key, {}).values():
                if rnap.is_occluding():
                    domains.append(rnap.domain)
                    blocked_promoters[promoter_index] += 1

            bound_domains[promoter_key] = set(domains)
            open_domains[promoter_key] = promoter_domains[promoter_key] - bound_domains[promoter_key]

        blocked_promoters = np.array(blocked_promoters)

        # Make the state for a gillespie simulation out of total number of each
        # promoter by copy number not blocked by initiated rnap,
        # concatenated with the number of each promoter that is bound by rnap.
        # These are the two states for each promoter the simulation
        # will operate on, essentially going back and forth between
        # bound and unbound states.
        copy_numbers = chromosome.promoter_copy_numbers()
        original_unbound_rnaps = proteins[UNBOUND_RNAP_KEY]
        monomer_limits = {
            monomer: molecules[monomer]
            for monomer in self.monomer_ids}
        unbound_rnaps = original_unbound_rnaps

        time = 0
        now = 0
        elongation = Elongation(
            self.sequences,
            chromosome.promoters,
            monomer_limits,
            self.symbol_to_monomer,
            self.elongation)

        initiation_affinity = self.build_affinity_vector(chromosome.promoters, factors)

        # print('initiation affinity - {}'.format(initiation_affinity))

        while time < timestep:
            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                copy_numbers - blocked_promoters,
                blocked_promoters,
                [unbound_rnaps]])

            log.debug('transcription substrate: {}'.format(substrate))
            log.debug('blocked promoters: {}'.format(blocked_promoters))

            # find number of monomers until next terminator
            distance = 1 / self.elongation_rate # chromosome.terminator_distance()

            # find interval of time that elongates to the point of the next terminator
            interval = min(distance, timestep - time)

            if interval == distance:
                # perform the elongation until the next event
                terminations, monomer_limits, chromosome.rnaps = elongation.step(
                    interval,
                    monomer_limits,
                    chromosome.rnaps)
                unbound_rnaps += terminations
            else:
                elongation.store_partial(interval)
                terminations = 0

            log.debug('time: {} --- interval: {}'.format(time, interval))
            log.debug('monomer limits: {}'.format(monomer_limits))
            log.debug('terminations: {}'.format(terminations))

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(
                interval,
                substrate,
                initiation_affinity)

            log.debug('result: {}'.format(result))

            # perform binding
            for now, event in zip(result['time'], result['events']):
                # RNAP has bound the promoter
                promoter_key = chromosome.promoter_order[event]
                promoter = chromosome.promoters[promoter_key]
                domains = open_domains[promoter_key]
                domain = choose_element(domains)

                blocked_promoters[event] += 1
                bound_domains[promoter_key].add(domain)
                open_domains[promoter_key].remove(domain)

                # create a new bound RNAP and add it to the chromosome.
                new_rnap = chromosome.bind_rnap(event, domain)
                new_rnap.start_polymerizing()

                log.debug('newly bound RNAP: {}'.format(new_rnap))

                unbound_rnaps -= 1

            # deal with occluding rnap
            for rnap in chromosome.rnaps:
                if rnap.is_unoccluding(self.polymerase_occlusion):
                    log.debug('RNAP unoccluding: {}'.format(rnap))

                    blocked_promoters[rnap.template_index] -= 1
                    bound_domains[rnap.template].remove(rnap.domain)
                    open_domains[rnap.template].add(rnap.domain)
                    rnap.unocclude()
                log.debug('rnap: {}'.format(rnap))

            log.debug('complete: {}'.format(elongation.complete_polymers))

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        proteins = {
            UNBOUND_RNAP_KEY: unbound_rnaps - original_unbound_rnaps}

        molecules = {
            key: count * -1
            for key, count in elongation.monomers.items()}

        chromosome_dict = chromosome.to_dict()
        update = {
            'chromosome': {
                key: chromosome_dict[key]
                for key in self.ports['chromosome']},
            'proteins': proteins,
            'molecules': molecules,
            'transcripts': elongation.complete_polymers}

        log.debug('molecules update: {}'.format(update['molecules']))

        return update


def test_transcription():
    parameters = {
        'elongation_rate': 10.0}

    chromosome = Chromosome(toy_chromosome_config)
    transcription = Transcription(parameters)

    states = {
        'chromosome': chromosome.to_dict(),
        'molecules': {},
        'proteins': {UNBOUND_RNAP_KEY: 10},
        'factors': {'tfA': 0.2, 'tfB': 0.7}}

    states['molecules'].update({
        nucleotide: 10
        for nucleotide in transcription.monomer_ids})

    update = transcription.next_update(1.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_transcription()
