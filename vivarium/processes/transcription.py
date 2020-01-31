import copy
import math
import numpy as np
from arrow import StochasticSystem

from vivarium.actor.process import Process
from vivarium.states.chromosome import Chromosome, Rnap, Promoter, frequencies, add_merge, test_chromosome_config
from vivarium.utils.polymerize import Elongation, build_stoichiometry, build_rates, all_products
from vivarium.data.nucleotides import nucleotides

def choose_element(elements):
    if elements:
        choice = np.random.choice(len(elements), 1)
        return list(elements)[int(choice)]

class Transcription(Process):
    def __init__(self, initial_parameters={}):
        '''
        The main parameters here are:
        * promoter_affinities - a dict where keys are tuples representing promoter
            states with respect to their transcription factor binding sites, and the
            values are the affinity RNAP have for that promoter state.
        * promoter_order - a list representing a canonical ordering of the promoters.
        * elongation_rate - elongation rate of polymerizing mRNA.
        * advancement_rate - affinity for RNAP to move from the 'bound' to the 
            'transcribing' state.
        '''

        print('inital_parameters: {}'.format(initial_parameters))

        # TODO: add monomer_mapping parameter for monomer names

        monomer_ids = list(nucleotides.values())
        self.unbound_rnap_key = 'RNA Polymerase'
        self.default_parameters = {
            'promoter_affinities': {
                'pA': 1.0,
                'pB': 1.0},
            'sequence': test_chromosome_config['sequence'],
            'templates': test_chromosome_config['promoters'],
            'genes': test_chromosome_config['genes'],
            'elongation_rate': 1.0,
            'advancement_rate': 1.0,
            'symbol_to_monomer': nucleotides,
            'monomer_ids': monomer_ids,
            'molecule_ids': monomer_ids + [self.unbound_rnap_key]}

        self.default_parameters['promoter_order'] = list(
            initial_parameters.get(
                'promoter_affinities',
                self.default_parameters['promoter_affinities']).keys())
        self.default_parameters['transcript_ids'] = all_products(
            initial_parameters.get(
                'templates',
                self.default_parameters['templates']))

        parameters = copy.deepcopy(self.default_parameters)
        parameters.update(initial_parameters)

        self.sequence = parameters['sequence']
        self.sequences = None # set when the chromosome first appears
        self.templates = parameters['templates']
        # {
        #     key: Promoter(config)
        #     for key, config in parameters['templates'].items()}

        self.genes = parameters['genes']
        self.symbol_to_monomer = parameters['symbol_to_monomer']

        print('chromosome sequence: {}'.format(self.sequence))

        self.promoter_affinities = parameters['promoter_affinities']
        self.promoter_order = parameters['promoter_order']
        self.promoter_count = len(self.promoter_order)

        self.molecule_ids = parameters['molecule_ids']
        self.monomer_ids = parameters['monomer_ids']
        self.transcript_ids = parameters['transcript_ids']
        self.elongation = 0
        self.elongation_rate = parameters['elongation_rate']
        self.advancement_rate = parameters['advancement_rate']

        self.affinity_vector = np.array([
            self.promoter_affinities[promoter_key]
            for promoter_key in self.promoter_order], dtype=np.float64)

        self.stoichiometry = build_stoichiometry(self.promoter_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        self.initiation = StochasticSystem(self.stoichiometry, self.rates)

        self.roles = {
            'chromosome': Chromosome({}).fields(),
            'molecules': self.molecule_ids,
            'transcripts': self.transcript_ids}

        print('transcription parameters: {}'.format(parameters))

        super(Transcription, self).__init__(self.roles, parameters)

    def default_settings(self):
        default_state = {
            'chromosome': {
                'sequence': self.sequence,
                'promoters': self.templates,
                'genes': self.genes,
                'rnaps': [],
                'domains': {
                    0: {
                        'id': 0,
                        'lead': 0,
                        'lag': 0,
                        'children': []}}},
            'molecules': {self.unbound_rnap_key: 10}}

        default_state['molecules'].update({
            nucleotide: 100
            for nucleotide in self.monomer_ids})

        chromosome = Chromosome(default_state['chromosome'])
        operons = [operon.id for operon in chromosome.operons()]

        default_state['transcripts'] = {
            operon: 0
            for operon in operons}

        default_emitter_keys = {
            'chromosome': ['rnaps'],
            'molecules': list(default_state['molecules'].keys()),
            'transcripts': operons}

        default_updaters = {
            'chromosome': {
                'sequence': 'set',
                'genes': 'set',
                'promoters': 'set',
                'domains': 'set',
                'root_domain': 'set',
                'promoter_order': 'set',
                'rnap_id': 'set',
                'rnaps': 'set'},
            'molecules': {},
            'transcripts': {}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'parameters': self.default_parameters}

    def next_update(self, timestep, states):
        chromosome = Chromosome(states['chromosome'])
        molecules = states['molecules']

        if self.sequences is None:
            self.sequences = chromosome.sequences()
            print('sequences: {}'.format(self.sequences))

        promoter_rnaps = chromosome.promoter_rnaps()
        promoter_domains = chromosome.promoter_domains()

        # Find out how many promoters are currently blocked by a
        # newly initiated rnap
        bound_rnap = []
        open_domains = {}
        bound_domains = {}
        for promoter_key in chromosome.promoter_order:
            bound_domains[promoter_key] = set([
                rnap.domain
                for rnap in promoter_rnaps.get(promoter_key, {}).values()
                if rnap.is_bound()])
            bound_rnap.append(len(bound_domains[promoter_key]))
            open_domains[promoter_key] = promoter_domains[promoter_key] - bound_domains[promoter_key]

        bound_rnap = np.array(bound_rnap)

        # Make the state for a gillespie simulation out of total number of each
        # promoter by copy number not blocked by initiated rnap,
        # concatenated with the number of each promoter that is bound by rnap.
        # These are the two states for each promoter the simulation
        # will operate on, essentially going back and forth between
        # bound and unbound states.
        copy_numbers = chromosome.promoter_copy_numbers()
        original_unbound_rnaps = states['molecules'][self.unbound_rnap_key]
        monomer_limits = {
            monomer: states['molecules'][monomer]
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

        while time < timestep:
            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                copy_numbers - bound_rnap,
                bound_rnap,
                [unbound_rnaps]])

            # find number of monomers until next terminator
            distance = chromosome.terminator_distance()

            # find interval of time that elongates to the point of the next terminator
            interval = distance / self.elongation_rate

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(interval, substrate)

            # go through each event in the simulation and update the state
            rnap_bindings = 0
            for now, event in zip(result['time'], result['events']):

                # perform the elongation until the next event
                terminations, monomer_limits, chromosome.rnaps = elongation.elongate(
                    time + now,
                    self.elongation_rate,
                    monomer_limits,
                    chromosome.rnaps)
                unbound_rnaps += terminations

                # RNAP has bound the promoter
                if event < self.promoter_count:
                    promoter_key = chromosome.promoter_order[event]
                    promoter = chromosome.promoters[promoter_key]
                    domains = open_domains[promoter_key]
                    domain = choose_element(domains)

                    bound_rnap[event] += 1
                    bound_domains[promoter_key].add(domain)
                    open_domains[promoter_key].remove(domain)

                    # create a new bound RNAP and add it to the chromosome.
                    new_rnap = chromosome.bind_rnap(promoter_key, domain)
                    promoter_rnaps[promoter_key][domain] = new_rnap

                    rnap_bindings += 1
                    unbound_rnaps -= 1
                # RNAP has begun polymerizing its transcript
                else:
                    promoter_index = event - self.promoter_count
                    promoter_key = chromosome.promoter_order[promoter_index]
                    domains = bound_domains[promoter_key]
                    domain = choose_element(domains)

                    open_domains[promoter_key].add(domain)
                    bound_domains[promoter_key].remove(domain)
                    bound_rnap[promoter_index] -= 1

                    # # find rnap on this domain
                    rnap = promoter_rnaps[promoter_key][domain]
                    rnap.start_transcribing()

                    del promoter_rnaps[promoter_key][domain]

            # now that all events have been accounted for, elongate
            # until the end of this interval.
            terminations, monomer_limits, chromosome.rnaps = elongation.elongate(
                time + interval,
                self.elongation_rate,
                monomer_limits,
                chromosome.rnaps)
            unbound_rnaps += terminations

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        molecules = {
            self.unbound_rnap_key: unbound_rnaps - original_unbound_rnaps}

        molecules.update({
            key: count * -1
            for key, count in elongation.monomers.items()})

        update = {
            'chromosome': chromosome.to_dict(),
            'molecules': molecules,
            'transcripts': elongation.complete_polymers}

        return update


def test_transcription():
    parameters = {
        'promoter_affinities': {
            'pA': 1.0,
            'pB': 1.0},
        'elongation_rate': 10.0,
        'advancement_rate': 10.0}

    chromosome = Chromosome(test_chromosome_config)
    transcription = Transcription(parameters)

    states = {
        'chromosome': chromosome.to_dict(),
        'molecules': {transcription.unbound_rnap_key: 10}}

    states['molecules'].update({
        nucleotide: 10
        for nucleotide in transcription.monomer_ids})

    update = transcription.next_update(1.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_transcription()
