import math
import numpy as np
from arrow import StochasticSystem

from vivarium.actor.process import Process
from vivarium.states.chromosome import Chromosome, Rnap, test_chromosome, frequencies

def build_stoichiometry(promoter_count):
    '''
    Builds a stoichiometry for the given promoters. There are two states per promoter,
    open and bound, and two reactions per promoter, binding and unbinding. In addition
    there is a single substrate for available RNAP in the final index.

    Here we are assuming
    '''
    stoichiometry = np.zeros((promoter_count * 2, promoter_count * 2 + 1), dtype=np.int64)
    for index in range(promoter_count):
        # forward reaction
        stoichiometry[index][index] = -1
        stoichiometry[index][index + promoter_count] = 1
        stoichiometry[index][-1] = -1 # forward reaction consumes RNAP also

        # reverse reaction
        stoichiometry[index + promoter_count][index] = 1
        stoichiometry[index + promoter_count][index + promoter_count] = -1

    return stoichiometry

def build_rates(affinities, advancement):
    return np.concatenate([
        affinities,
        np.repeat(advancement, len(affinities))])

def choose_element(elements):
    if elements:
        choice = np.random.choice(len(elements), 1)
        return list(elements)[int(choice)]

class Elongation(object):
    def __init__(self, elongation=0):
        self.time = 0
        self.monomers = ''
        self.complete_transcripts = []
        self.previous_elongations = elongation
        self.elongation = elongation

    def elongate(self, chromosome, now, rate):
        '''
        Track increments of time and accumulate partial elongations, emitting the full elongation
        once a unit is attained.

        Returns number of RNAP that terminated transcription this step.
        '''

        progress = rate * (now - self.time)
        self.elongation += progress
        elongations = int(self.elongation) - self.previous_elongations
        self.time = now
        terminated = 0

        if elongations:
            iterations, monomers, complete = chromosome.next_polymerize(elongations)
            self.monomers += monomers
            self.complete_transcripts.extend(complete)
            self.previous_elongations = int(self.elongation)
            terminated += len(complete)

        return terminated

    def complete(self):
        return len(self.complete_transcripts)

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

        self.promoter_affinities = initial_parameters.get('promoter_affinities', {})
        self.promoter_order = initial_parameters.get('promoter_order', [])
        self.promoter_count = len(self.promoter_order)
        self.affinity_vector = np.array([
            self.promoter_affinities[promoter_key]
            for promoter_key in self.promoter_order])

        self.elongation = 0
        self.elongation_rate = initial_parameters.get('elongation_rate', 1.0)
        self.advancement_rate = initial_parameters.get('advancement_rate', 1.0)

        self.stoichiometry = build_stoichiometry(self.promoter_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        self.initiation = StochasticSystem(self.stoichiometry, self.rates)


    def next_update(self, timestep, states):
        chromosome = Chromosome(states['chromosome'])
        molecules = states['molecules']

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
        original_unbound_rnaps = states['molecules']['unbound_rnaps']
        unbound_rnaps = original_unbound_rnaps

        time = 0
        elongation = Elongation(self.elongation)

        while time < timestep:
            print('time: {} --------------------------------------------------------'.format(time))

            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                copy_numbers - bound_rnap,
                bound_rnap,
                [unbound_rnaps]])

            print('state: {}'.format(substrate))
            print('unbound rnaps: {}'.format(unbound_rnaps))

            # find number of monomers until next terminator
            distance = chromosome.terminator_distance()

            print('distance: {}'.format(distance))

            # find interval of time that elongates to the point of the next terminator
            interval = distance / self.elongation_rate

            print('interval: {}'.format(interval))

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(interval, substrate)

            # go through each event in the simulation and update the state
            rnap_bindings = 0
            for now, event in zip(result['time'], result['events']):

                print('event {}: {}'.format(now, event))

                # perform the elongation until the next event
                unbound_rnaps += elongation.elongate(
                    chromosome,
                    time + now,
                    self.elongation_rate)

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

                    print('{}: RNAP binding {}'.format(time, new_rnap))

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

                    print('{}: RNAP commencing {}'.format(time, rnap))

            # now that all events have been accounted for, elongate
            # until the end of this interval.
            elongation.elongate(
                chromosome,
                time + interval,
                self.elongation_rate)

            print('bound rnaps: {}'.format(chromosome.rnaps))
            print('complete_transcripts: {}'.format(elongation.complete_transcripts))

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        molecules = {
            'unbound_rnaps': unbound_rnaps - original_unbound_rnaps}

        molecules.update({
            key: count * -1
            for key, count in frequencies(elongation.monomers).items()})

        transcripts = frequencies(elongation.complete_transcripts)

        return {
            'chromosome': chromosome.to_dict(),
            'molecules': molecules,
            'transcripts': transcripts}


def test_transcription():
    parameters = {
        'promoter_affinities': {
            'pA': 1.0,
            'pB': 1.0},
        'promoter_order': ['pA', 'pB'],
        'elongation_rate': 10.0,
        'advancement_rate': 10.0}

    chromosome = test_chromosome()
    transcription = Transcription(parameters)

    states = {
        'chromosome': chromosome.to_dict(),
        'molecules': {
            'unbound_rnaps': 10}}

    update = transcription.next_update(1.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_transcription()
