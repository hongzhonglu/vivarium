import numpy as np
from arrow import StochasticSystem

from vivarium.actor.process import Process
from vivarium.states.chromosome import Chromosome, Rnap, test_chromosome

def build_stoichiometry(promoter_count):
    stoichiometry = np.zeros((promoter_count * 2, promoter_count * 2), dtype=np.int64)
    for index in range(promoter_count):
        # forward reaction
        stoichiometry[index][index] = -1
        stoichiometry[index][index + promoter_count] = 1

        # reverse reaction
        stoichiometry[index + promoter_count][index] = 1
        stoichiometry[index + promoter_count][index + promoter_count] = -1

    return stoichiometry

def build_rates(affinities, advancement):
    return np.concatenate([
        affinities,
        np.repeat(advancement, len(affinities))])

class Transcription(Process):
    def __init__(self, initial_parameters={}):
        self.promoter_affinities = initial_parameters.get('promoter_affinities', {})
        self.promoter_order = initial_parameters.get('promoter_order', [])
        self.promoter_count = len(self.promoter_affinities)
        self.affinity_vector = np.array([
            self.promoter_affinities[promoter_key]
            for promoter_key in self.promoter_order])

        self.elongation_rate = initial_parameters.get('elongation_rate', 50.0)
        self.advancement_rate = initial_parameters.get('advancement_rate', 1.0)

        self.stoichiometry = build_stoichiometry(self.promoter_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        self.initiation = StochasticSystem(self.stoichiometry, self.rates)

    def next_update(self, timestep, states):
        chromosome = Chromosome(states['chromosome'])
        molecules = states['molecules']

        # Find out how many promoters are currently blocked by a
        # newly initiated rnap
        promoter_rnaps = chromosome.promoter_rnaps()
        promoter_domains = chromosome.promoter_domains()

        bound_rnap = []
        open_domains = {}
        bound_domains = {}
        for promoter_key in chromosome.promoter_order:
            bound_domains[promoter_key] = set([
                rnap.domain
                for rnap in promoter_rnaps.get(promoter_key, [])
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
        substrate = np.concatenate([
            copy_numbers - bound_rnap,
            bound_rnap])

        # run simulation
        result = self.initiation.evolve(timestep, substrate)

        def choose_element(elements):
            if elements:
                choice = np.random.choice(len(elements), 1)
                return list(elements)[int(choice)]

        rnap_bindings = 0
        complete_transcripts = []
        time = 0
        for elapse, event in zip(result['time'], result['events']):
            complete = chromosome.polymerize(self.elongation_rate * (elapse - time))
            complete_transcripts.extend(complete)
            time = elapse

            if event < self.promoter_count:
                # RNAP has bound the promoter
                promoter_key = chromosome.promoter_order[event]
                promoter = chromosome.promoters[promoter_key]
                domains = open_domains[promoter_key]
                domain = choose_element(domains)

                bound_domains[promoter_key].add(domain)
                open_domains[promoter_key].remove(domain)

                # create a new bound RNAP and add it to the chromosome.
                new_rnap = chromosome.bind_rnap(promoter_key, domain)
                promoter_rnaps[promoter_key].append(new_rnap)

                print('{}: RNAP binding {}'.format(time, new_rnap))

                rnap_bindings += 1
            else:
                # RNAP has begun polymerizing its transcript
                promoter_index = event - self.promoter_count
                promoter_key = chromosome.promoter_order[promoter_index]
                domains = bound_domains[promoter_key]
                domain = choose_element(domains)

                open_domains[promoter_key].add(domain)
                bound_domains[promoter_key].remove(domain)

                # find rnap on this domain
                # TODO(Ryan): don't search through a list for this
                for rnap in promoter_rnaps[promoter_key]:
                    if rnap.domain == domain:
                        break

                rnap.start_transcribing()

                print('{}: RNAP commencing {}'.format(time, rnap))

        remaining = timestep - time
        complete = chromosome.polymerize(self.elongation_rate * remaining)
        complete_transcripts.extend(complete)

        return chromosome.to_dict()


def test_transcription():
    parameters = {
        'promoter_affinities': {
            'pA': 2,
            'pB': 5},
        'promoter_order': ['pA', 'pB']}

    chromosome = test_chromosome()
    transcription = Transcription(parameters)

    states = {
        'chromosome': chromosome.to_dict(),
        'molecules': {}}

    update = transcription.next_update(1.0, states)
    
    print(update)



if __name__ == '__main__':
    test_transcription()
