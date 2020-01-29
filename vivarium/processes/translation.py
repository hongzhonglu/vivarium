import copy

from vivarium.actor.process import Process
from vivarium.data.amino_acids import amino_acid_records
from vivarium.utils.datum import Datum
from vivarium.utils.polymerize import Polymerase, Template

class Ribosome(Polymerase):
    pass

class Transcript(Template):
    pass

class Translation(Process):
    def __init__(self, initial_parameters={}):
        self.monomer_ids = [record['symbol'] for record in amino_acid_records]
        self.default_parameters = {
            'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
            'elongation_rate': 1.0,
            'advancement_rate': 1.0,
            'monomer_ids': self.monomer_ids}
        self.default_parameters['transcript_order'] = list(
            self.default_parameters['transcript_affinities'].keys())

        parameters = copy.deepcopy(self.default_parameters)
        parameters.update(initial_parameters)

        self.transcript_affinities = parameters['transcript_affinities']
        self.transcript_order = parameters['transcript_order']
        self.transcript_count = len(self.transcript_order)

        self.unbound_ribosomes_key = 'unbound_ribosomes'
        self.molecule_ids = parameters['molecule_ids']
        self.monomer_ids = parameters['monomer_ids']
        self.transcript_ids = parameters['transcript_ids']
        self.elongation = 0
        self.elongation_rate = parameters['elongation_rate']
        self.advancement_rate = parameters['advancement_rate']

        self.affinity_vector = np.array([
            self.promoter_affinities[promoter_key]
            for promoter_key in self.promoter_order], dtype=np.float64)

        self.stoichiometry = build_stoichiometry(self.transcript_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        print('stoichiometry: {}'.format(self.stoichiometry))
        print('rates: {}'.format(self.rates))
        self.initiation = StochasticSystem(self.stoichiometry, self.rates)

        self.ribosome_id = 0

        self.roles = {
            'ribosomes': ['ribosomes'],
            'molecules': self.molecule_ids,
            'transcripts': self.transcript_ids}

        super(Transcription, self).__init__(self.roles, parameters)

    def default_settings(self):
        default_state = {
            'ribosomes': {
                'ribosomes': []},
            'molecules': dict({
                self.unbound_ribosomes_key: 10}),
            'transcripts': {
                'oA': 10,
                'oAZ': 15,
                'oB': 25,
                'oBY': 40}}

        default_state['molecules'].update({
            monomer_id: 100
            for monomer_id in self.monomer_ids})
        chromosome = Chromosome(default_state['chromosome'])
        operons = [operon.id for operon in chromosome.operons()]

        default_emitter_keys = {
            'ribosomes': ['ribosomes'],
            'molecules': self.monomer_ids + [self.unbound_ribosomes_key],
            'transcripts': operons}

        default_updaters = {
            'ribosomes': {},
            'molecules': {},
            'transcripts': {}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'parameters': self.default_parameters}

    def next_update(self, timestep, states):
        ribosomes = states['ribosomes']
        molecules = states['molecules']
        transcripts = states['transcripts']
        transcript_counts = np.array([
            transcripts.get(transcript_key, 0)
            for transcript_key in self.transcript_order])

        # Find out how many transcripts are currently blocked by a
        # newly initiated ribosome
        bound_transcripts = np.zeros(len(transcript_order))
        ribosomes_by_transcript = {}
        for ribosome in ribosomes:
            if not ribosome.template in ribosomes_by_transcript:
                ribosomes_by_transcript[ribosome.template] = []
            ribosomes_by_transcript[ribosome.template].append(ribosome)
        for index, transcript in enumerate(self.transcript_order):
            if transcript in ribosomes_by_transcript:
                bound_transcripts[index] = len(ribosomes_by_transcript[transcript])

        # for promoter_key in chromosome.promoter_order:
        #     bound_domains[promoter_key] = set([
        #         rnap.domain
        #         for rnap in promoter_rnaps.get(promoter_key, {}).values()
        #         if rnap.is_bound()])
        #     bound_rnap.append(len(bound_domains[promoter_key]))
        #     open_domains[promoter_key] = promoter_domains[promoter_key] - bound_domains[promoter_key]

        # bound_rnap = np.array(bound_rnap)

        # Make the state for a gillespie simulation out of total number of each
        # promoter by copy number not blocked by initiated rnap,
        # concatenated with the number of each promoter that is bound by rnap.
        # These are the two states for each promoter the simulation
        # will operate on, essentially going back and forth between
        # bound and unbound states.

        original_unbound_ribosomes = states['molecules']['unbound_ribosomes']
        monomer_limits = {
            monomer: states['molecules'][monomer]
            for monomer in self.monomer_ids}
        unbound_ribosomes = original_unbound_ribosomes

        time = 0
        now = 0
        elongation = Elongation(
            sequences,
            templates,
            monomer_limits,
            self.elongation)

        while time < timestep:
            print('time: {} -----------======--------------------------'.format(time))

            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                transcript_counts - bound_transcripts,
                bound_transcripts,
                [unbound_ribosomes]])

            print('state: {}'.format(substrate))
            print('unbound rnaps: {}'.format(unbound_rnaps))

            # find number of monomers until next terminator
            # distance = chromosome.terminator_distance()
            distance = 1

            print('distance: {}'.format(distance))

            # find interval of time that elongates to the point of the next terminator
            interval = distance / self.elongation_rate

            print('interval: {}'.format(interval))
            print('substrates: {}'.format(substrate))

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(interval, substrate)

            # go through each event in the simulation and update the state
            ribosome_bindings = 0
            for now, event in zip(result['time'], result['events']):

                print('event {}: {}'.format(now, event))

                # perform the elongation until the next event
                terminations, monomer_limits, ribosomes = elongation.elongate(
                    time + now,
                    self.elongation_rate,
                    monomer_limits,
                    ribosomes)
                unbound_ribosomes += terminations

                # RNAP has bound the promoter
                if event < self.transcript_count:
                    transcript_key = self.transcript_order[event]
                    transcript = self.templates[transcript_key]
                    bound_ribosomes[event] += 1

                    self.ribosome_id += 1
                    new_ribosome = Ribosome({
                        'id': self.ribosome_id,
                        'template': promoter_key,
                        'position': 0})
                    new_ribosome.bind()
                    ribosomes_by_transcript[transcript_key] = new_ribosome

                    print('{}: ribosome binding {}'.format(time, new_ribosome))

                    ribosome_bindings += 1
                    unbound_ribosomes -= 1
                # RNAP has begun polymerizing its transcript
                else:
                    transcript_index = event - self.transcript_count
                    transcript_key = self.transcript_order[transcript_index]

                    bound_ribosomes[transcript_index] -= 1

                    # # find rnap on this domain
                    ribosome = ribosomes_by_transcript[transcript_key].pop()
                    ribosome.start_transcribing()

                    print('{}: ribosome commencing {}'.format(time, ribosome))

            # now that all events have been accounted for, elongate
            # until the end of this interval.
            terminations, monomer_limits, ribosomes = elongation.elongate(
                time + interval,
                self.elongation_rate,
                monomer_limits,
                ribosomes)
            unbound_ribosomes += terminations

            print('bound ribosomes: {}'.format(ribosomes))
            print('complete transcripts: {}'.format(elongation.complete_polymers))
            print('monomer limits: {}'.format(monomer_limits))

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        molecules = {
            'unbound_ribosomes': unbound_ribosomes - original_unbound_ribosomes}

        molecules.update({
            key: count * -1
            for key, count in elongation.monomers.items()})

        update = {
            'ribosomes': {
                'ribosomes': [ribosome.to_dict() for ribosome in ribosomes]},
            'molecules': molecules,
            'transcripts': elongation.complete_polymers}

        print('molecules update: {}'.format(molecules))

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

    translation = Translation(parameters)

    states = {
        'ribosomes': [],
        'molecules': {translation.unbound_ribosomes_key: 10}}
    states['molecules'].update({
        molecule_id: 100
        for molecule_id in translation.molecule_ids})

    update = transcription.next_update(1.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_translation()
        
