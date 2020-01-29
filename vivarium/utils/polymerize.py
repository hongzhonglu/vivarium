from vivarium.utils.datum import Datum

class Polymerize(object):
    def __init__(self, elongation=0, limits={}):
        self.time = 0
        self.monomers = ''
        self.complete_polymers = []
        self.previous_elongations = int(elongation)
        self.elongation = elongation
        self.limits = limits

    def elongate(self, now, rate, limits):
        '''
        Track increments of time and accumulate partial elongations, emitting the full elongation
        once a unit is attained.

        Returns number of polymerases that terminated this step.
        '''

        progress = rate * (now - self.time)
        self.elongation += progress
        elongations = int(self.elongation) - self.previous_elongations
        self.time = now
        terminated = 0

        if elongations:
            iterations, monomers, complete, limits = self.next_polymerize(
                elongations, limits)
            self.monomers += monomers
            self.complete_polymers.extend(complete)
            self.previous_elongations = int(self.elongation)
            terminated += len(complete)

            print('iterations: {}, monomers: {}, complete: {}'.format(
                iterations, monomers, complete))
            print('terminated: {}'.format(terminated))

        return terminated, limits

    def next_polymerize(self, polymerases, elongation_limit=INFINITY, monomer_limits={}):
        elongate_to = min(elongation_limit, distance)
        complete_polymers = []
        monomers = ''

        for step in range(elongate_to):
            for polymerase in polymerases:
                if polymerase.is_transcribing():
                    promoter = self.promoters[polymerase.promoter]
                    extent = promoter.direction
                    projection = polymerase.position + extent

                    monomer = self.sequence[projection]
                    if monomer_limits[monomer] > 0:
                        monomer_limits[monomer] -= 1
                        monomers += monomer
                        polymerase.position = projection

                        terminator = promoter.terminators[polymerase.terminator]
                        if terminator.position == polymerase.position:
                            if promoter.terminates_at(polymerase.terminator):
                                polymerase.complete()
                                complete_transcripts.append(terminator.operon)

        polymerases = [
            polymerase
            for polymerase in self.polymerases
            if not polymerase.is_complete()]

        return elongate_to, monomers, complete_transcripts, monomer_limits, polymerases

    def complete(self):
        return len(self.complete_polymers)

