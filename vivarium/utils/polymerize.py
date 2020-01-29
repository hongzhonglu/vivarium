import random

from vivarium.utils.datum import Datum

INFINITY = float('inf')

def add_merge(ds):
    '''
    Given a list of dicts, sum the values of each key.
    '''

    result = {}
    for d in ds:
        for key, value in d.items():
            if not key in result:
                result[key] = 0
            result[key] += value
    return result

class Polymerase(Datum):
    defaults = {
        'id': 0,
        'state': None, # other states: ['bound', 'transcribing', 'complete']
        'position': 0,
        'template': '',
        'terminator': 0}

    def __init__(self, config, defaults=defaults):
        super(Polymerase, self).__init__(config, self.defaults)

    def bind(self):
        self.state = 'bound'

    def start_transcribing(self):
        self.state = 'transcribing'

    def complete(self):
        self.state = 'complete'
        print('completing polymerization: {}'.format(self.to_dict()))

    def is_bound(self):
        return self.state == 'bound'

    def is_transcribing(self):
        return self.state == 'transcribing'

    def is_complete(self):
        return self.state == 'complete'

class BindingSite(Datum):
    defaults = {
        'position': 0,
        'length': 0,
        'thresholds': []} # list of pairs, (factor, threshold)

    def __init__(self, config):
        super(BindingSite, self).__init__(config, self.defaults)

    def state_when(self, levels):
        '''
        Provide the binding state for the given levels of factors. 
        '''

        state = None
        for factor, threshold in thresholds:
            if levels[factor] >= threshold:
                state = factor
                break
        return state

class Terminator(Datum):
    defaults = {
        'position': 0,
        'strength': 0,
        'product': ''}

    def __init__(self, config, defaults=defaults):
        super(Terminator, self).__init__(config, self.defaults)

    def between(self, before, after):
        return before < self.position < after or after < self.position < before

class Template(Datum):
    schema = {
        'sites': BindingSite,
        'terminators': Terminator}

    defaults = {
        'id': 0,
        'position': 0,
        'direction': 1,
        'sites': [],
        'terminators': []}

    def __init__(self, config, defaults=defaults):
        super(Template, self).__init__(config, self.defaults)

        self.terminator_strength = 0
        for terminator in self.terminators:
            self.terminator_strength += terminator.strength

    def binding_state(self, levels):
        state = [
            site.state_when(levels)
            for site in self.sites]

        return tuple([self.id] + state)

    def strength_from(self, terminator_index):
        total = 0
        for index in range(terminator_index, len(self.terminators)):
            total += self.terminators[index].strength
        return total

    def next_terminator(self, position):
        for index, terminator in enumerate(self.terminators):
            if terminator.position * self.direction > position * self.direction:
                break
        return index

    def terminates_at(self, index=0):
        if len(self.terminators[index:]) > 1:
            choice = random.random() * self.strength_from(index)
            return choice <= self.terminators[index].strength
        else:
            return True

    def choose_terminator(self, index=0):
        if len(self.terminators[index:]) > 1:
            choice = random.random() * self.strength_from(index)
            for terminator in self.terminators[index:]:
                if choice <= terminator.strength:
                    break
                else:
                    choice -= terminator.strength
            return terminator
        else:
            return self.terminators[index]

    def choose_product(self):
        terminator = self.choose_terminator()
        return terminator.product

    def products(self):
        return [
            terminator.product
            for terminator in self.terminators]


def polymerize_to(
        sequence,
        polymerases,
        templates,
        additions,
        monomer_limits={}):

    complete_polymers = {}
    monomers = {
        monomer: 0
        for monomer in monomer_limits.keys()}

    for step in range(additions):
        for polymerase in polymerases:
            if polymerase.is_transcribing():
                template = templates[polymerase.template]
                extent = template.direction
                projection = polymerase.position + extent

                monomer = sequence[projection]
                if monomer_limits[monomer] > 0:
                    monomer_limits[monomer] -= 1
                    monomers[monomer] += 1
                    polymerase.position = projection

                    terminator = template.terminators[polymerase.terminator]
                    if terminator.position == polymerase.position:
                        if template.terminates_at(polymerase.terminator):
                            polymerase.complete()

                            if not terminator.product in complete_polymers:
                                complete_polymers[terminator.product] = 0
                            complete_polymers[terminator.product] += 1

    polymerases = [
        polymerase
        for polymerase in polymerases
        if not polymerase.is_complete()]

    return monomers, monomer_limits, complete_polymers, polymerases


class Elongation(object):
    def __init__(self, sequence, templates, elongation=0, limits={}):
        self.sequence = sequence
        self.templates = templates
        self.time = 0
        self.monomers = {}
        self.complete_polymers = {}
        self.previous_elongations = int(elongation)
        self.elongation = elongation
        self.limits = limits

    def elongate(self, now, rate, limits, polymerases):
        '''
        Track increments of time and accumulate partial elongations, emitting the full
        elongation once a unit is attained.

        Returns number of polymerases that terminated this step, and the updated 
        monomer limits after all elongations.
        '''

        progress = rate * (now - self.time)
        self.elongation += progress
        elongations = int(self.elongation) - self.previous_elongations
        self.time = now
        terminated = 0

        if elongations:
            monomers, limits, complete, polymerases = polymerize_to(
                self.sequence,
                polymerases,
                self.templates,
                elongations,
                limits)
            self.monomers = add_merge([self.monomers, monomers])
            self.complete_polymers = add_merge([
                self.complete_polymers, complete])
            self.previous_elongations = int(self.elongation)
            terminated += len(complete)

        return terminated, limits, polymerases

    def complete(self):
        return len(self.complete_polymers)

