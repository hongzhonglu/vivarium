import random
import copy
import numpy as np

from vivarium.data.chromosome import test_chromosome_config

INFINITY = float('inf')

def first(l):
    if l:
        return l[0]

def first_value(d):
    if d:
        return d[list(d.keys())[0]]

def flatten(l):
    '''
    Flatten a list by one level:
        [[1, 2, 3], [[4, 5], 6], [7]] --> [1, 2, 3, [4, 5], 6, 7]
    '''

    return [
        item
        for sublist in l
        for item in sublist]

def frequencies(l):
    '''
    Return number of times each item appears in the list.
    '''

    result = {}
    for item in l:
        if not item in result:
            result[item] = 0
        result[item] += 1
    return result

def traverse(tree, key, f, combine):
    '''
    Traverse the given tree starting using the `key` node as the root and calling `f` on each leaf,
    combining values with `combine` at each subsequent level to create new leaves for `f`.
    '''

    node = tree[key]
    if node.children:
        eldest = traverse(tree, node.children[0], f, combine)
        youngest = traverse(tree, node.children[1], f, combine)
        outcome = combine(eldest, youngest)
        return f(node, outcome)
    else:
        return f(node)


class Datum(object):
    '''
    The Datum class enables functions to be defined on dicts of a certain schema. 
    It provides two class level variables:
      * `defaults`: a dictionary of keys to default values this Datum will have if 
           none is provided to __init__
      * `schema`: a dictionary of keys to constructors which invoke subdata. 

    Once these are defined, a Datum subclass can be constructed with a dict that provides any
    values beyond the defaults, and then all of the defined methods for that Datum subclass
    are available to operate on its values. Once the modifications are complete, it can be
    rendered back into a dict using the `to_dict()` method.
    '''

    schema = {}
    defaults = {}

    def __init__(self, config, default):
        self.keys = list(set(list(config.keys()) + list(default.keys()))) # a dance
        for key in self.keys:
            value = config.get(key, default[key])
            if value and key in self.schema:
                realize = self.schema[key]
                if isinstance(value, list):
                    value = [realize(item) for item in value]
                elif isinstance(value, dict):
                    value = {inner: realize(item) for inner, item in value.items()}
                else:
                    value = realize(item)
            setattr(self, key, value)

    def fields(self):
        return list(self.defaults.keys())

    def to_dict(self):
        to = {}
        for key in self.keys:
            value = getattr(self, key)
            if isinstance(value, Datum):
                value = value.to_dict()
            elif value and isinstance(value, list) and isinstance(first(value), Datum):
                value = [datum.to_dict() for datum in value]
            elif value and isinstance(value, dict) and isinstance(first_value(value), Datum):
                value = {inner: datum.to_dict() for inner, datum in value.items()}
            to[key] = value
        return to

    def __repr__(self):
        return str(self.to_dict())

class Operon(Datum):
    defaults = {
        'id': '',
        'position': 0,
        'direction': 1,
        'length': 0,
        'genes': []}

    def __init__(self, config):
        super(Operon, self).__init__(config, self.defaults)

class Domain(Datum):
    defaults = {
        'id': 0,
        'lead': 0,
        'lag': 0,
        'children': []}

    def __init__(self, config):
        super(Domain, self).__init__(config, self.defaults)

    def contains(self, position):
        if position < 0:
            return position < self.lag
        else:
            return position > self.lead

    def surpassed(self, position, lead, lag):
        if position < 0:
            return position <= self.lag and position > self.lag + lag
        else:
            return position >= self.lead and position < self.lead + lead

    def random_child(self):
        return random.choice(self.children)

    def descendants(self, tree):
        return [self] + [tree[child].descendants(tree) for child in self.children]

class BindingSite(Datum):
    defaults = {
        'position': 0,
        'length': 0,
        'thresholds': []} # list of pairs, (TF, threshold)

    def __init__(self, config):
        super(BindingSite, self).__init__(config, self.defaults)

    def state_when(self, levels):
        '''
        Provide the binding state for the given levels of transcription factors. 
        '''

        state = None
        for tf, threshold in thresholds:
            if levels[tf] >= threshold:
                state = tf
                break
        return state

class Terminator(Datum):
    defaults = {
        'position': 0,
        'strength': 0,
        'operon': ''}

    def __init__(self, config):
        super(Terminator, self).__init__(config, self.defaults)

    def operon_from(self, genes, promoter):
        return Operon({
            'id': self.operon,
            'position': promoter.position,
            'direction': promoter.direction,
            'length': (self.position - promoter.position) * promoter.direction,
            'genes': genes.get(self.operon, [])})

    def between(self, before, after):
        return before < self.position < after or after < self.position < before

class Promoter(Datum):
    '''
    Promoters are the main unit of expression. They define a direction of polymerization, 
    contain binding sites for transcription factors, and declare some number of terminators,
    each of which has a strength and corresponds to a particular operon if chosen.
    '''

    schema = {
        'sites': BindingSite,
        'terminators': Terminator}

    defaults = {
        'id': 0,
        'position': 0,
        'direction': 1,
        'sites': [],
        'terminators': []}

    def __init__(self, config):
        super(Promoter, self).__init__(config, self.defaults)

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

    def choose_operon(self, genes):
        terminator = self.choose_terminator()
        return terminator.operon_from(genes, self)

    def operons(self, genes):
        return [
            terminator.operon_from(genes, self)
            for terminator in self.terminators]

class Rnap(Datum):
    defaults = {
        'id': 0,
        'promoter': '',
        'terminator': 0,
        'domain': 0,
        'state': None, # other states: ['bound', 'transcribing', 'complete']
        'position': 0}

    def __init__(self, config):
        super(Rnap, self).__init__(config, self.defaults)

    def bind(self):
        self.state = 'bound'

    def start_transcribing(self):
        self.state = 'transcribing'

    def complete(self):
        self.state = 'complete'
        print('RNAP completing transcription: {}'.format(self.to_dict()))

    def is_bound(self):
        return self.state == 'bound'

    def is_transcribing(self):
        return self.state == 'transcribing'

    def is_complete(self):
        return self.state == 'complete'

class Chromosome(Datum):
    schema = {
        'promoters': Promoter,
        'domains': Domain,
        'rnaps': Rnap}

    defaults = {
        'sequence': '',
        'genes': {},
        'promoters': {},
        'domains': {},
        'root_domain': 0,
        'rnaps': []}

    def operons(self):
        return [
            operon
            for promoter in self.promoters.values()
            for operon in promoter.operons(self.genes)]

    def copy_number(self, position, domain_key=None):
        if not domain_key:
            domain_key = self.root_domain
        domain = self.domains[domain_key]
        if domain.contains(position):
            return 1
        else:
            return sum([
                self.copy_number(position, child)
                for child in domain.children])

    def promoter_copy_numbers(self):
        copy_numbers = [
            self.copy_number(self.promoters[promoter_key].position)
            for promoter_key in self.promoter_order]
        return np.array(copy_numbers)

    def promoter_rnaps(self):
        by_promoter = {
            promoter_key: {}
            for promoter_key in self.promoter_order}

        for rnap in self.rnaps:
            if rnap.is_bound():
                by_promoter[rnap.promoter][rnap.domain] = rnap

        return by_promoter

    def promoter_domains(self):
        return {
            promoter_key: self.position_domains(
                self.root_domain,
                self.promoters[promoter_key].position)
            for promoter_key in self.promoter_order}

    def position_domains(self, domain_index, position):
        domain = self.domains[domain_index]
        if len(domain.children) == 0 or (position < 0 and domain.lag >= position) or (position >= 0 and domain.lead <= position):
            return set([domain_index])
        else:
            return set.union(*[
                self.position_domains(child, position)
                for child in domain.children])

    def bind_rnap(self, promoter_key, domain):
        self.rnap_id += 1
        new_rnap = Rnap({
            'id': self.rnap_id,
            'promoter': promoter_key,
            'domain': domain,
            'position': self.promoters[promoter_key].position})
        new_rnap.bind()
        self.rnaps.append(new_rnap)
        return new_rnap

    def terminator_distance(self):
        distance = INFINITY
        for rnap in self.rnaps:
            if rnap.is_transcribing():
                promoter = self.promoters[rnap.promoter]
                terminator_index = promoter.next_terminator(rnap.position)
                rnap.terminator = terminator_index
                terminator = promoter.terminators[terminator_index]
                span = abs(terminator.position - rnap.position)
                if span < distance:
                    distance = span
        if distance == INFINITY:
            distance = 1
        return distance

    def sequence_monomers(self, begin, end):
        if begin < end:
            return self.sequence[begin:end]
        else:
            return self.sequence[end:begin]

    def next_polymerize(self, elongation_limit=INFINITY, monomer_limits={}):
        distance = self.terminator_distance()
        elongate_to = min(elongation_limit, distance)

        complete_transcripts = []
        monomers = ''

        for step in range(elongate_to):
            for rnap in self.rnaps:
                if rnap.is_transcribing():
                    promoter = self.promoters[rnap.promoter]
                    extent = promoter.direction
                    projection = rnap.position + extent

                    monomer = self.sequence[projection]
                    if monomer_limits[monomer] > 0:
                        monomer_limits[monomer] -= 1
                        monomers += monomer
                        rnap.position = projection

                        terminator = promoter.terminators[rnap.terminator]
                        if terminator.position == rnap.position:
                            if promoter.terminates_at(rnap.terminator):
                                rnap.complete()
                                complete_transcripts.append(terminator.operon)

        self.rnaps = [
            rnap
            for rnap in self.rnaps
            if not rnap.is_complete()]

        return elongate_to, monomers, complete_transcripts, monomer_limits

    def polymerize(self, elongation, monomer_limits):
        iterations = 0
        attained = 0
        all_monomers = ''
        complete_transcripts = []

        while attained < elongation:
            elongated, monomers, complete, monomer_limits = self.next_polymerize(
                elongation_limit=elongation - attained,
                monomer_limits=monomer_limits)
            attained += elongated
            all_monomers += monomers
            complete_transcripts += complete
            iterations += 1

        return iterations, frequencies(all_monomers), complete_transcripts, monomer_limits

    def initiate_replication(self):
        leaves = [leaf for leaf in self.domains.values() if not leaf.children]
        next_id = max([leaf.id for leaf in leaves]) + 1
        for leaf in leaves:
            for child in [0, 1]:
                domain = Domain({'id': next_id + child})
                self.domains[domain.id] = domain
            leaf.children = [next_id, next_id + 1]
            next_id += 2

    def advance_replisomes(self, distances):
        '''
        distances is a dictionary of domain ids to tuples of how far each strand advances
        of the form (lead, lag)
        '''
        for domain_key, distance in distances.items():
            domain = self.domains[domain_key]
            lead, lag = distances[domain_key]

            for rnap in self.rnaps:
                if rnap.domain == domain_key:
                    promoter = self.promoters[rnap.promoter]
                    position = promoter.position
                    position += rnap.position * promoter.direction
                    if domain.surpassed(position, lead, -lag):
                        rnap.domain = domain.random_child()

            domain.lead += lead
            domain.lag -= lag

    def divide_chromosome(self, domain, division=None):
        if not division:
            division = {
                'sequence': self.sequence,
                'promoters': {id: promoter.to_dict() for id, promoter in self.promoters.items()},
                'domains': {domain.id: domain.to_dict()},
                'root_domain': domain.id,
                'rnaps': [
                    rnap.to_dict()
                    for rnap in self.rnaps
                    if rnap.domain == domain.id]}

        else:
            division['domains'][domain.id] = domain.to_dict()
            for rnap in self.rnaps:
                if rnap.domain == domain.id:
                    division['rnaps'].append(rnap.to_dict())

        return division

    def combine_state(self, a, b):
        merge = copy.deepcopy(a)
        merge['domains'].update(b['domains'])
        merge['rnaps'].extend(b['rnaps'])
        return merge

    def terminate_replication(self):
        children = self.domains[self.root_domain].children
        divided = [
            traverse(
                self.domains,
                child,
                self.divide_chromosome,
                self.combine_state)
            for child in children]

        return [Chromosome(fork) for fork in divided]

    def __init__(self, config):
        super(Chromosome, self).__init__(config, self.defaults)
        self.promoter_order = list(self.promoters.keys())
        self.rnap_id = 0



def test_chromosome():
    chromosome = Chromosome(test_chromosome_config)
    print(chromosome.promoters['pA'].terminators[0].operon)
    print(chromosome)

    print('operons:')
    print(chromosome.operons())

    chromosome.initiate_replication()
    print(chromosome.domains)
    assert len(chromosome.domains) == 3

    chromosome.advance_replisomes({0: (5, 7)})
    print('replisomes:')
    print(chromosome)

    # chromosome.initiate_replication()
    # print(chromosome.domains)
    # assert len(chromosome.domains) == 7

    # chromosome.advance_replisomes({0: (7, 5), 1: (3, 4), 2: (8, 9)})
    # print('replisomes:')
    # print(chromosome)

    operon = chromosome.promoters['pA'].choose_operon(chromosome.genes)
    print('operon')
    print(operon)

    print('copy numbers 1 4 7 -9 11')
    print(chromosome.copy_number(1))
    print(chromosome.copy_number(4))
    print(chromosome.copy_number(7))
    print(chromosome.copy_number(-9))
    print(chromosome.copy_number(11))

    print('promoter copy numbers')
    print(chromosome.promoter_copy_numbers())

    print('promoter rnaps')
    print(chromosome.promoter_rnaps())

    print('promoter domains')
    print(chromosome.promoter_domains())

    print('rnaps')
    print([rnap.to_dict() for rnap in chromosome.rnaps])

    print('completed after advancing 5')
    print(chromosome.polymerize(5, {
        'A': 100, 'T': 100, 'G': 100, 'C': 100}))

    print('rnaps after polymerizing')
    print(chromosome.rnaps)
    
    children = chromosome.terminate_replication()
    print('termination:')
    print(children)

    return chromosome

if __name__ == '__main__':
    test_chromosome()
