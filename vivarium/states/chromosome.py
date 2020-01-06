import random
import copy

def first(l):
    if l:
        return l[0]

def first_value(d):
    if d:
        return d[list(d.keys())[0]]

def traverse(tree, key, f, combine):
    node = tree[key]
    if node.children:
        eldest = traverse(tree, node.children[0], f, combine)
        youngest = traverse(tree, node.children[1], f, combine)
        outcome = combine(eldest, youngest)
        return f(node, outcome)
    else:
        return f(node)


class Datum(object):
    schema = {}

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

class Promoter(Datum):
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

    def choose_terminator(self):
        if len(self.terminators) > 1:
            choice = random.random() * self.terminator_strength
            for terminator in self.terminators:
                if choice <= terminator.strength:
                    break
                else:
                    choice -= terminator.strength
            return terminator
        else:
            return self.terminators[0]

    def choose_operon(self, genes):
        terminator = self.choose_terminator()
        return terminator.operon_from(genes, self)

    def operons(self, genes):
        return [
            terminator.operon_from(genes, self)
            for terminator in self.terminators]

class Rnap(Datum):
    defaults = {
        'promoter': '',
        'domain': 0,
        'position': 0}

    def __init__(self, config):
        super(Rnap, self).__init__(config, self.defaults)

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
        'rnaps': []}

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
        root = min(self.domains.keys())
        children = self.domains[root].children
        divided = [
            traverse(
                self.domains,
                child,
                self.divide_chromosome,
                self.combine_state)
            for child in children]

        return [Chromosome(fork) for fork in divided]

    def operons(self):
        return [
            operon
            for promoter in self.promoters.values()
            for operon in promoter.operons(self.genes)]
                
    def __init__(self, config):
        super(Chromosome, self).__init__(config, self.defaults)


def test_chromosome():
    chromosome_config = {
        'sequence': 'ATACGGCACGTGACCGTCAACTTA',
        'genes': {
            'oAZ': ['A', 'Z'],
            'oA': ['A'],
            'oB': ['B']},
        'promoters': {
            'pA': {
                'id': 'pA',
                'position': 3,
                'direction': 1,
                'sites': [{
                    'position': 0,
                    'length': 3,
                    'thresholds': [
                        ('tfX', 0.3)]}],
                'terminators': [
                    {
                        'position': 6,
                        'strength': 0.5,
                        'operon': 'oA'},
                    {
                        'position': 11,
                        'strength': 1.0,
                        'operon': 'oAZ'}]},
            'pB': {
                'id': 'pB',
                'position': -3,
                'direction': -1,
                'sites': [{
                    'position': 0,
                    'length': 3,
                    'thresholds': [
                        ('tfX', 0.5)]}],
                'terminators': [
                    {
                        'position': -6,
                        'strength': 0.5,
                        'operon': 'oB'},
                    {
                        'position': -11,
                        'strength': 1.0,
                        'operon': 'oBY'}]}},
        'domains': {
            0: {
                'id': 0,
                'lead': 0,
                'lag': 0,
                'children': []}},
        'rnaps': [
            {
                'promoter': 'pA',
                'domain': 0,
                'position': 3},
            {
                'promoter': 'pA',
                'domain': 0,
                'position': 6},
            {
                'promoter': 'pA',
                'domain': 0,
                'position': 0}]}

    chromosome = Chromosome(chromosome_config)
    print(chromosome.promoters['pA'].terminators[0].operon)
    print(chromosome.to_dict())

    assert chromosome.to_dict() == chromosome_config

    print('operons:')
    print([operon.to_dict() for operon in chromosome.operons()])

    chromosome.initiate_replication()
    print(chromosome.to_dict()['domains'])
    assert len(chromosome.domains) == 3

    chromosome.advance_replisomes({0: (5, 7)})
    print('replisomes:')
    print(chromosome.to_dict())

    chromosome.initiate_replication()
    print(chromosome.to_dict()['domains'])
    assert len(chromosome.domains) == 7

    chromosome.advance_replisomes({0: (7, 5), 1: (3, 4), 2: (8, 9)})
    print('replisomes:')
    print(chromosome.to_dict())

    operon = chromosome.promoters['pA'].choose_operon(chromosome.genes)
    print('operon')
    print(operon.to_dict())

    children = chromosome.terminate_replication()
    print('termination:')
    print([child.to_dict() for child in children])


if __name__ == '__main__':
    test_chromosome()
