from vivarium.utils.fasta import read_sequence
from vivarium.utils.polymerize import generate_template
from vivarium.data.knowledge_base import KnowledgeBase

knowledge_base = KnowledgeBase()

ECOLI_GENOME_PATH = 'vivarium/data/flat/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa'

ecoli_sequence = read_sequence(ECOLI_GENOME_PATH)

class FlagellaChromosome(object):
    def __init__(self):
        self.config = {
            'sequence': ecoli_sequence,
            'genes': {
                'flhDC': ['flhD', 'flhC'],
                'fliL': ['fliL', 'fliM', 'fliN', 'fliO', 'fliP', 'fliQ', 'fliR'],
                'fliE': ['fliE'],
                'fliF': ['fliF', 'fliG', 'fliH', 'fliI', 'fliJ', 'fliK'],
                'flgA': ['flgA', 'flgM', 'flgN'],
                'flgB': ['flgB', 'flgC', 'flgD', 'flgE', 'flgF', 'flgG', 'flgH', 'flgI', 'flgJ'],
                'flhB': ['flhB', 'flhA', 'flhE'],
        'fliA': ['fliA', 'fliZ'], # ignore 'tcyJ' for now
                'fliD': ['fliD', 'fliS', 'fliT'],
                'flgK': ['flgK', 'flgL'],
                'fliC': ['fliC'],
                # 'meche': [],
                # 'mocha': [],
                'flgM': ['flgM', 'flgN']},
            'promoters': {
                'flhDp': {
                    'id': 'flhDp',
                    'position': 1978197,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('CRP', 0.1)]}],
                    'terminators': [
                        {
                            'position': 1977266,
                            'strength': 1.0,
                            'products': ['flhDC']}]},
                'fliLp1': {
                    'id': 'fliLp1',
                    'position': 2019618,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.1)]}],
                    'terminators': [
                        {
                            'position': 2023678,
                            'strength': 1.0,
                            'products': ['fliL']}]},
                'fliEp1': {
                    'id': 'fliEp1',
                    'position': 2013014,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.15)]}],
                    'terminators': [
                        {
                            'position': 2012700,
                            'strength': 1.0,
                            'products': ['fliE']}]},
                'fliFp1': {
                    'id': 'fliFp1',
                    'position': 2013229,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.2)]}],
                    'terminators': [
                        {
                            'position': 2019513,
                            'strength': 1.0,
                            'products': ['fliF']}]},
                'flgAp': {
                    'id': 'flgAp',
                    'position': 1130863,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.25)]}],
                    'terminators': [
                        {
                            'position': 1129414,
                            'strength': 1.0,
                            'products': ['flgA']}]},
                'flgBp': {
                    'id': 'flgBp',
                    'position': 1131018,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.3)]}],
                    'terminators': [
                        {
                            'position': 1138312,
                            'strength': 1.0,
                            'products': ['flgB']}]},
                'flhBp': {
                    'id': 'flhBp',
                    'position': 1966191,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.35)]}],
                    'terminators': [
                        {
                            'position': 1962580,
                            'strength': 1.0,
                            'products': ['flhB']}]},
                'fliAp1': {
                    'id': 'fliAp1',
                    'position': 2001789,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.4)]}],
                    'terminators': [
                        {
                            'position': 1999585,
                            'strength': 1.0,
                            'products': ['fliA']}]},
                'fliDp': {
                    'id': 'fliDp',
                    'position': 2003872,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.45)]}],
                    'terminators': [
                        {
                            'position': 2006078,
                            'strength': 1.0,
                            'products': ['fliD']}]},
                'flgKp': {
                    'id': 'flgKp',
                    'position': 1138378,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('flhD', 0.5)]}],
                    'terminators': [
                        {
                            'position': 1140986,
                            'strength': 1.0,
                            'products': ['flgK']}]},
                'fliCp': {
                    'id': 'fliCp',
                    'position': 2002110,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('GadE', 0.55),
                                ('H-NS', 0.6)]}],
                    'terminators': [
                        {
                            'position': 2003606,
                            'strength': 1.0,
                            'products': ['fliC']}]},
                # 'meche': {
                #     'id': 'meche',
                #     'position': 0,
                #     'direction': 1,
                #     'sites': [],
                #     'terminators': []},
                # 'mocha': {
                #     'id': 'mocha',
                #     'position': 0,
                #     'direction': 1,
                #     'sites': [],
                #     'terminators': []},
                'flgMp': {
                    'id': 'flgMp',
                    'position': 1130128,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('CsgD', 0.1)]}],
                    'terminators': [
                        {
                            'position': 1129414,
                            'strength': 1.0,
                            'products': ['flgM']}]}},
            'domains': {
                0: {
                    'id': 0,
                    'lead': 0,
                    'lag': 0,
                    'children': []}},
            'rnaps': []}

        self.promoters = [
            'flhDp',
            'fliLp1',
            'fliEp1',
            'fliFp1',
            'flgAp',
            'flgBp',
            'flhBp',
            'fliAp1',
            'fliDp',
            'flgKp',
            'fliCp',
            # 'meche',
            # 'mocha',
            'flgMp']

        self.flhD_activated = [
            'fliLp1',
            'fliEp1',
            'fliFp1',
            'flgAp',
            'flgBp',
            'flhBp',
            'fliAp1',
            'fliDp',
            'flgKp']

        self.promoter_affinities = {
            ('flhDp', None): 0.01,
            ('flhDp', 'CRP'): 0.5}

        self.promoter_affinities[('fliCp', None)] = 0.0
        self.promoter_affinities[('flgMp', None)] = 0.0

        for promoter in self.flhD_activated:
            self.promoter_affinities[(promoter, None)] = 0.0
            self.promoter_affinities[(promoter, 'flhD')] = 0.5

        self.transcripts = [
            product
            for products in self.config['genes'].values()
            for product in products]

        self.protein_sequences = {
            symbol: knowledge_base.proteins[knowledge_base.genes[symbol]['id']]
            for symbol in self.transcripts}

        self.operon_sequences = {
            operon: knowledge_base.concatenate_sequences(units)
            for operon, units in self.config['genes'].items()}

        self.transcript_templates = {
            operon: generate_template(
                operon,
                len(sequence),
                self.config['genes'][operon])
            for operon, sequence in self.operon_sequences.items()}

        self.transcript_affinities = {
            operon: 1.0
            for operon in self.config['genes'].keys()}

        self.transcription_factors = ['flhD', 'CsgD', 'CRP', 'GadE', 'H-NS']
