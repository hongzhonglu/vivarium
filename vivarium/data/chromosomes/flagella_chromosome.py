
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
                'tarp': ['tar', 'tap', 'cheR', 'cheB', 'cheY', 'cheZ'],
                'motAp': ['motA', 'motB', 'cheA', 'cheW'],
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
                                ('CRP', 1e-05)]}],
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
                                ('flhDC', 1e-05)]}],
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
                                ('flhDC', 1.5e-05)]}],
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
                                ('flhDC', 2e-05)]}],
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
                                ('flhDC', 2.5e-05)]}],
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
                                ('flhDC', 3e-05)]}],
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
                                ('flhDC', 3.5e-05)]}],
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
                                ('flhDC', 4e-05)]}],
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
                                ('flhDC', 4.5e-05)]}],
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
                                ('flhDC', 5e-05)]}],
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
                'tarp': {
                    'id': 'tarp',
                    'position': 1972691,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('Fnr', 1e-5)]}],
                    'terminators': [
                        {
                            'position': 1971030,
                            'strength': 1.0,
                            'products': ['tarp']}]},
                'motAp': {
                    'id': 'motAp',
                    'position': 1977139,
                    'direction': -1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('CpxR', 1e-5)]}],
                    'terminators': [
                        {
                            'position': 1976252,
                            'strength': 1.0,
                            'products': ['motAp']}]},
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
            'tarp',
            'motAp'
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
            ('flhDp', None): 0.0,
            ('flhDp', 'CRP'): 0.01}

        self.promoter_affinities[('fliCp', None)] = 0.0
        self.promoter_affinities[('flgMp', None)] = 0.0
        self.promoter_affinities[('tarp', None)] = 0.0
        self.promoter_affinities[('tarp', 'Fnr')] = 0.0
        self.promoter_affinities[('motAp', None)] = 0.0
        self.promoter_affinities[('motAp', 'CpxR')] = 1.0

        for promoter in self.flhD_activated:
            self.promoter_affinities[(promoter, None)] = 0.0
            self.promoter_affinities[(promoter, 'flhDC')] = 1.0

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
            operon: 0.01
            for operon in self.config['genes'].keys()}

        self.transcription_factors = [
            'flhDC', 'CsgD', 'CRP', 'GadE', 'H-NS', 'CpxR', 'Fnr']

        self.complexation_monomer_ids = [
            'fliG', 'fliM', 'fliN', 'flhA', 'flhB', 'flhD', 'flhC', 'fliO', 'fliP', 'fliQ', 'fliR', 'fliJ', 'fliI', 'fliH', 'fliL', 'flgH', 'motA', 'motB', 'flgB', 'flgC', 'flgF', 'flgG', 'flgI', 'fliF', 'fliE','fliC','flgL','flgK','fliD','flgE']

        self.complexation_complex_ids = [
            'flhDC',
            'flagellar motor switch',
            'flagellum',
            'flagellar export apparatus',
            'flagellar motor']

        self.complexation_stoichiometry = {
            'flhDC': {
                'flhD': -4.0,
                'flhC': -2.0,
                'flhDC': 1.0},
            'flagellar motor switch reaction': {
                'flagellar motor switch': 1.0,
                'fliG': -26.0,
                'fliM': -34.0,
                'fliN': -1.0},
            'flagellar export apparatus reaction': {
                'flagellar export apparatus': 1.0,
                'flhA': -1.0,
                'flhB': -1.0,
                'fliO': -1.0,
                'fliP': -1.0,
                'fliQ': -1.0,
                'fliR': -1.0,
                'fliJ': -1.0,
                'fliI': -6.0,
                'fliH': -12.0},
            'flagellar motor reaction': {
                'flagellar motor': 1.0,
                'flagellar motor switch': -1.0,
                'fliL': -2.0,
                'flgH': -1.0,
                'motA': -1.0,
                'motB': -1.0,
                'flgB': -1.0,
                'flgC': -1.0,
                'flgF': -1.0,
                'flgG': -1.0,
                'flgI': -1.0,
                'fliF': -1.0,
                'fliE': -1.0},
            'flagellum reaction': {
                'flagellum': 1.0,
                'flagellar export apparatus': -1.0,
                'flagellar motor': -1.0,
                'fliC': -1.0,
                'flgL': -1.0,
                'flgK': -1.0,
                'fliD': -5.0,
                'flgE': -120.0}}

        reaction_default = 1e-5
        self.complexation_rates = {
            'flhDC': reaction_default,
            'flagellar motor switch reaction': reaction_default,
            'flagellar export apparatus reaction': reaction_default,
            'flagellar motor reaction': reaction_default,
            'flagellum reaction': reaction_default}
