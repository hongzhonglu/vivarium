from vivarium.utils.fasta import read_sequence
from vivarium.utils.polymerize import generate_template
from vivarium.data.knowledge_base import KnowledgeBase
from vivarium.states.chromosome import Chromosome

ECOLI_GENOME_PATH = 'vivarium/data/flat/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa'


class LacChromosome(object):
    ecoli_sequence = None
    knowledge_base = None

    def __init__(self, parameters={}):
        if self.knowledge_base is None:
            self.knowledge_base = KnowledgeBase()
        if self.ecoli_sequence is None:
            self.ecoli_sequence = read_sequence(ECOLI_GENOME_PATH)

        self.factor_thresholds = {
            ('lacZp1', 'Glucose'): 1e-05,
        }

        self.factor_thresholds.update(parameters.get('thresholds', {}))

        self.chromosome_config = {
            'sequence': self.ecoli_sequence,
            'genes': {
                'lac': ['lacZ', 'lacY', 'lacA'],
            },
            'promoters': {
                'lacZp1': {
                    'id': 'lacZp1',
                    'position': 366343,
                    'direction': -1,
                    'sites': [{'thresholds': {'Glucose': 1e-05}}],
                    'terminators': [
                        {
                            'position': 361209,
                            'strength': 1.0,
                            'products': ['lac']}]},
            },
            'domains': {
                0: {
                    'id': 0,
                    'lead': 0,
                    'lag': 0,
                    'children': []}},
            'rnaps': []}

        # build chromosome and apply thresholds
        self.chromosome = Chromosome(self.chromosome_config)
        self.chromosome.apply_thresholds(self.factor_thresholds)
        self.chromosome_config['promoters'] = {
            key: promoter.to_dict()
            for key, promoter in self.chromosome.promoters.items()}

        self.promoters = self.chromosome_config['promoters'].keys()
        self.promoter_affinities = {
            ('lacZp1', 'Glucose'): 0.01}
        self.promoter_affinities.update(
            parameters.get('promoter_affinities', {}))

        self.transcripts = [
            (operon, product)
            for operon, products in self.chromosome_config['genes'].items()
            for product in products]

        self.protein_sequences = {
            (operon, product): self.knowledge_base.proteins[
                self.knowledge_base.genes[product]['id']]['seq']
            for operon, product in self.transcripts}

        self.transcript_templates = {
            key: generate_template(
                key,
                len(sequence),
                [key[1]])
            for key, sequence in self.protein_sequences.items()}

        self.transcript_affinities = {
            operon: 0.01
            for operon in self.transcripts}

        self.transcript_affinities.update(
            parameters.get('transcript_affinities', {}))

        self.transcription_factors = ['Glucose']
        self.complexation_monomer_ids = []
        self.complexation_complex_ids = []
        self.complexation_stoichiometry = {}
        self.complexation_rates = {}

