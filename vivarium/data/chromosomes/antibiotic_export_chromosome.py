from vivarium.data.knowledge_base import KnowledgeBase
from vivarium.utils.fasta import read_sequence
from vivarium.utils.constants import ECOLI_GENOME_PATH
from vivarium.utils.polymerize import generate_template


class AntibioticExportChromosome(object):
    def __init__(self):
        knowledge_base = KnowledgeBase()
        ecoli_sequence = read_sequence(ECOLI_GENOME_PATH)
        self.config = {
            'sequence': ecoli_sequence,
            'genes': {
                'acrR': ['acrR'],
                'acrAB': ['acrA', 'acrB'],
                'marRAB': ['marR', 'marA', 'marB'],
                'tolC-ygiBC': ['tolC', 'ygiB', 'ygiC'],
            },
            'promoters': {
                # EcoCyc warns no high-quality experimental evidence
                # for this promoter
                'acrR': {
                    'id': 'acrR',
                    # Promoter positions are set at TSS
                    'position': 485699,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('acrR', 1.0e-5),
                            ],
                        },
                    ],
                    'terminators': [
                        {
                            'position': 486408,
                            'strength': 1.0,
                            'products': ['acrR'],
                        },
                    ],
                },
                'acrAB': {
                    'id': 'acrAB',
                    'position': 485698,
                    'direction': -1,
                    'sites': [
                        {
                            # Promoter binding site location information
                            # is present in EcoCyc but excluded here
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                # Tuples are (TF, threshold). When the
                                # concentration of the transcription
                                # factor (TF) is above threshold, the TF
                                # binds this site.
                                ('marA', 1.0e-5),
                                ('acrR', 1.0e-5),
                            ],
                        },
                    ],
                    'terminators': [
                        {
                            # Terminator set inaccurately at end of
                            # operon because downstream genes not
                            # modeled
                            'position': 481254,
                            'strength': 1.0,
                            # The transcripts produced by ending at this
                            # terminator
                            'products': ['acrAB'],
                        },
                    ],
                },
                'marRAB': {
                    'id': 'marRAB',
                    'position': 1619093,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('marA', 1.0e-5),
                                ('marR', 1.0e-5),
                                ('acrR', 1.0e-5),
                            ],
                        }
                    ],
                    'terminators': [
                        {
                            # Inaccurately set at end of operon
                            'position': 1620207,
                            'strength': 1.0,
                            'products': ['marRAB'],
                        },
                    ],
                },
                'tolC-ygiBC': {
                    # There are actually 4 TSSs for this operon. Only
                    # one is modeled
                    'id': 'tolC-ygiBC',
                    'position': 3178074,
                    'direction': 1,
                    'sites': [
                        {
                            'position': 0,
                            'length': 0,
                            'thresholds': [
                                ('marA', 1.0e-5),
                            ],
                        },
                    ],
                    'terminators': [
                        {
                            # Set innacurately at end of operon
                            'position': 3183581,
                            'strength': 1.0,
                            'products': ['tolC-ygiBC'],
                        },
                    ],
                },
            },
            'domains': {
                0: {
                    'id': 0,
                    'lead': 0,
                    'lag': 0,
                    'children': [],
                },
            },
            'rnaps': [],
        }

        self.promoters = self.config['promoters'].keys()
        self.promoter_affinities = {
            ('acrR', None): 1.0,
            ('acrR', 'acrR'): 0.1,
            ('acrAB', None): 1.0,
            ('acrAB', 'marA'): 2.0,
            ('acrAB', 'acrR'): 0.1,
            ('marRAB', None): 1.0,
            ('marRAB', 'marA'): 2.0,
            ('marRAB', 'marR'): 0.1,
            ('marRAB', 'acrR'): 0.1,
            ('tolC-ygiBC', None): 1.0,
            ('tolC-ygiBC', 'marA'): 2.0,
        }
        self.transcripts = [
            product
            for products in self.config['genes'].values()
            for product in products
        ]
        self.protein_sequences = {
            symbol: knowledge_base.proteins[
                knowledge_base.genes[symbol]['id']]
            for symbol in self.transcripts
        }
        self.operon_sequences = {
            operon: knowledge_base.concatenate_sequences(units)
            for operon, units in self.config['genes'].items()
        }
        self.transcript_templates = {
            operon: generate_template(
                operon,
                len(sequence),
                self.config['genes'][operon]
            )
            for operon, sequence in self.operon_sequences.items()
        }
        self.transcript_affinities = {
            operon: 0.01
            for operon in self.config['genes']
        }
        self.transcription_factors = ['marA', 'acrR', 'marR']
        self.complexation_monomer_ids = ['acrA', 'acrB', 'tolC']
        self.complexation_complex_ids = ['acrAB-tolC']
        self.complexation_stoichiometry = {
            'acrAB-tolC Assembly': {
                # From EcoCyc, ID TRANS-CPLX-201
                'acrAB-tolC': 1.0,
                'acrA': -6.0,
                'acrB': -3.0,
                'tolC': -3.0
            },
        }
        self.complexation_rates = {
            'acrAB-tolC Assembly': 1.0
        }
