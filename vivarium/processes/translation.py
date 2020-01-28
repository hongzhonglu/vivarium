from vivarium.actor.process import Process
from vivarium.data.chromosome import amino_acids

class Translation(Process):
    def __init__(self, initial_parameters={}):
        monomer_ids = list(amino_acids.keys())
        self.default_parameters = {
            'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
            'elongation_rate': 1.0,
            'advancement_rate': 1.0,
            'monomer_ids': monomer_ids}
