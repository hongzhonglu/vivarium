import copy

from vivarium.actor.process import Process
from vivarium.data.chromosome import amino_acid_records

class Translation(Process):
    def __init__(self, initial_parameters={}):
        monomer_ids = [record.symbol for record in amino_acid_records]
        self.default_parameters = {
            'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
            'elongation_rate': 1.0,
            'advancement_rate': 1.0,
            'monomer_ids': monomer_ids}
        self.default_parameters['transcript_order'] = list(
            self.default_parameters['transcript_affinities'].keys())

        parameters = copy.deepcopy(self.default_parameters)
        parameters.update(initial_parameters)

        self.transcript_affinities = parameters['transcript_affinities']
        self.transcript_order = parameters['transcript_order']
        self.transcript_count = len(self.transcript_order)

        self.molecule_ids = parameters['molecule_ids']
        self.monomer_ids = parameters['monomer_ids']
        self.transcript_ids = parameters['transcript_ids']
        self.elongation = 0
        self.elongation_rate = parameters['elongation_rate']
        self.advancement_rate = parameters['advancement_rate']
