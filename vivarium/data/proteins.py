from vivarium.utils.datum import Datum

class Protein(Datum):
    defaults = {
        'id': '',
        'sequence': ''}

    def __init__(self, config, defaults=defaults):
        super(Protein, self).__init__(config, self.defaults)

GFP = Protein({
    'id': 'GFP',
    'sequence': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'})

