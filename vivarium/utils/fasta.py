def read_sequence(path):
    with open(path, 'r') as fasta:
        header = False
        sequence = ''
        for line in fasta:
            if not header:
                header = True
            else:
                sequence += line.strip()
    return sequence
