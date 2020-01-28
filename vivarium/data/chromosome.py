amino_acid_records = [
    {'name': 'Alanine', 'abbreviation': 'Ala', 'symbol': 'A'},
    {'name': 'Arginine', 'abbreviation': 'Arg', 'symbol': 'R'},
    {'name': 'Asparagine', 'abbreviation': 'Asn', 'symbol': 'N'},
    {'name': 'Aspartic acid', 'abbreviation': 'Asp', 'symbol': 'D'},
    {'name': 'Cysteine', 'abbreviation': 'Cys', 'symbol': 'C'},
    {'name': 'Glutamic acid', 'abbreviation': 'Glu', 'symbol': 'E'},
    {'name': 'Glutamine', 'abbreviation': 'Gln', 'symbol': 'Q'},
    {'name': 'Glycine', 'abbreviation': 'Gly', 'symbol': 'G'},
    {'name': 'Histidine', 'abbreviation': 'His', 'symbol': 'H'},
    {'name': 'Isoleucine', 'abbreviation': 'Ile', 'symbol': 'I'},
    {'name': 'Leucine', 'abbreviation': 'Leu', 'symbol': 'L'},
    {'name': 'Lysine', 'abbreviation': 'Lys', 'symbol': 'K'},
    {'name': 'Methionine', 'abbreviation': 'Met', 'symbol': 'M'},
    {'name': 'Phenylalanine', 'abbreviation': 'Phe', 'symbol': 'F'},
    {'name': 'Proline', 'abbreviation': 'Pro', 'symbol': 'P'},
    {'name': 'Serine', 'abbreviation': 'Ser', 'symbol': 'S'},
    {'name': 'Threonine', 'abbreviation': 'Thr', 'symbol': 'T'},
    {'name': 'Tryptophan', 'abbreviation': 'Trp', 'symbol': 'W'},
    {'name': 'Tyrosine', 'abbreviation': 'Tyr', 'symbol': 'Y'},
    {'name': 'Valine', 'abbreviation': 'Val', 'symbol': 'V'}]

amino_acids = {
    record['symbol']: record
    for record in amino_acid_records}

test_chromosome_config = {
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
                    'position': -8,
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
            'state': 'transcribing',
            'position': 3},
        {
            'promoter': 'pA',
            'domain': 0,
            'state': 'transcribing',
            'position': 6},
        {
            'promoter': 'pA',
            'domain': 0,
            'state': 'transcribing',
            'position': 0}]}


gfp_plasmid_config = {
    'sequence': 'TAATACGACTCACTATAGGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG',
    'genes': {
        'GFP': ['GFP']},
    'promoters': {
        'T7': {
            'id': 'T7',
            'position': 19,
            'direction': 1,
            'sites': [],
            'terminators': [
                {
                    'position': 736,
                    'strength': 1.0,
                    'operon': 'GFP'}]},
        },
    'domains': {
        0: {
            'id': 0,
            'lead': 0,
            'lag': 0,
            'children': []}},
    'rnaps': []}
