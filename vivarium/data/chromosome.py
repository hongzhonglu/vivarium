from vivarium.data.knowledge_base import KnowledgeBase

knowledge_base = KnowledgeBase()

test_chromosome_config = {
    'sequence': 'ATACGGCACGTGACCGTCAACTTA',
    'genes': {
        'oAZ': ['A', 'Z'],
        'oA': ['A'],
        'oB': ['B']},
    'promoter_order': ['pA', 'pB'],
    'promoters': {
        'pA': {
            'id': 'pA',
            'position': 3,
            'direction': 1,
            'sites': [{
                'position': 0,
                'length': 3,
                'thresholds': [
                    ('tfA', 0.3)]}],
            'terminators': [
                {
                    'position': 6,
                    'strength': 0.5,
                    'products': ['oA']},
                {
                    'position': 12,
                    'strength': 1.0,
                    'products': ['oAZ']}]},
        'pB': {
            'id': 'pB',
            'position': -3,
            'direction': -1,
            'sites': [{
                'position': 0,
                'length': 3,
                'thresholds': [
                    ('tfB', 0.5)]}],
            'terminators': [
                {
                    'position': -9,
                    'strength': 0.5,
                    'products': ['oB']},
                {
                    'position': -12,
                    'strength': 1.0,
                    'products': ['oBY']}]}},
    'domains': {
        0: {
            'id': 0,
            'lead': 0,
            'lag': 0,
            'children': []}},
    'rnaps': []}



gfp_plasmid_config = {
    'sequence': 'TAATACGACTCACTATAGGATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG',
    'genes': {
        'GFP_RNA': ['GFP_RNA']},
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
                    'products': ['GFP_RNA']}]}},
    'domains': {
        0: {
            'id': 0,
            'lead': 0,
            'lag': 0,
            'children': []}},
    'rnaps': []}


flagella_genes = knowledge_base.wcEcoli_genes
flagella_sequence = ''

flagella_config = {
    'sequence': flagella_sequence,
    'genes': {},
    'promoters': {},
    'domains': {
        0: {
            'id': 0,
            'lead': 0,
            'lag': 0,
            'children': []}},
    'rnaps': []}
