from vivarium.data.knowledge_base import KnowledgeBase

knowledge_base = KnowledgeBase()

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
                    'product': ['oA']},
                {
                    'position': 11,
                    'strength': 1.0,
                    'product': ['oAZ']}]},
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
                    'product': ['oB']},
                {
                    'position': -11,
                    'strength': 1.0,
                    'product': ['oBY']}]}},
    'domains': {
        0: {
            'id': 0,
            'lead': 0,
            'lag': 0,
            'children': []}},
    'rnaps': [
        {
            'template': 'pA',
            'domain': 0,
            'state': 'transcribing',
            'position': 3},
        {
            'template': 'pA',
            'domain': 0,
            'state': 'transcribing',
            'position': 6},
        {
            'template': 'pA',
            'domain': 0,
            'state': 'transcribing',
            'position': 0}]}


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
                    'product': ['GFP_RNA']}]}},
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
