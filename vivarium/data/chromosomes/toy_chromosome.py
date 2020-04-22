toy_chromosome_config = {
    'sequence': 'ATACGGCACGTGACCGTCAACTTA',
    'genes': {
        'oA': ['eA'],
        'oAZ': ['eA', 'eZ'],
        'oB': ['eB'],
        'oBY': ['eB', 'eY']},
    'promoter_order': ['pA', 'pB'],
    'promoters': {
        'pA': {
            'id': 'pA',
            'position': 3,
            'direction': 1,
            'sites': [{
                'position': 0,
                'length': 3,
                'thresholds': {'tfA': 0.3}}],
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
                'thresholds': {'tfB': 0.5}}],
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
