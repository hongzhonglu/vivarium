
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