def invert_index(mapping):
    inverse = {}
    for canon, synonyms in mapping.items():
        for synonym in synonyms:
            inverse[synonym] = canon
    return inverse

synonyms = {
    'rATP': ['A', 'ATP'],
    'rGTP': ['G', 'GTP']}

index = invert_index(synonyms)

def convert_ids(ids):
    return [
        index[id]
        for id in ids]

def test_synonyms():
    pass
