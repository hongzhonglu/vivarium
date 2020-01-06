from __future__ import absolute_import, division, print_function


def merge_dicts(dicts):
    merge = {}
    for d in dicts:
        merge.update(d)
    return merge

def flatten_role_dicts(dicts):
    '''
    Input:
        dicts (dict): embedded state dictionaries with the {'role_id': {'state_id': state_value}}
    Return:
        merge (dict): flattened dictionary with {'state_id_role_id': value}
    '''
    merge = {}
    for role, states_dict in dicts.items():
        for state, value in states_dict.items():
            merge.update({state + '_' + role: value})
    return merge
