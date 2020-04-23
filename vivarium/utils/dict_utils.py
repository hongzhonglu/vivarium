from __future__ import absolute_import, division, print_function

import collections
import copy

tuple_separator = '___'

def merge_dicts(dicts):
    merge = {}
    for d in dicts:
        merge.update(d)
    return merge

def deep_merge_check(dct, merge_dct):
    '''
    Recursive dict merge, which throws exceptions for conflicting values
    This mutates dct - the contents of merge_dct are added to dct (which is also returned).
    If you want to keep dct you could call it like deep_merge(dict(dct), merge_dct)'''

    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            try:
                deep_merge_check(dct[k], merge_dct[k])
            except:
                raise Exception('dict merge mismatch: key "{}" has values {} AND {}'.format(k, dct[k], merge_dct[k]))
        elif k in dct and (dct[k] is not merge_dct[k]):
            raise Exception('dict merge mismatch: key "{}" has values {} AND {}'.format(k, dct[k], merge_dct[k]))
        else:
            dct[k] = merge_dct[k]
    return dct

def deep_merge(dct, merge_dct):
    '''
    Recursive dict merge
    This mutates dct - the contents of merge_dct are added to dct (which is also returned).
    If you want to keep dct you could call it like deep_merge(dict(dct), merge_dct)'''
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            deep_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct

def flatten_port_dicts(dicts):
    '''
    Input:
        dicts (dict): embedded state dictionaries with the {'port_id': {'state_id': state_value}}
    Return:
        merge (dict): flattened dictionary with {'state_id_port_id': value}
    '''
    merge = {}
    for port, states_dict in dicts.items():
        for state, value in states_dict.items():
            merge.update({state + '_' + port: value})
    return merge

def tuplify_port_dicts(dicts):
    '''
    Input:
        dicts (dict): embedded state dictionaries with the {'port_id': {'state_id': state_value}}
    Return:
        merge (dict): tuplified dictionary with {(port_id','state_id'): value}
    '''
    merge = {}
    for port, states_dict in dicts.items():
        for state, value in states_dict.items():
            merge.update({(port, state): value})
    return merge

def flatten_timeseries(timeseries):
    '''Flatten a timeseries in the style of flatten_port_dicts'''
    flat = {}
    for port, store_dict in timeseries.items():
        if port == 'time':
            flat[port] = timeseries[port]
            continue
        for variable_name, values in store_dict.items():
            key = "{}_{}".format(port, variable_name)
            flat[key] = values
    return flat

def tuple_to_str_keys(dictionary):
    # take a dict with tuple keys, and convert them to strings with tuple_separator as a delimiter
    new_dict = copy.deepcopy(dictionary)
    make_str_dict(new_dict)
    return new_dict

def make_str_dict(dictionary):
    # get down to the leaves first
    for k, v in dictionary.items():
        if isinstance(v, dict):
            make_str_dict(v)

        # convert tuples in lists
        if isinstance(v, list):
            for idx, var in enumerate(v):
                if isinstance(var, tuple):
                    v[idx] = tuple_separator.join(var)
                if isinstance(var, dict):
                    make_str_dict(var)

    # which keys are tuples?
    tuple_ks = [k for k in dictionary.keys() if isinstance(k, tuple)]
    for tuple_k in tuple_ks:
        str_k = tuple_separator.join(tuple_k)
        dictionary[str_k] = dictionary[tuple_k]
        del dictionary[tuple_k]

    return dictionary

def str_to_tuple_keys(dictionary):
    # take a dict with keys that have tuple_separator, and convert them to tuples

    # get down to the leaves first
    for k, v in dictionary.items():
        if isinstance(v, dict):
            str_to_tuple_keys(v)

        # convert strings in lists
        if isinstance(v, list):
            for idx, var in enumerate(v):
                if isinstance(var, str) and tuple_separator in var:
                    v[idx] = tuple(var.split(tuple_separator))
                if isinstance(var, dict):
                    str_to_tuple_keys(var)

    # which keys are tuples?
    str_ks = [k for k in dictionary.keys() if isinstance(k, str) and tuple_separator in k]
    for str_k in str_ks:
        tuple_k = tuple(str_k.split(tuple_separator))
        dictionary[tuple_k] = dictionary[str_k]
        del dictionary[str_k]

    return dictionary
