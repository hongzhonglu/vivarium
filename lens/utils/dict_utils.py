from __future__ import absolute_import, division, print_function


def merge_dicts(dicts):
    merge = {}
    for d in dicts:
        merge.update(d)
    return merge
