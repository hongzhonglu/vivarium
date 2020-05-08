import json
import numpy as np

def _json_serialize(elem):
    if isinstance(elem, np.int64):
        return int(elem)
    raise TypeError(
        "Objectof type {} is not JSON serializable".format(type(elem))
    )

def format_dict(d, sort_keys=True):
    '''Format a dict as a pretty string

    Aside from the normal JSON-serializable data types, data of type
    ``numpy.int64`` are supported.

    For example:

    >>> import numpy as np
    >>> d = {
    ...     'foo': {
    ...         'bar': 1,
    ...         '3.0': np.int64(5),
    ...     },
    ...     'a': 'hi!',
    ... }
    >>> print(format_dict(d))
    {
        "a": "hi!",
        "foo": {
            "3.0": 5,
            "bar": 1
        }
    }

    Arguments:
        d: The dictionary to format
        sort_keys: Whether to sort the dictionary keys. This is useful
            for reproducible output.

    Returns:
        A string of the prettily-formatted dictionary

    '''
    return json.dumps(
        d, indent=4, default=_json_serialize, sort_keys=sort_keys)
