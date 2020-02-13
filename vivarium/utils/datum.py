def first(l):
    if l:
        return l[0]

def first_value(d):
    if d:
        return d[list(d.keys())[0]]

class Datum(object):
    '''
    The Datum class enables functions to be defined on dicts of a certain schema. 
    It provides two class level variables:
      * `defaults`: a dictionary of keys to default values this Datum will have if 
           none is provided to __init__
      * `schema`: a dictionary of keys to constructors which invoke subdata. 

    Once these are defined, a Datum subclass can be constructed with a dict that provides any
    values beyond the defaults, and then all of the defined methods for that Datum subclass
    are available to operate on its values. Once the modifications are complete, it can be
    rendered back into a dict using the `to_dict()` method.
    '''

    schema = {}
    defaults = {}

    def __init__(self, config, default):
        self.keys = list(set(list(config.keys()) + list(default.keys()))) # a dance

        for key in self.keys:
            value = config.get(key, default[key])
            if value and key in self.schema:
                realize = self.schema[key]
                if isinstance(value, list):
                    value = [realize(item) for item in value]
                elif isinstance(value, dict):
                    value = {inner: realize(item) for inner, item in value.items()}
                else:
                    value = realize(item)
            setattr(self, key, value)

    def fields(self):
        return list(self.defaults.keys())

    def to_dict(self):
        to = {}
        for key in self.keys:
            value = getattr(self, key)
            if isinstance(value, Datum):
                value = value.to_dict()
            elif value and isinstance(value, list) and isinstance(first(value), Datum):
                value = [datum.to_dict() for datum in value]
            elif value and isinstance(value, dict) and isinstance(first_value(value), Datum):
                value = {inner: datum.to_dict() for inner, datum in value.items()}
            to[key] = value
        return to

    def __repr__(self):
        return str(type(self)) + ': ' + str(self.to_dict())

