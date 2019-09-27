"""
Subclasses of DictWriter and DictReader that parse plaintext as JSON strings,
allowing for basic type parsing and fields that are dictionaries or lists.
"""

import os
import csv
import json
import re
import numpy as np
from itertools import ifilter

from lens.utils.units import units

TSV_DIALECT = csv.excel_tab

def array_to_list(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()

    return value


class JsonWriter(csv.DictWriter):
    def __init__(self, *args, **kwargs):
        csv.DictWriter.__init__(
            self, quotechar = "'", quoting = csv.QUOTE_MINIMAL, lineterminator="\n", *args, **kwargs
            )

    def _dict_to_list(self, rowdict):
        return csv.DictWriter._dict_to_list(self, {
            key:json.dumps(array_to_list(value))
            for key, value in rowdict.viewitems()
            })


class JsonReader(csv.DictReader):
    def __init__(self, *args, **kwargs):
        csv.DictReader.__init__(
            self, quotechar = "'", quoting = csv.QUOTE_MINIMAL, *args, **kwargs
            )

        # This is a hack to strip extra quotes from the field names
        # Not proud of it, but it works.
        self.fieldnames # called for side effect

        self._fieldnames = [
            fieldname.strip('"') for fieldname in self._fieldnames
            ]

    def next(self):
        attributeDict = {}
        for key, raw_value in csv.DictReader.next(self).viewitems():
            try:
                value = json.loads(raw_value) if raw_value else ""

            except (ValueError, TypeError) as e:
                repr(e)
                raise Exception("failed to parse json string:{}".format(raw_value))

            try:
                attribute = re.search('(.*?) \(', key).group(1)
                value_units =  eval(re.search('\((.*?)\)',key).group(1))
                attributeDict[attribute] = value * value_units
            except AttributeError:
                attributeDict[key] = value
        return attributeDict

        # return {
        # 	key:json.loads(value) if value else "" # catch for empty field
        # 	for key, value in csv.DictReader.next(self).viewitems()
        # 	}

def load_tsv(dir_name, file_name):
    file_path = os.path.join(dir_name, file_name)
    with open(file_path, 'rU') as tsvfile:
        reader = JsonReader(
            ifilter(lambda x: x.lstrip()[0] != "#", tsvfile),  # Strip comments
            dialect=TSV_DIALECT)
        attr_list = []
        fieldnames = reader.fieldnames
        for row in reader:
            attr_list.append({field: row[field] for field in fieldnames})
    return attr_list