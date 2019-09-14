
import os
import csv

from lens.reconstruction.spreadsheets import JsonReader
from itertools import ifilter

CSV_DIALECT = csv.excel_tab

ENV_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "environment")


LIST_OF_ENV_FILENAMES = (
    os.path.join("condition", "environment_molecules.tsv"),
    os.path.join("condition", "timelines_def.tsv"),
    os.path.join("condition", "media_recipes.tsv"),
    os.path.join("condition", "media", "wcEcoli_base.tsv"),
    os.path.join("condition", "media", "M9.tsv"),
    os.path.join("condition", "media", "M9_GLC.tsv"),
    os.path.join("condition", "media", "5X_supplement_EZ.tsv"),
    os.path.join("condition", "media", "GLC_G6P.tsv"),
    os.path.join("condition", "media", "GLC_LCT.tsv"),
    )

class DataStore(object):
    def __init__(self):
        pass

class KnowledgeBase(object):
    """ KnowledgeBase """

    def __init__(self):
        # Load raw data from TSV files

        for filename in LIST_OF_ENV_FILENAMES:
            self._load_tsv(ENV_DIR, os.path.join(ENV_DIR, filename))

    def _load_tsv(self, dir_name, file_name):
        path = self
        for subPath in file_name[len(dir_name) + 1 : ].split(os.path.sep)[:-1]:
            if not hasattr(path, subPath):
                setattr(path, subPath, DataStore())
            path = getattr(path, subPath)
        attrName = file_name.split(os.path.sep)[-1].split(".")[0]
        setattr(path, attrName, [])

        with open(file_name, 'rU') as csvfile:
            reader = JsonReader(
                ifilter(lambda x: x.lstrip()[0] != "#", csvfile), # Strip comments
                dialect = CSV_DIALECT)
            setattr(path, attrName, [row for row in reader])