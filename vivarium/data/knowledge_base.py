import os
import csv

from vivarium.data.spreadsheets import load_tsv

FLAT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "flat")

LIST_OF_FLAT_FILENAMES = (
    os.path.join("wcEcoli_environment_molecules.tsv"),
    os.path.join("timelines_def.tsv"),
    os.path.join("media_recipes.tsv"),
    os.path.join("media", "wcEcoli_base.tsv"),
    os.path.join("media", "M9.tsv"),
    os.path.join("media", "M9_GLC.tsv"),
    os.path.join("media", "5X_supplement_EZ.tsv"),
    os.path.join("media", "GLC_G6P.tsv"),
    os.path.join("media", "GLC_LCT.tsv"),
    os.path.join("media", "ecoli_core_GLC.tsv"),
)

class DataStore(object):
    def __init__(self):
        pass

class KnowledgeBase(object):
    """ KnowledgeBase """

    def __init__(self):
        # Load raw data from TSV files

        for filename in LIST_OF_FLAT_FILENAMES:
            self._load_tsv(FLAT_DIR, os.path.join(FLAT_DIR, filename))

    def _load_tsv(self, dir_name, file_name):
        path = self
        for subPath in file_name[len(dir_name) + 1 : ].split(os.path.sep)[:-1]:
            if not hasattr(path, subPath):
                setattr(path, subPath, DataStore())
            path = getattr(path, subPath)
        attrName = file_name.split(os.path.sep)[-1].split(".")[0]
        setattr(path, attrName, [])

        file_path = os.path.join(dir_name, file_name)
        rows = load_tsv(file_path)
        setattr(path, attrName, [row for row in rows])
