from __future__ import absolute_import, division, print_function

from cobra import Model, Reaction, Metabolite

from lens.actor.process import Process, deep_merge

def build_model(stoichiometry):
    metabolite_keys = {}
    for reaction_key, chemistry in stoichiometry.items():
        metabolite_keys.update(chemistry)

    metabolites = {
        metabolite: Metabolite(metabolite, name=metabolite, compartment='c')
        for metabolite in metabolite_keys.keys()}

    reactions = []
    for reaction_key, chemistry in stoichiometry.items():
        reaction = Reaction(reaction_key, name=reaction_key)
        chemistry_model = {
            metabolites[metabolite]: value
            for metabolite, value in chemistry.items()}
        reaction.add_metabolites(chemistry_model)
        reactions.append(reaction)

    model = Model('metabolism')
    model.add_reactions(reactions)

    return model

class CobraMetabolism(Process):
    def __init__(self, initial_parameters={}):
        pass


def test_metabolism():
    stoichiometry = {
        "R1": {"A": -1, "ATP": -1, "B": 1},
        "R2a": {"B": -1, "ATP": 2, "NADH": 2, "C": 1},
        "R2b": {"C": -1, "ATP": -2, "NADH": -2, "B": 1},
        "R3": {"B": -1, "F": 1},
        "R4": {"C": -1, "G": 1},
        "R5": {"G": -1, "C": 0.8, "NADH": 2},
        "R6": {"C": -1, "ATP": 2, "D": 3},
        "R7": {"C": -1, "NADH": -4, "E": 3},
        "R8a": {"G": -1, "ATP": -1, "NADH": -2, "H": 1},
        "R8b": {"G": 1, "ATP": 1, "NADH": 2, "H": -1},
        "Rres": {"NADH": -1, "O2": -1, "ATP": 1}}

    model = build_model(stoichiometry)
    return model

if __name__ == '__main__':
    model = test_metabolism()
    print(model)
    print(model.reactions)
    print(model.metabolites)
    print(model.genes)
