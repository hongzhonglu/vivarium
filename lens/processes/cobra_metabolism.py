from __future__ import absolute_import, division, print_function

from scipy import constants
from cobra import Model, Reaction, Metabolite

from lens.actor.process import Process, deep_merge
from lens.utils.units import units

def get_reverse(stoichiometry, reversible_reactions, reverse_key):
    '''
    stoichiometry (dict) -- {mol_id (str): stoich (dict)}
    reversible_reactions (list) -- reactions that need a reverse stoichiometry
    '''
    reverse_stoichiometry = {}
    for rxn_id in reversible_reactions:
        forward_stoich = stoichiometry[rxn_id]
        reverse_stoichiometry[rxn_id + reverse_key] = {
            mol_id: -1 * coeff for mol_id, coeff in forward_stoich.iteritems()}
    return reverse_stoichiometry

def build_model(stoichiometry, objective):
    metabolite_keys = {}
    for reaction_key, chemistry in stoichiometry.items():
        metabolite_keys.update(chemistry)

    metabolites = {
        metabolite: Metabolite(metabolite, name=metabolite, compartment='c')
        for metabolite in metabolite_keys.keys()}

    reactions = {}
    for reaction_key, chemistry in stoichiometry.items():
        reaction = Reaction(reaction_key, name=reaction_key)
        chemistry_model = {
            metabolites[metabolite]: value
            for metabolite, value in chemistry.items()}
        reaction.add_metabolites(chemistry_model)
        reactions[reaction_key] = reaction

    model = Model('metabolism')
    model.add_reactions(reactions.values())

    model.objective = {
        reactions[reaction_key]: value
        for reaction_key, value in objective.items()}

    return model

class CobraMetabolism(Process):
    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1/units.mol
        self.density = 1100 * units.g/units.L

        self.initial_state = initial_parameters['initial_state']
        self.stoichiometry = initial_parameters['stoichiometry']

        self.reverse_key = '_reverse'
        self.reversible_reactions = initial_parameters.get('reversible_reactions', [])
        reverse_stoichiometry = get_reverse(self.stoichiometry, self.reversible_reactions, self.reverse_key)
        self.stoichiometry.update(reverse_stoichiometry)

        self.objective = initial_parameters['objective']
        self.stoichiometry.update(self.objective)

        self.objective_declaration = {reaction: 1.0 for reaction in self.objective.keys()}
        self.model = build_model(self.stoichiometry, self.objective_declaration)


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

    objective = {"v_biomass": {"C": 1, "F": 1, "H": 1, "ATP": 10}}

    initial_state = {
        'internal': {
            'mass': 1339,
            'volume': 1E-15},
        'external': {
            "A": 21.0,
            "F": 5.0,
            "D": 12.0,
            "E": 12.0,
            "H": 5.0,
            "O2": 100.0}}

    metabolism = CobraMetabolism({
        'stoichiometry': stoichiometry,
        'objective': objective,
        'initial_state': initial_state})

    return metabolism

if __name__ == '__main__':
    metabolism = test_metabolism()
    print(metabolism.model)
    print(metabolism.model.reactions)
    print(metabolism.model.metabolites)
    print(metabolism.model.genes)
