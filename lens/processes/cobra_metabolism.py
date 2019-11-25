from __future__ import absolute_import, division, print_function

import cobra
import cobra.test

from scipy import constants
from cobra import Model, Reaction, Metabolite, Configuration

from lens.actor.process import Process, deep_merge
from lens.utils.units import units

EXTERNAL_SUFFIX = '_external'

def get_reverse(stoichiometry, reversible_reactions, reverse_key):
    '''
    stoichiometry (dict) -- {mol_id (str): stoich (dict)}
    reversible_reactions (list) -- reactions that need a reverse stoichiometry
    '''
    reverse_stoichiometry = {}
    for rxn_id in reversible_reactions:
        forward_stoich = stoichiometry[rxn_id]
        reverse_stoichiometry[rxn_id + reverse_key] = {
            mol_id: -1 * coeff for mol_id, coeff in forward_stoich.items()}
    return reverse_stoichiometry

def build_model(stoichiometry, objective, external_molecules):
    model = Model('metabolism')
    model.compartments = {'c': 'cytoplasm'}

    metabolite_keys = {}
    for reaction_key, chemistry in stoichiometry.items():
        metabolite_keys.update(chemistry)

    metabolites = {
        metabolite: Metabolite(metabolite, name=metabolite, compartment='c')
        for metabolite in metabolite_keys.keys()}

    model.add_metabolites(metabolites.values())

    reactions = {}
    for reaction_key, chemistry in stoichiometry.items():
        reaction = Reaction(reaction_key, name=reaction_key)
        reaction_model = {
            metabolites[metabolite]: value
            for metabolite, value in chemistry.items()}
        reaction.add_metabolites(reaction_model)
        # reaction.gene_reaction_rule = 'X'
        reactions[reaction_key] = reaction

    for external in external_molecules:
        external_key = external + EXTERNAL_SUFFIX
        reaction = Reaction(external_key, name=external_key)
        reaction_model = {metabolites[external]: -1}
        reaction.add_metabolites(reaction_model)
        reactions[external_key] = reaction

    model.add_reactions(reactions.values())

    model.objective = {
        reactions[reaction_key]: value
        for reaction_key, value in objective.items()}

    return model

class CobraMetabolism(Process):
    cobra_configuration = Configuration()

    def __init__(self, initial_parameters={}):
        self.nAvogadro = constants.N_A * 1/units.mol
        self.density = 1100 * units.g/units.L

        self.initial_state = initial_parameters['initial_state']
        self.stoichiometry = initial_parameters['stoichiometry']
        self.external_molecules = initial_parameters['external_molecules']

        self.reverse_key = '_reverse'
        self.reversible_reactions = initial_parameters.get('reversible_reactions', [])
        reverse_stoichiometry = get_reverse(self.stoichiometry, self.reversible_reactions, self.reverse_key)
        self.stoichiometry.update(reverse_stoichiometry)

        self.objective = initial_parameters['objective']
        self.stoichiometry.update(self.objective)

        self.objective_declaration = {reaction: 1.0 for reaction in self.objective.keys()}
        self.model = build_model(
            self.stoichiometry,
            self.objective_declaration,
            self.external_molecules)

    def set_external_levels(self, levels):
        for external, level in levels.items():
            reaction = self.model.reactions.get_by_id(external + EXTERNAL_SUFFIX)
            reaction.lower_bound = -level


def test_minimal():
    stoichiometry = {
        'R1': {'A': -1, 'B': 1}}

    objective = {
        'EB': {'B': -1}}

    external_molecules = ['A']

    initial_state = {
        'A': 5}

    metabolism = CobraMetabolism({
        'stoichiometry': stoichiometry,
        'objective': objective,
        'external_molecules': external_molecules,
        'initial_state': initial_state})

    metabolism.set_external_levels(initial_state)

    return metabolism

def test_metabolism():
    stoichiometry = {
        'R1': {'A': -1, 'ATP': -1, 'B': 1},
        'R2a': {'B': -1, 'ATP': 2, 'NADH': 2, 'C': 1},
        'R2b': {'C': -1, 'ATP': -2, 'NADH': -2, 'B': 1},
        'R3': {'B': -1, 'F': 1},
        'R4': {'C': -1, 'G': 1},
        'R5': {'G': -1, 'C': 0.8, 'NADH': 2},
        'R6': {'C': -1, 'ATP': 2, 'D': 3},
        'R7': {'C': -1, 'NADH': -4, 'E': 3},
        'R8a': {'G': -1, 'ATP': -1, 'NADH': -2, 'H': 1},
        'R8b': {'G': 1, 'ATP': 1, 'NADH': 2, 'H': -1},
        'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2', 'BIOMASS']

    objective = {'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'BIOMASS': 1}}

    initial_state = {
        'internal': {
            'mass': 1339,
            'volume': 1E-15},
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0}}

    metabolism = CobraMetabolism({
        'stoichiometry': stoichiometry,
        'objective': objective,
        'external_molecules': external_molecules,
        'initial_state': initial_state})

    metabolism.set_external_levels(initial_state['external'])

    return metabolism

class JsonFBA(object):
    def __init__(self, path):
        self.model = cobra.io.load_json_model(path)

def test_canonical():
    metabolism = JsonFBA('models/e_coli_core.json')
    return metabolism

def test_test():
    metabolism = JsonFBA('models/minimal_model.json')
    return metabolism

def test_demo():
    model = Model('example_model')

    reaction = Reaction('3OAS140')
    reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default

    ACP_c = Metabolite(
        'ACP_c',
        formula='C11H21N2O7PRS',
        name='acyl-carrier-protein',
        compartment='c')
    omrsACP_c = Metabolite(
        '3omrsACP_c',
        formula='C25H45N2O9PRS',
        name='3-Oxotetradecanoyl-acyl-carrier-protein',
        compartment='c')
    co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
    malACP_c = Metabolite(
        'malACP_c',
        formula='C14H22N2O10PRS',
        name='Malonyl-acyl-carrier-protein',
        compartment='c')
    h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
    ddcaACP_c = Metabolite(
        'ddcaACP_c',
        formula='C23H43N2O8PRS',
        name='Dodecanoyl-ACP-n-C120ACP',
        compartment='c')    

    reaction.add_metabolites({
        malACP_c: -1.0,
        h_c: -1.0,
        ddcaACP_c: -1.0,
        co2_c: 1.0,
        ACP_c: 1.0,
        omrsACP_c: 1.0
    })

    print(reaction.reaction)  # This gives a string representation of the reaction

    reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
    print(reaction.genes)

    model.add_reactions([reaction])

    # Now there are things in the model
    print('%i reaction' % len(model.reactions))
    print('%i metabolites' % len(model.metabolites))
    print('%i genes' % len(model.genes))

    # Iterate through the the objects in the model
    print('Reactions')
    print('---------')
    for x in model.reactions:
        print('%s : %s' % (x.id, x.reaction))

    print('')
    print('Metabolites')
    print('-----------')
    for x in model.metabolites:
        print('%9s : %s' % (x.id, x.formula))

    print('')
    print('Genes')
    print('-----')
    for x in model.genes:
        associated_ids = (i.id for i in x.reactions)
        print('%s is associated with reactions: %s' %
              (x.id, '{' + ', '.join(associated_ids) + '}'))

    model.objective = '3OAS140'

    print(model.objective.expression)
    print(model.objective.direction)

    class DemoFBA(object):
        def __init__(self, model):
            self.model = model

    return DemoFBA(model)

if __name__ == '__main__':
    metabolism = test_metabolism()
    # metabolism = test_minimal()
    # metabolism = test_canonical()
    # metabolism = test_demo()
    # metabolism = test_test()

    # cobra.io.save_json_model(metabolism.model, 'demo_model.json')

    print(metabolism.model)
    print(metabolism.model.reactions)
    print(metabolism.model.metabolites)
    print(metabolism.model.genes)
    print(metabolism.model.compartments)
    print(metabolism.model.solver)
    print(metabolism.model.objective.expression)

    print(metabolism.model.optimize().objective_value)
