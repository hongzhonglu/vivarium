from __future__ import absolute_import, division, print_function

import cobra
import cobra.test

from scipy import constants
from cobra import Model, Reaction, Metabolite, Configuration

from lens.actor.process import Process, deep_merge
from lens.utils.units import units

EXTERNAL_SUFFIX = '_external'
REVERSE_SUFFIX = '_REVERSE'

def build_model(stoichiometry, objective, external_molecules):
    model = Model('fba')
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

class CobraFBA(object):
    cobra_configuration = Configuration()

    def __init__(self, config={}):
        self.nAvogadro = constants.N_A * 1/units.mol
        self.density = 1100 * units.g/units.L

        self.initial_state = config['initial_state']
        self.stoichiometry = config['stoichiometry']
        self.external_molecules = config['external_molecules']

        self.objective = config['objective']

        self.model = build_model(
            self.stoichiometry,
            self.objective,
            self.external_molecules)

        self.solution = None

    def set_external_levels(self, levels):
        for external, level in levels.items():
            reaction = self.model.reactions.get_by_id(external + EXTERNAL_SUFFIX)
            reaction.lower_bound = -level

    def objective_value(self):
        return self.solution.objective_value if self.solution else float('nan')

    def optimize(self):
        self.solution = self.model.optimize()
        return self.objective_value()

    def external_reactions(self):
        return [
            molecule + EXTERNAL_SUFFIX
            for molecule in self.external_molecules]

    def internal_reactions(self):
        all_reactions = set(self.reaction_ids())
        return all_reactions - set(self.external_reactions())

    def read_levels(self, molecules, scale=1.0):
        return {
            molecule: self.solution.fluxes[molecule] * scale
            for molecule in molecules}

    def read_internal_levels(self):
        return self.read_levels(self.internal_reactions())

    def read_external_levels(self, scale=1.0):
        external = self.external_reactions()
        levels = self.read_levels(external)
        return {
            molecule[:len(molecule) - len(EXTERNAL_SUFFIX)]: level
            for molecule, level in levels.items()}

    def reaction_ids(self):
        return [reaction.id for reaction in self.model.reactions]

    def get_reactions(self, reactions=[]):
        if not reactions:
            reactions = self.reaction_ids()

        return {
            reaction: self.model.reactions.get_by_id(reaction)
            for reaction in reactions}

    def set_reaction_bounds(self, reaction_bounds):
        reactions = self.get_reactions(reaction_bounds.keys())
        for reaction, bounds in reaction_bounds.items():
            reaction = reactions[reaction]
            reaction.lower_bound, reaction.upper_bound = bounds

    def get_reaction_bounds(self, reactions=[]):
        return {
            reaction_key: (reaction.lower_bound, reaction.upper_bound)
            for reaction_key, reaction in self.get_reactions(reactions).items()}
            

def test_minimal():
    stoichiometry = {
        'R1': {'A': -1, 'B': 1},
        'EB': {'B': -1}}

    objective = {'EB': 1.0}

    external_molecules = ['A']

    initial_state = {
        'A': 5}

    fba = CobraFBA({
        'stoichiometry': stoichiometry,
        'objective': objective,
        'external_molecules': external_molecules,
        'initial_state': initial_state})

    fba.set_external_levels(initial_state)

    return fba

def test_fba():
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
        'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1},
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'BIOMASS': 1}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2', 'BIOMASS']

    objective = {'v_biomass': 1.0}

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

    fba = CobraFBA({
        'stoichiometry': stoichiometry,
        'objective': objective,
        'external_molecules': external_molecules,
        'initial_state': initial_state})

    fba.set_external_levels(initial_state['external'])

    return fba

class JsonFBA(object):
    def __init__(self, path):
        self.model = cobra.io.load_json_model(path)

def test_canonical():
    fba = JsonFBA('models/e_coli_core.json')
    return fba

def test_test():
    fba = JsonFBA('models/minimal_model.json')
    return fba

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
    fba = test_fba()
    # fba = test_minimal()
    # fba = test_canonical()
    # fba = test_demo()
    # fba = test_test()

    # cobra.io.save_json_model(fba.model, 'demo_model.json')

    print(fba.model)
    print(fba.model.reactions)
    print(fba.model.metabolites)
    print(fba.model.genes)
    print(fba.model.compartments)
    print(fba.model.solver)
    print(fba.model.objective.expression)

    print(fba.optimize())
    print(fba.model.summary())
    print('internal: {}'.format(fba.internal_reactions()))
    print('external: {}'.format(fba.external_reactions()))
    print(fba.reaction_ids())
    print(fba.get_reactions())
    print(fba.get_reaction_bounds())
    print(fba.read_external_levels())
