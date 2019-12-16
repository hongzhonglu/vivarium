from __future__ import absolute_import, division, print_function

import cobra
import cobra.test

from cobra import Model, Reaction, Metabolite, Configuration


EXTERNAL_SUFFIX = '_exchange'

def build_model(stoichiometry, reversible, objective, external_molecules, default_upper_bound=1000):
    model = Model('fba')
    model.compartments = {'c': 'cytoplasm'}

    metabolite_keys = {}
    for reaction_key, chemistry in stoichiometry.items():
        metabolite_keys.update(chemistry)

    metabolites = {
        metabolite: Metabolite(metabolite, name=metabolite, compartment='c')
        for metabolite in metabolite_keys.keys()}

    model.add_metabolites(metabolites.values())

    # make reactions
    reactions = {}
    for reaction_key, chemistry in stoichiometry.items():
        reaction = Reaction(reaction_key, name=reaction_key)

        # set reaction bounds
        reaction.upper_bound = default_upper_bound
        if reaction_key in reversible:
            reaction.lower_bound = -reaction.upper_bound

        # make stoichiometry
        reaction_model = {
            metabolites[metabolite]: value
            for metabolite, value in chemistry.items()}
        reaction.add_metabolites(reaction_model)

        reactions[reaction_key] = reaction

    # make exchange reactions for all external_molecules
    for external in external_molecules:
        external_key = external + EXTERNAL_SUFFIX
        reaction = Reaction(external_key, name=external_key)

        # set reaction bounds
        reaction.upper_bound = default_upper_bound
        reaction.lower_bound = -default_upper_bound  # TODO -- should exchanges have symmetric bounds by default?

        # make stoichiometry
        reaction_model = {metabolites[external]: -1}
        reaction.add_metabolites(reaction_model)

        reactions[external_key] = reaction

    model.add_reactions(reactions.values())

    model.objective = {
        reactions[reaction_key]: value
        for reaction_key, value in objective.items()}

    return model

def extract_model(model):
    reactions = model.reactions
    metabolites = model.metabolites
    boundary = model.boundary
    objective_expression = model.objective.expression.args

    # get stoichiometry
    stoichiometry = {}
    flux_bounds = {}
    for reaction in reactions:
        reaction_metabolites = reaction.metabolites
        stoichiometry[reaction.id] = {
            metabolite.id: coeff for metabolite, coeff in reaction_metabolites.items()}
        # get flux bounds
        flux_bounds[reaction.id] = list(reaction.bounds)

    # get external molecules
    external_molecules = []
    external_state = {}
    exchange_bounds = {}
    for reaction in boundary:
        reaction_metabolites = reaction.metabolites.keys()
        assert len(reaction_metabolites) == 1  # only 1 molecule in the exchange reaction
        metabolite_id = reaction_metabolites[0].id
        external_molecules.append(metabolite_id)
        external_state[metabolite_id] = 0.0

        # get exchange bounds
        exchange_bounds[metabolite_id] = list(reaction.bounds)

    # get molecular weights
    molecular_weights = {}
    for metabolite in metabolites:
        molecular_weights[metabolite.id] = metabolite.formula_weight

    # get objective
    objective = {}
    for expression in objective_expression:
        exp_str = str(expression)
        coeff, reaction_id = exp_str.split('*')
        try:
            reactions.get_by_id(reaction_id)
            objective[reaction_id] = float(coeff)
        except:
            pass

    # initial state
    initial_state = {
        'internal': {
            'mass': 1339,  # fg
            'volume': 1E-15},  # fL
        'external': external_state}

    # adjustments
    for exchange_id, [lb, ub] in exchange_bounds.items():
        exchange_bounds[exchange_id] = [lb/100, ub/100]

    config = {
        'stoichiometry': stoichiometry,
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state,
        'flux_bounds': flux_bounds,
        'exchange_bounds': exchange_bounds,
        'molecular_weights': molecular_weights}

    return config
    

class CobraFBA(object):
    cobra_configuration = Configuration()

    def __init__(self, config={}):
        if config.get('model_path'):
            self.model = cobra.io.load_json_model(model_path)
            extract = extract_model(self.model)

            self.stoichiometry = extract['stoichiometry']
            self.reversible = extract.get('reversible', [])
            self.external_molecules = extract['external_molecules']
            self.objective = extract['objective']
            self.flux_bounds = extract['flux_bounds']
            self.exchange_bounds = extract['exchange_bounds']
            self.molecular_weights = extract['molecular_weights']
        else:
            self.stoichiometry = config['stoichiometry']
            self.reversible = config.get('reversible', [])
            self.external_molecules = config['external_molecules']
            self.objective = config['objective']
            self.flux_bounds = config.get('flux_bounds', {})
            self.exchange_bounds = config.get('exchange_bounds', {})
            self.molecular_weights = config.get('molecular_weights', {})

            default_upper_bound = config.get('default_upper_bound', 1000.0)

            self.model = build_model(
                self.stoichiometry,
                self.reversible,
                self.objective,
                self.external_molecules,
                default_upper_bound)

            self.constrain_reaction_bounds(self.flux_bounds)
            self.constrain_exchange_flux(self.exchange_bounds)

        self.solution = None

    def constrain_exchange_flux(self, levels):
        for external, level in levels.items():
            reaction = self.model.reactions.get_by_id(external + EXTERNAL_SUFFIX)

            if type(level) is list:
                reaction.upper_bound = -level[0]
                reaction.lower_bound = -level[1]

            elif isinstance(level, int) or isinstance(level, float):
                reaction.lower_bound = -level

    def constrain_flux(self, levels):
        for external, level in levels.items():
            reaction = self.model.reactions.get_by_id(external)
            reaction.upper_bound = level

    def constrain_reaction_bounds(self, reaction_bounds):
        reactions = self.get_reactions(reaction_bounds.keys())
        for reaction, bounds in reaction_bounds.items():
            reaction = reactions[reaction]
            reaction.lower_bound, reaction.upper_bound = bounds

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

    def read_fluxes(self, molecules):
        return {
            molecule: self.solution.fluxes[molecule]
            for molecule in molecules}

    def read_internal_fluxes(self):
        return self.read_fluxes(self.internal_reactions())

    def read_exchange_fluxes(self):
        external = self.external_reactions()
        levels = self.read_fluxes(external)
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
        'reversible': stoichiometry.keys(),
        'objective': objective,
        'external_molecules': external_molecules,
        'initial_state': initial_state,
        })

    fba.constrain_exchange_flux(initial_state)

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
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'BIOMASS': 1}
    }

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
        'reversible': stoichiometry.keys(),
        'objective': objective,
        'external_molecules': external_molecules})

    fba.constrain_exchange_flux(initial_state['external'])

    return fba

class JsonFBA(object):
    def __init__(self, path):
        self.model = cobra.io.load_json_model(path)

def test_canonical():
    fba = JsonFBA('lens/data/json_files/e_coli_core.json')
    return fba

# def test_test():
#     fba = JsonFBA('models/minimal_model.json')
#     return fba

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
    print('MODEL: {}'.format(fba.model))
    print('REACTIONS: {}'.format(fba.model.reactions))
    print('METABOLITES: {}'.format(fba.model.metabolites))
    print('GENES: {}'.format(fba.model.genes))
    print('COMPARTMENTS: {}'.format(fba.model.compartments))
    print('SOLVER: {}'.format(fba.model.solver))
    print('EXPRESSION: {}'.format(fba.model.objective.expression))

    print(fba.optimize())
    print(fba.model.summary())
    print('internal: {}'.format(fba.internal_reactions()))
    print('external: {}'.format(fba.external_reactions()))
    print(fba.reaction_ids())
    print(fba.get_reactions())
    print(fba.get_reaction_bounds())
    print(fba.read_exchange_fluxes())
