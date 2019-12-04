from __future__ import absolute_import, division, print_function

import os
import cobra

from lens.processes.metabolism import Metabolism

DATA_FILE = os.path.join('lens', 'data', 'json_files', 'e_coli_core.json')


def EcoliCoreMetabolism(parameters):
    '''load in flux_targets and target_key through parameters'''
    model = cobra.io.load_json_model(DATA_FILE)

    reactions = model.reactions
    boundary = model.boundary
    objective_expression = model.objective.expression.args

    # get stoichiometry
    stoichiometry = {}
    for reaction in reactions:
        metabolites = reaction.metabolites
        stoichiometry[reaction.id] = {metabolite.id: coeff for metabolite, coeff in metabolites.items()}

    # get external molecules
    external_molecules = []
    exchange_bounds = {}
    for reaction in boundary:
        metabolites = reaction.metabolites.keys()
        assert len(metabolites) == 1
        metabolite_id = metabolites[0].id
        external_molecules.append(metabolite_id)

        # get exchange bounds
        exchange_bounds[metabolite_id] = list(reaction.bounds)

    # get objective
    objective = {}
    for expression in objective_expression:
        exp_str = str(expression)
        coeff, reaction_id = exp_str.split('*')

        # make sure reaction is in model
        try:
            reactions.get_by_id(reaction_id)
            objective[reaction_id] = float(coeff)
        except:
            pass


    config = {
        'stoichiometry': stoichiometry,
    #     'reversible': reversible,
        'external_molecules': external_molecules,
        'objective': objective,
    #     'initial_state': initial_state,
        'exchange_bounds': exchange_bounds,
    #     'default_upper_bound': default_reaction_bounds,
    #     'molecular_weights': molecular_weights,
        }

    return Metabolism(config)


# tests and analyses of process
def test_core_ecoli():
    # configure process
    metabolism = EcoliCoreMetabolism({})

    print('MODEL: {}'.format(metabolism.fba.model))
    print('REACTIONS: {}'.format(metabolism.fba.model.reactions))
    print('METABOLITES: {}'.format(metabolism.fba.model.metabolites))
    print('GENES: {}'.format(metabolism.fba.model.genes))
    print('COMPARTMENTS: {}'.format(metabolism.fba.model.compartments))
    print('SOLVER: {}'.format(metabolism.fba.model.solver))
    print('EXPRESSION: {}'.format(metabolism.fba.model.objective.expression))

    print(metabolism.fba.optimize())
    print(metabolism.fba.model.summary())
    print('internal: {}'.format(metabolism.fba.internal_reactions()))
    print('external: {}'.format(metabolism.fba.external_reactions()))
    print(metabolism.fba.reaction_ids())
    print(metabolism.fba.get_reactions())
    print(metabolism.fba.get_reaction_bounds())
    print(metabolism.fba.read_exchange_fluxes())





if __name__ == '__main__':
    ## test model
    test_core_ecoli()
