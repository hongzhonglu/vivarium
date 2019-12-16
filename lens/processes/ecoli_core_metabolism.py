from __future__ import absolute_import, division, print_function

import os
import cobra

from lens.environment.make_media import Media

from lens.processes.metabolism import Metabolism

DATA_FILE = os.path.join('models', 'e_coli_core.json')
# DATA_FILE = os.path.join('lens', 'data', 'json_files', 'e_coli_core.json')


def EcoliCoreMetabolism(parameters):
    # '''load in flux_targets and target_key through parameters'''
    # model = cobra.io.load_json_model(DATA_FILE)

    # reactions = model.reactions
    # metabolites = model.metabolites
    # boundary = model.boundary
    # objective_expression = model.objective.expression.args

    # # get stoichiometry
    # stoichiometry = {}
    # flux_bounds = {}
    # for reaction in reactions:
    #     reaction_metabolites = reaction.metabolites
    #     stoichiometry[reaction.id] = {
    #         metabolite.id: coeff for metabolite, coeff in reaction_metabolites.items()}
    #     # get flux bounds
    #     flux_bounds[reaction.id] = list(reaction.bounds)

    # # get external molecules
    # external_molecules = []
    # exchange_bounds = {}
    # for reaction in boundary:
    #     reaction_metabolites = reaction.metabolites.keys()
    #     assert len(reaction_metabolites) == 1  # only 1 molecule in the exchange reaction
    #     metabolite_id = reaction_metabolites[0].id
    #     external_molecules.append(metabolite_id)

    #     # get exchange bounds
    #     exchange_bounds[metabolite_id] = list(reaction.bounds)

    # # get molecular weights
    # molecular_weights = {}
    # for metabolite in metabolites:
    #     molecular_weights[metabolite.id] = metabolite.formula_weight

    # # get objective
    # objective = {}
    # for expression in objective_expression:
    #     exp_str = str(expression)
    #     coeff, reaction_id = exp_str.split('*')
    #     try:
    #         reactions.get_by_id(reaction_id)
    #         objective[reaction_id] = float(coeff)
    #     except:
    #         pass

    # # initial state
    # make_media = Media()
    # ecoli_core_GLC_media = make_media.get_saved_media('ecoli_core_GLC')
    # initial_state = {
    #     'internal': {
    #         'mass': 1339,  # fg
    #         'volume': 1E-15},  # fL
    #     'external': ecoli_core_GLC_media,
    #     }

    # # adjustments
    # for exchange_id, [lb, ub] in exchange_bounds.items():
    #     exchange_bounds[exchange_id] = [lb/100, ub/100]

    # config = {
    #     'stoichiometry': stoichiometry,
    #     'external_molecules': external_molecules,
    #     'objective': objective,
    #     'initial_state': initial_state,
    #     'exchange_bounds': exchange_bounds,
    #     'flux_bounds': flux_bounds,
    # #     'default_upper_bound': default_reaction_bounds,
    #     'molecular_weights': molecular_weights,
    #     }

    config = {'model_path': parameters.get('model_path', DATA_FILE)}

    return Metabolism(config)

def kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def toy_transport_kinetics():
    transport_kinetics = {
        "GLCpts": kinetic_rate('glc__D_e', 1e1, 5),  # glucose
        "GLUt2r": kinetic_rate('glu__L_e', 1e1, 5),  # glucose
        # "PYRt2": kinetic_rate('pyr_e', 1e2, 5),
    }
    return transport_kinetics

def test_ecoli_core():
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
    from lens.processes.metabolism import simulate_metabolism, plot_output, save_network

    out_dir = os.path.join('out', 'tests', 'e_coli_core_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## test model
    test_ecoli_core()

    ## set up metabolism with a toy configuration
    metabolism = EcoliCoreMetabolism({})

    ## simulate model
    simulation_config = {
        'total_time': 100,
        'transport_kinetics': toy_transport_kinetics(),
        'environment_volume': 1e-13}
    saved_data = simulate_metabolism(metabolism, simulation_config)
    plot_output(saved_data, out_dir)

    ## make flux network from toy model
    save_network(metabolism, 10, out_dir)
