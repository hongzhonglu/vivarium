from __future__ import absolute_import, division, print_function

import os

from lens.actor.process import deep_merge
from lens.data.spreadsheets import load_tsv
from lens.data.helper import get_mols_from_stoich, get_mols_from_reg_logic

from lens.environment.lattice_compartment import add_str_in_list, remove_str_in_list, add_str_to_keys
from lens.environment.make_media import Media
import lens.utils.regulation_logic as rl

from lens.processes.metabolism import Metabolism


DATA_DIR = os.path.join('lens', 'data', 'flat')
LIST_OF_FILENAMES = (
   "covert2002_reactions.tsv",
   "covert2002_transport.tsv",
   "covert2002_exchange_fluxes.tsv",
   "covert2002_maintenance_biomass_fluxes.tsv",
   "covert2002_GLC_G6P_flux_bounds.tsv",
   "covert2002_ecoli_metabolism_met_mw.tsv",
   )

def Covert2002Metabolism(parameters):
    '''load in flux_targets and target_key through parameters'''
    config = load_data(DATA_DIR, LIST_OF_FILENAMES)
    config.update(parameters)

    return Metabolism(config)

def load_data(data_dir, filenames):
    '''Load raw data from TSV files'''

    external_key = '[e]'

    data = {}
    for filename in filenames:
        attrName = filename.split(os.path.sep)[-1].split(".")[0]
        data[attrName] = load_tsv(data_dir, filename)

    ## Stoichiometry
    stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_reactions']}
    transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_transport']}
    maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_maintenance_biomass_fluxes']}
    exchange_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_exchange_fluxes']}
    stoichiometry.update(transport_stoichiometry)
    stoichiometry.update(maintenance_stoichiometry)
    stoichiometry.update(exchange_stoichiometry)

    # list of reversible reactions
    reversible_reactions = [reaction['Reaction'] for reaction in data['covert2002_reactions'] if reaction['Reversible']]

    # get all molecules
    metabolites = get_mols_from_stoich(stoichiometry)
    enzymes = [reaction['Protein'] for reaction in data['covert2002_reactions'] if reaction['Protein'] is not '']
    transporters = [reaction['Protein'] for reaction in data['covert2002_transport'] if reaction['Protein'] is not '']
    regulation_molecules = get_mols_from_reg_logic(data['covert2002_reactions'])

    all_molecules = set(metabolites + enzymes + transporters + regulation_molecules)
    all_molecules.remove('mass')
    # external_molecules = remove_str_in_list([mol_id for mol_id in all_molecules if external_key in mol_id], external_key)
    external_molecules = [mol_id for mol_id in all_molecules if external_key in mol_id]
    internal_molecules = [mol_id for mol_id in all_molecules if external_key not in mol_id]

    # molecular weights
    molecular_weights = {molecule['molecule id']: molecule['molecular weight']
         for molecule in data['covert2002_ecoli_metabolism_met_mw']}

    ## Objective
    objective = {'VGRO': 1}

    ## Flux bounds on reactions
    reaction_bounds = {flux['flux']: [flux['lower'], flux['upper']]
                   for flux in data['covert2002_GLC_G6P_flux_bounds']}
    default_reaction_bounds = reaction_bounds.pop('default')

    ## Regulatory logic functions
    regulation_logic = {}
    for reaction in data['covert2002_reactions']:
        reaction_id = reaction['Reaction']
        rule = rl.build_rule(reaction['Regulatory Logic'])
        if rule({}):
            regulation_logic[reaction_id] = rule

    ## Initial state
    # external
    make_media = Media()
    external_state = make_media.get_saved_media('GLC_G6P')
    external_state = add_str_to_keys(external_state, external_key)

    # # internal
    internal_state = {}
    mass_volume = {'mass': 1339,'volume': 1}
    mols = {mol_id: 0 for mol_id in internal_molecules}  # TODO -- are initial states for regulation molecules known?
    rxns = {rxn_id: 0.0 for rxn_id in stoichiometry.keys()}
    internal_state.update(mass_volume)
    internal_state.update(mols)
    internal_state.update(rxns)

    initial_state = {
        'external': external_state,
        'internal': internal_state}

    return {
        'stoichiometry': stoichiometry,
        'reversible': reversible_reactions,
        # 'reaction_bounds': reaction_bounds,  # TODO -- option for reaction_bounds in metabolism
        # 'default_reaction_bounds': default_reaction_bounds,
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state,
        'molecular_weights': molecular_weights,}


# tests and analyses of process
def test_covert2002():
    # configure process
    metabolism = Covert2002Metabolism({})

    print(metabolism.fba.model)
    print(metabolism.fba.model.reactions)
    print(metabolism.fba.model.metabolites)
    print(metabolism.fba.model.genes)
    print(metabolism.fba.model.compartments)
    print(metabolism.fba.model.solver)
    print(metabolism.fba.model.objective.expression)

    print(metabolism.fba.optimize())
    print(metabolism.fba.model.summary())
    print('internal: {}'.format(metabolism.fba.internal_reactions()))
    print('external: {}'.format(metabolism.fba.external_reactions()))
    print(metabolism.fba.reaction_ids())
    print(metabolism.fba.get_reactions())
    print(metabolism.fba.get_reaction_bounds())
    print(metabolism.fba.read_external_fluxes())


if __name__ == '__main__':
    from lens.processes.metabolism import simulate_metabolism, plot_output, save_network

    saved_state = test_covert2002()
    out_dir = os.path.join('out', 'tests', 'CovertPalsson2002_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## set up metabolism with a toy configuration
    # get_config = get_toy_configuration
    metabolism = Covert2002Metabolism({})

    ## test model
    test_covert2002()

    ## simulate model
    saved_data = simulate_metabolism(metabolism, 600)
    plot_output(saved_data, out_dir)

    ## make flux network from toy model
    save_network(metabolism, 10, out_dir)