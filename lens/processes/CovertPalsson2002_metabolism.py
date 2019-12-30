from __future__ import absolute_import, division, print_function

import os

from lens.processes.metabolism import Metabolism
from lens.environment.make_media import Media
from lens.utils.units import units
import lens.utils.regulation_logic as rl

from lens.data.spreadsheets import load_tsv
from lens.data.helper import get_mols_from_stoich, get_mols_from_reg_logic
from lens.environment.lattice_compartment import add_str_in_list, remove_str_in_list, add_str_to_keys


DATA_DIR = os.path.join('lens', 'data', 'flat')
LIST_OF_FILENAMES = (
   "covert2002_reactions.tsv",
   "covert2002_transport.tsv",
   "covert2002_exchange_fluxes.tsv",
   "covert2002_maintenance_biomass_fluxes.tsv",
   "covert2002_GLC_G6P_flux_bounds.tsv",
   "covert2002_ecoli_metabolism_met_mw.tsv",
   )

def CovertMetabolism(parameters):
    config = load_data(DATA_DIR, LIST_OF_FILENAMES)
    config.update(parameters)

    return Metabolism(config)


def load_data(data_dir, filenames):
    '''Load raw data from TSV files'''

    external_key = '[e]'
    rxn_key = '__RXN'

    data = {}
    for filename in filenames:
        attrName = filename.split(os.path.sep)[-1].split(".")[0]
        data[attrName] = load_tsv(data_dir + '/'+ filename)

    ## Stoichiometry
    stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_reactions']}
    transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_transport']}
    maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
         for reaction in data['covert2002_maintenance_biomass_fluxes']}
    stoichiometry.update(transport_stoichiometry)
    stoichiometry.update(maintenance_stoichiometry)  # TODO -- in standard FBA, this should only be needed in objective

    # list of reversible reactions
    reversible_reactions = [reaction['Reaction'] for reaction in data['covert2002_reactions'] if reaction['Reversible']]

    # add rxn_key to all entries of stoichiometry and transport_stoichiometry
    # this helps identify reactions in later analyses.
    stoichiometry = add_str_to_keys(stoichiometry, rxn_key)
    reversible_reactions = add_str_in_list(reversible_reactions, rxn_key)

    # get all molecules
    metabolites = get_mols_from_stoich(stoichiometry)
    enzymes = [reaction['Protein'] for reaction in data['covert2002_reactions'] if reaction['Protein'] is not '']
    transporters = [reaction['Protein'] for reaction in data['covert2002_transport'] if reaction['Protein'] is not '']
    regulation_molecules = get_mols_from_reg_logic(data['covert2002_reactions'])

    all_molecules = set(metabolites + enzymes + transporters + regulation_molecules)
    all_molecules.remove('mass')
    external_molecules = remove_str_in_list([mol_id for mol_id in all_molecules if external_key in mol_id], external_key)
    internal_molecules = [mol_id for mol_id in all_molecules if external_key not in mol_id]

    ## Objective
    # objective = maintenance_stoichiometry['VGRO']
    objective = {'mass': 1}  # TODO -- this is non-standard

    ## Flux bounds on reactions
    flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
                   for flux in data['covert2002_GLC_G6P_flux_bounds']}
    default_flux_bounds = flux_bounds['default']

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
        'reversible': stoichiometry.keys(),  #reversible_reactions,
        'external_molecules': external_molecules,  # external molecules are for lattice environment
        # 'external_key': external_key,
        'objective': objective,
        'regulation': regulation_logic,
        # 'transport_limits': transport_limits,
        'flux_bounds': flux_bounds,
        'default_upper_bound': default_flux_bounds,
        'initial_state': initial_state}




def test_metabolism():
    # configure process
    metabolism = CovertMetabolism({})

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
    from lens.processes.metabolism import simulate_metabolism, save_network
    from lens.actor.process import convert_to_timeseries, plot_simulation_output
    out_dir = os.path.join('out', 'tests', 'covert_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # add toy transport to config
    toy_config = {}

    # get ecoli core metabolism model
    covert_metabolism = CovertMetabolism(toy_config)

    # simulate model
    simulation_config = {
        'process': covert_metabolism,
        'total_time': 1000,
        'environment_volume': 5e-13}

    plot_settings = {}
        # 'skip_roles': ['exchange'],
        # 'overlay': {
        #     'reactions': 'flux_bounds'}}

    saved_data = simulate_metabolism(simulation_config)
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # make flux network from model
    save_network(covert_metabolism, 10, out_dir)
