from __future__ import absolute_import, division, print_function

import os

from lens.actor.process import dict_merge
from lens.data.spreadsheets import load_tsv
from lens.data.helper import mols_from_reg_logic

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
   rxn_key = '__RXN'

   data = {}
   for filename in filenames:
       attrName = filename.split(os.path.sep)[-1].split(".")[0]
       data[attrName] = load_tsv(data_dir, filename)

   # compose stoichiometry
   stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
                    for reaction in data['covert2002_reactions']}
   transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
                              for reaction in data['covert2002_transport']}
   maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
                                for reaction in data['covert2002_maintenance_biomass_fluxes']}
   stoichiometry.update(transport_stoichiometry)

   # objective
   # objective = maintenance_stoichiometry['VGRO']
   # stoichiometry.update({'ATPM': maintenance_stoichiometry['ATPM']})
   stoichiometry.update(maintenance_stoichiometry)
   objective = {'mass': 1}

   # list of reversible reactions
   reversible_reactions = [reaction['Reaction'] for reaction in data['covert2002_reactions'] if reaction['Reversible']]

   # add rxn_key to all entries of stoichiometry and transport_stoichiometry
   stoichiometry = add_str_to_keys(stoichiometry, rxn_key)
   reversible_reactions = add_str_in_list(reversible_reactions, rxn_key)

   # get all molecules
   metabolites = get_molecules_from_stoich(stoichiometry)
   enzymes = [reaction['Protein'] for reaction in data['covert2002_reactions'] if reaction['Protein'] is not '']
   transporters = [reaction['Protein'] for reaction in data['covert2002_transport'] if reaction['Protein'] is not '']
   regulation_molecules = mols_from_reg_logic(data['covert2002_reactions'])

   all_molecules = list(set(metabolites + enzymes + transporters + regulation_molecules))
   external_molecules = remove_str_in_list([mol_id for mol_id in all_molecules if external_key in mol_id], external_key)
   internal_molecules = [mol_id for mol_id in all_molecules if external_key not in mol_id]

   # flux bounds on reactions
   flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
                  for flux in data['covert2002_GLC_G6P_flux_bounds']}
   # default_flux_bounds = flux_bounds['default']

   # # make regulatory logic functions from data
   # regulation_logic = {}
   # for reaction in data['covert2002_reactions']:
   #     reaction_id = reaction['Reaction']
   #     rule = rl.build_rule(reaction['Regulatory Logic'])
   #     if rule({}):
   #         regulation_logic[reaction_id] = rule


   # initial state
   make_media = Media()
   external_state = make_media.get_saved_media('GLC_G6P')
   external_state = add_str_to_keys(external_state, external_key)
   internal = {'mass': 1339,'volume': 1}
   rxns = {rxn_id: 0.0 for rxn_id in stoichiometry.keys()}
   internal_state = dict_merge(dict(internal), rxns)
   initial_state = {'external': external_state, 'internal': internal_state}

   # TODO -- add regulation_logic, flux_bounds
   return {
       'stoichiometry': stoichiometry,
       'reversible_reactions': stoichiometry.keys(),  #reversible_reactions,
       'external_molecules': external_molecules,  # external molecules are for lattice environment
       'external_key': external_key,
       'objective': objective,
       # 'transport_limits': transport_limits,
       'initial_state': initial_state}

def get_molecules_from_stoich(stoichiometry):
   molecules = set()
   for reaction, stoich in stoichiometry.iteritems():
       molecules.update(stoich.keys())
   return list(molecules)



def test_covert2002(total_time=3600):
    import numpy as np
    from lens.utils.units import units

    target_key = '__target'
    # make kinetic rate laws to mimic transport kinetics
    transport_rates = {}

    # configure process
    metabolism = Covert2002Metabolism({})
    target_rxn_ids = metabolism.flux_targets

    # get initial state and parameters
    state = metabolism.default_state()
    density = metabolism.density
    nAvogadro = metabolism.nAvogadro

    saved_state = {
        'internal': {state_id: [value] for state_id, value in state['internal'].iteritems()},
        'external': {state_id: [value] for state_id, value in state['external'].iteritems()},
        'time': [0]}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < total_time:
        mass_t0 = state['internal']['mass']
        volume_t0 = (mass_t0 * units.fg).to('g') / density
        time += timestep

        # get update
        update = metabolism.next_update(timestep, state)
        state['internal'] = update['internal']
        growth_rate = metabolism.fba.getObjectiveValue()

        # apply external update
        mmolToCount = (nAvogadro.to('1/mmol') * volume_t0).to('L/mmol').magnitude
        for mol_id, exchange in update['external'].iteritems():
            exchange_rate = exchange / mmolToCount  # TODO -- per second?
            delta_conc = exchange_rate / growth_rate * mass_t0 * (np.exp(growth_rate * timestep) - 1)

            if np.isnan(delta_conc):
                delta_conc = 0.0

            state['external'][mol_id] += delta_conc * 0.001  # TODO -- scaling?
            if state['external'][mol_id] < 0:  # this shouldn't be needed
                state['external'][mol_id] = 0

        # # get new flux targets
        # target_fluxes = {}
        # for rxn_id, rate_law in transport_rates.iteritems():
        #     target_flux = rate_law(state['external'])
        #     target_fluxes[rxn_id + target_key] = target_flux
        # state['internal'].update(target_fluxes)

        # save state
        saved_state['time'].append(time)
        saved_state['internal']['volume'].append(volume_t0.magnitude)  # TODO -- get new volume
        for state_id, value in state['internal'].iteritems():
            saved_state['internal'][state_id].append(value)
        for state_id, value in state['external'].iteritems():
            saved_state['external'][state_id].append(value)

    # check that mass values are strictly increasing
    mass_ts = saved_state['internal']['mass'][2:]
    assert all(i < j for i, j in zip(mass_ts, mass_ts[1:]))

    data = {
        'saved_state': saved_state,
        'target_rxn_ids': target_rxn_ids,
        'target_key': target_key}
    return data

def plot_environment_output(data):
    # plot glc/mass only
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    saved_state = data['saved_state']
    time_vec = [t/60./60. for t in saved_state['time']]  # convert to hours
    external_data = data['saved_state']['external']
    mass_series = data['saved_state']['internal']['mass']


    # make figure
    fig, ax1 = plt.subplots(figsize=(5,5))
    color = 'black'
    ax1.plot(time_vec, mass_series, color=color, linewidth=3.0, label='mass')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_xlabel('time (hr)')
    ax1.set_ylabel('mass (fg)', color=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    for state_id, series in external_data.iteritems():
        ax2.plot(time_vec, series, label=state_id)
    ax2.set_ylabel('env concentrations (mmol/L)', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.legend(prop={'size': 6}, loc=7)

    # make figure output directory and save figure
    fig_path = os.path.join('out', 'covert2002_metabolism_environment')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')

def plot_metabolism_output(data):
    import os
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import math

    saved_state = data['saved_state']
    target_rxn_ids = data.get('target_rxn_ids')
    target_key = data.get('target_key')
    targeted_states = remove_str_in_list(target_rxn_ids, target_key)

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [t/60./60. for t in saved_state['time']]  # convert to hours
    n_data = [len(saved_state[key].keys()) for key in data_keys]

    n_col = 3
    n_rows = math.ceil(sum(n_data) / n_col) + 1

    # make figure
    fig = plt.figure(figsize=(n_col * 8, n_rows * 2.5))
    plot_idx = 1
    for key in data_keys:
        for state_id, series in sorted(saved_state[key].iteritems()):
            if state_id not in target_rxn_ids:
                ax = fig.add_subplot(n_rows, n_col, plot_idx)
                ax.plot(time_vec, series)
                ax.title.set_text(str(key) + ': ' + state_id)
                ax.ticklabel_format(style='sci',axis='y')
                ax.set_xlabel('time (hr)')
                plot_idx += 1

                # if state has target, plot series in red
                if state_id in targeted_states:
                    target_id = state_id + target_key
                    target_series = saved_state[key][target_id]
                    ax.plot(time_vec, target_series, 'r', label='target')
                    ax.legend()

    # make figure output directory and save figure
    fig_path = os.path.join('out', 'covert2002_metabolism_all')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')

def save_metabolic_network():
    # TODO -- make this function into an analysis
    import math
    from lens.utils.make_network import make_network, save_network

    flux_simstep = 10 # the simulation step used to set edge weights with reaction flux

    # initialize process
    metabolism = Covert2002Metabolism({})
    stoichiometry = metabolism.stoichiometry
    reaction_ids = stoichiometry.keys()
    external_mol_ids = metabolism.fba.getExternalMoleculeIDs()
    objective = metabolism.objective

    # run test to get simulation output
    run_time = 60
    data = test_covert2002(run_time)
    saved_state = data['saved_state']
    external = saved_state['external']
    internal = saved_state['internal']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        flux = internal[rxn_id][flux_simstep]
        reaction_fluxes[rxn_id] = math.log(1000 * flux + 1.1)

    # transport node type
    node_types = {rxn_id: 'reaction' for rxn_id in reaction_ids}
    node_types.update({mol_id: 'external_mol' for mol_id in external_mol_ids})
    node_types.update({mol_id: 'objective' for mol_id in objective.keys()})
    info = {
        'node_types': node_types,
        'reaction_fluxes': reaction_fluxes}

    nodes, edges = make_network(stoichiometry, info)
    save_network(nodes, edges, 'out/covert2002_metabolism_network')

if __name__ == '__main__':
   saved_state = test_covert2002()
   plot_metabolism_output(saved_state)
   plot_environment_output(saved_state)
   save_metabolic_network()
