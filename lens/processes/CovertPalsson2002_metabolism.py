from __future__ import absolute_import, division, print_function

import os
import re
from scipy import constants
import numpy as np

from lens.actor.process import Process, dict_merge
from lens.data.spreadsheets import load_tsv
from lens.data.helper import mols_from_reg_logic
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.environment.make_media import Media
import lens.utils.regulation_logic as rl


DATA_DIR = os.path.join('lens', 'data', 'flat')
LIST_OF_FILENAMES = (
    "covert2002_reactions.tsv",
    "covert2002_maintenance_biomass_fluxes.tsv",
    "covert2002_transport.tsv",
    "covert2002_exchange_fluxes.tsv",
    "covert2002_GLC_G6P_flux_bounds.tsv",
    )

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS


# helper functions
def get_reverse(reactions):
    reverse_stoichiometry = {}
    for reaction in reactions:
        if reaction['Reversible']:
            reaction_id = reaction['Reaction']
            stoich = {mol_id: -1 * coeff
                      for mol_id, coeff in reaction['Stoichiometry'].iteritems()}
            reverse_stoichiometry[reaction_id + '_reverse'] = stoich
    return reverse_stoichiometry

def get_molecules_from_stoich(stoichiometry):
    # TODO -- merge with rate_law_utilities.get_molecules_from_stoich
    molecules = set()
    for reaction, stoich in stoichiometry.iteritems():
        molecules.update(stoich.keys())
    return list(molecules)


class Metabolism(Process):
    def __init__(self, initial_parameters={}):
        self.e_key = '[e]'
        self.nAvogadro = constants.N_A * 1/units.mol

        # load data from files
        data = self.load_data()

        ## Initialize FBA
        objective = {"mass": 1.0}
        external_mol_ids_e = self.add_e_key(self.external_molecule_ids)

        self.fba = FluxBalanceAnalysis(
            reactionStoich=data['stoichiometry'],
            externalExchangedMolecules=external_mol_ids_e,
            objective=objective,
            objectiveType="standard",
            solver="glpk-linear",
        )

        # assign internal, external state ids
        roles = {
            'external': data['external_state_ids'],
            'internal': data['internal_state_ids']}
        parameters = {}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
            - external (dict) -- external states with default initial values, will be overwritten by environment
            - internal (dict) -- internal states with default initial values
        '''
        glc_g6p = True
        glc_lct = False

        make_media = Media()
        if glc_g6p:
            external = make_media.get_saved_media('GLC_G6P')
        elif glc_lct:
            external = make_media.get_saved_media('GLC_LCT')
        internal = {'mass': 0.032, 'volume': 1}

        # add reaction fluxes to internal state
        rxns = {rxn_id: 0.0 for rxn_id in self.reaction_ids}
        internal_state = dict_merge(internal, rxns)

        return {
            'external': external,
            'internal': internal_state}

    def default_emitter_keys(self):
        keys = {
            'internal': ['mass', 'lacI'] + self.transport_ids, #self.reaction_ids,
            'external': ['GLC', 'LAC', 'ACET']
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta'''

        updater_types = {
            'internal': {rxn_id: 'set' for rxn_id in self.reaction_ids},  # reactions set values directly
            'external': {mol_id: 'accumulate' for mol_id in self.external_molecule_ids}}  # all external values use default 'delta' udpater

        return updater_types

    def next_update(self, timestep, states):

        internal_state = states['internal']
        external_state = self.add_e_to_dict(states['external'])
        volume = internal_state['volume'] * units.fL

        # get the regulatory state of the reactions
        # TODO -- what to do about reversible regulated reactions?
        total_state = dict_merge(internal_state, external_state)
        boolean_state = {mol_id: (value>0) for mol_id, value in total_state.iteritems()}
        regulatory_state = {mol_id: regulatory_logic(boolean_state)
                            for mol_id, regulatory_logic in self.regulation_logic.iteritems()}

        # get exchange_molecule ids from FBA, remove self.e_key, look up in external_state
        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [external_state.get(molID, 0.0)
                                   for molID in exchange_molecules] * CONC_UNITS

        # set external_concentrations in FBA (mmol/L)
        self.fba.setExternalMoleculeLevels(external_concentrations.to(CONC_UNITS).magnitude)

        # set reaction flux bounds, based on default bounds and regulatory_state
        flux_bounds = np.array([self.default_flux_bounds for rxn in self.reaction_ids])
        for rxn_index, rxn_id in enumerate(self.reaction_ids):
            if regulatory_state.get(rxn_id) is False:
                flux_bounds[rxn_index] = [0.0, 0.0]
            # set lower bound to 0
            flux_bounds[rxn_index, 0] = 0 # TODO -- if negative, lower bounds should be set bound on reverse reactions

        self.fba.setReactionFluxBounds(
            self.reaction_ids,
            lowerBounds=flux_bounds[:,0],
            upperBounds=flux_bounds[:,1])

        # TODO -- should units have time in denominator?
        exchange_fluxes = CONC_UNITS * self.fba.getExternalExchangeFluxes() * timestep
        objective_value = CONC_UNITS * self.fba.getObjectiveValue() * timestep

        growth = {'mass': objective_value.to(CONC_UNITS).magnitude}

        # update state based on internal and external concentrations
        mmolToCount = self.nAvogadro.to('1/mmol') * volume.to('L')  # convert volume fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = (mmolToCount * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.external_molecule_ids, delta_exchange_counts))  # TODO -- refactor self.external_molecule_ids

        # get delta reaction flux
        rxn_ids = self.fba.getReactionIDs()
        rxn_fluxes = self.fba.getReactionFluxes()
        rxn_dict = dict(zip(rxn_ids, rxn_fluxes))

        update = {
            'internal': dict_merge(growth, rxn_dict),
            'external': environment_deltas,
        }

        return update

    def add_e_key(self, molecule_ids):
        return [mol_id + self.e_key for mol_id in molecule_ids]

    def remove_e_key(self, molecule_ids):
        return [mol_id.replace(self.e_key, '') for mol_id in molecule_ids]

    def add_e_to_dict(self, molecules_dict):
        ''' convert external state to compatible format by adding e_key'''
        e_dict = {}
        for key, value in molecules_dict.iteritems():
            if self.e_key in key:
                e_dict[key] = value
            else:
                e_dict[key + self.e_key] = value
        return e_dict

    def load_data(self):
        '''Load raw data from TSV files,
        save to data dictionary and then assign to class variables

        TODO -- what is covert2002_exchange_fluxes doing besides providing external molecule ids?
        TODO -- remove Growth reaction from covert2002_exchange_fluxes?
        '''

        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        # compose stoichiometry
        stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_reactions']}

        # additional stoichiomteries
        transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_transport']}
        maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_maintenance_biomass_fluxes']}
        reverse_stoichiometry = get_reverse(data['covert2002_reactions'])
        reverse_transport_stoichiometry = get_reverse(data['covert2002_transport'])

        # add to stoichiometry
        stoichiometry.update(transport_stoichiometry)
        stoichiometry.update(maintenance_stoichiometry)
        stoichiometry.update(reverse_stoichiometry)
        stoichiometry.update(reverse_transport_stoichiometry)

        # get all molecules
        metabolites = get_molecules_from_stoich(stoichiometry)
        enzymes = [reaction['Protein'] for reaction in data['covert2002_reactions'] if reaction['Protein'] is not '']
        transporters = [reaction['Protein'] for reaction in data['covert2002_transport'] if reaction['Protein'] is not '']
        regulation_molecules = mols_from_reg_logic(data['covert2002_reactions'])

        all_molecules = list(set(metabolites + enzymes + transporters + regulation_molecules))
        external_molecules = self.remove_e_key([mol_id for mol_id in all_molecules if self.e_key in mol_id])
        internal_molecules = [mol_id for mol_id in all_molecules if self.e_key not in mol_id]

        # reaction ids for tracking fluxes
        self.reaction_ids = stoichiometry.keys()
        self.transport_ids = transport_stoichiometry.keys()  # transport_ids are used by default_emitter

        # save external molecule ids, for use in update
        self.external_molecule_ids = external_molecules

        # make regulatory logic functions
        self.regulation_logic = {}
        for reaction in data['covert2002_reactions']:
            reaction_id = reaction['Reaction']
            rule = rl.build_rule(reaction['Regulatory Logic'])
            if rule({}):
                self.regulation_logic[reaction_id] = rule


        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.external_molecule_ids}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['covert2002_GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']

        return {
            'internal_state_ids': internal_molecules + self.reaction_ids + ['volume'],
            'external_state_ids': external_molecules,
            'stoichiometry': stoichiometry,
        }
