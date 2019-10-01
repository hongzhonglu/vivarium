from __future__ import absolute_import, division, print_function

import os
import re
from scipy import constants

from lens.actor.process import Process
from lens.data.spreadsheets import load_tsv
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.environment.make_media import Media
from lens.utils.regulation_logic import build_rule


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

def merge_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


class Metabolism(Process):
    def __init__(self, initial_parameters={}):
        self.exchange_key = initial_parameters['exchange_key']
        self.e_key = '[e]'
        parameters = {'nAvogadro': constants.N_A}

        # load data from files
        self.load_data()

        # get internal molecules from stoichiometry
        all_molecule_ids = get_molecules_from_stoich(self.stoichiometry)

        # remove external molecules
        self.internal_molecule_ids = [mol_id
            for mol_id in all_molecule_ids if mol_id not in self.external_molecule_ids + ['mass']]

        # add internal regulation_molecules
        self.internal_molecule_ids = list(set(self.internal_molecule_ids) | self.regulation_molecules)

        # add e_key to external_molecules to match stoichiometry
        external_molecule_ids_e = [mol_id + self.e_key for mol_id in self.external_molecule_ids]

        # add reaction fluxes to internal states
        self.reaction_ids = self.stoichiometry.keys()  # [rxn_id + '_RXN' for rxn_id in self.stoichiometry.keys()]
        internal_state = self.internal_molecule_ids + self.reaction_ids + ['volume']

        # initialize FBA
        objective = {"mass": 1.0}
        self.fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=external_molecule_ids_e,
            objective=objective,
            objectiveType="standard",
            solver="glpk-linear",
        )

        roles = {
            'external': self.external_molecule_ids,
            'internal': internal_state}
        parameters.update(initial_parameters)

        super(Metabolism, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
            - environment_deltas (list) -- external molecule ids with added self.exchange_key string, for use to accumulate deltas in state
            - environment_ids (list) -- unmodified external molecule ids for use to accumulate deltas in state
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
        internal = {'mass': 0.032,
                    'volume': 1}

        environment_deltas = [key + self.exchange_key for key in external.keys()]

        # declare the states
        external_molecules = merge_dicts(external,{key: 0 for key in environment_deltas})

        # add reaction fluxes to internal state
        rxns = {rxn_id: 0.0 for rxn_id in self.reaction_ids}
        internal_state = merge_dicts(internal, rxns)

        return {
            'environment_deltas': environment_deltas,
            'external': external_molecules,
            'internal': internal_state}

    def default_emitter_keys(self):
        keys = {
            'internal': self.reaction_ids,
            'external': []
        }
        return keys

    def next_update(self, timestep, states):

        internal_state = states['internal']
        external_state = states['external']
        volume = internal_state['volume']
        # TODO -- set transport bounds

        # get exchange_molecule ids from FBA, remove self.e_key, and look up in external_state
        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [external_state.get(molID.replace(self.e_key, ''), 0.0)
            for molID in exchange_molecules]
        self.fba.setExternalMoleculeLevels(external_concentrations)

        # TODO -- check units on exchange flux
        exchange_fluxes = self.fba.getExternalExchangeFluxes() # TODO * timestep

        # update state based on internal and external concentrations
        countsToMolar = 1 / (self.parameters['nAvogadro'] * volume * 1e-15)  # convert volume fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.external_molecule_ids, delta_exchange_counts))

        # get delta reaction flux
        # TODO -- direct update flux, rather than passing delta
        rxn_ids = self.fba.getReactionIDs()
        rxn_fluxes = self.fba.getReactionFluxes()
        rxn_dict = dict(zip(rxn_ids, rxn_fluxes))
        rxn_delta = {rxn_id: new - internal_state[rxn_id] for rxn_id, new in rxn_dict.iteritems()}

        # TODO -- update internal state mass
        update = {
            'internal': rxn_delta,
            'external': {mol_id + self.exchange_key: delta
                for mol_id, delta in environment_deltas.iteritems()},
        }

        return update

    def load_data(self):
        # Load raw data from TSV files, save to data dictionary and then assign to class variables
        data = {}
        for filename in LIST_OF_FILENAMES:
            attrName = filename.split(os.path.sep)[-1].split(".")[0]
            data[attrName] = load_tsv(DATA_DIR, filename)

        # compose stoichiometry from files
        self.stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_reactions']}

        # get additional stoichiomteries
        transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_transport']}
        maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_maintenance_biomass_fluxes']}
        reverse_stoichiometry = get_reverse(data['covert2002_reactions'])
        reverse_transport_stoichiometry = get_reverse(data['covert2002_transport'])

        # add to stoichiometry
        self.stoichiometry.update(transport_stoichiometry)
        self.stoichiometry.update(maintenance_stoichiometry)
        self.stoichiometry.update(reverse_stoichiometry)
        self.stoichiometry.update(reverse_transport_stoichiometry)

        # make regulatory logic functions
        self.regulation_functions = {
            reaction['Reaction']: build_rule(reaction['Regulatory Logic'])
            for reaction in data['covert2002_reactions']}

        # get all molecules listed in "Regulatory Logic"
        self.regulation_molecules = set()
        regex_split = 'IF |not |or |and |\(|\)| '
        for reaction in data['covert2002_reactions']:
            reg_logic = reaction['Regulatory Logic']
            reg_logic_parsed = re.split(regex_split, reg_logic)  # split string to remove patterns in regex_split
            reg_logic_parsed = list(filter(None, reg_logic_parsed))  # remove empty strings
            self.regulation_molecules.update(reg_logic_parsed)

        # remove external molecules from regulation_molecules and remove the self.e_key
        external_regulation_molecules = set([mol_id
            for mol_id in self.regulation_molecules if self.e_key in mol_id])
        self.regulation_molecules.difference_update(external_regulation_molecules)
        external_regulation_molecules = set([mol_id.replace(self.e_key, '')
            for mol_id in external_regulation_molecules])

        # get list of external molecules.
        # TODO -- what is covert2002_exchange_fluxes doing besides providing external molecule ids?
        # TODO -- remove Growth reaction from covert2002_exchange_fluxes?
        self.external_molecule_ids = [reaction['Stoichiometry'].keys()[0].replace(self.e_key, '')
            for reaction in data['covert2002_exchange_fluxes'] if reaction['Reaction'] != 'Growth']

        # check that all external regulation molecules are in external_molecule_ids
        set_ex = set(self.external_molecule_ids)
        in_reg = external_regulation_molecules.difference(set_ex)
        assert len(in_reg) == 0

        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.external_molecule_ids}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['covert2002_GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']
