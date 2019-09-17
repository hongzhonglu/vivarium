from __future__ import absolute_import, division, print_function

import os
import csv
from scipy import constants

from lens.actor.process import Process
from lens.environment.make_media import Media
from lens.environment.look_up import LookUp
from lens.utils.rate_law_utilities import load_reactions
from lens.utils.rate_law_utilities import get_reactions_from_exchange
from lens.utils.rate_law_utilities import get_molecules_from_reactions
from lens.data.spreadsheets import JsonReader
from itertools import ifilter


EXTERNAL_MOLECULES_FILE = os.path.join('lens', 'data', 'flat', 'environment_molecules.tsv')
TRANSPORT_IDS_FILE = os.path.join('lens', 'data', 'flat', 'transport_reactions.tsv')

TSV_DIALECT = csv.excel_tab

amino_acids = [
    'L-ALPHA-ALANINE',
    'ARG',
    'ASN',
    'L-ASPARTATE',
    'CYS',
    'GLT',
    'GLN',
    'GLY',
    'HIS',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'L-SELENOCYSTEINE',
    'VAL'
]
additional_exchange = ['OXYGEN-MOLECULE', 'GLC']
external_molecule_ids = additional_exchange + amino_acids

# add [p] label. TODO (Eran) -- fix this
external_molecule_ids_p = [mol_id + '[p]' for mol_id in external_molecule_ids]


class TransportLookup(Process):
    def __init__(self, initial_parameters={}):
        self.exchange_key = initial_parameters['exchange_key']
        self.media_id = 'minimal' # initial_parameters.get('media_id', 'minimal')
        self.lookup_type = 'average' # initial_parameters.get('lookup', 'average')
        self.nAvogadro = constants.N_A
        self.external_molecule_ids = external_molecule_ids

        # load all reactions and maps
        self.load_data()

        # external_molecule_ids declares which molecules' exchange will be applied
        self.transport_reaction_ids = get_reactions_from_exchange(self.all_transport_reactions, external_molecule_ids_p)
        all_molecule_ids = get_molecules_from_reactions(self.transport_reaction_ids, self.all_transport_reactions)
        internal_molecule_ids = [mol_id for mol_id in all_molecule_ids if mol_id not in external_molecule_ids_p]

        # make look up object
        self.look_up = LookUp()

        # get initial fluxes
        self.transport_fluxes = self.look_up.look_up(
            self.lookup_type,
            self.media_id,
            self.transport_reaction_ids)

        roles = {
            'external': external_molecule_ids,
            'internal': internal_molecule_ids + ['volume']}
        parameters = {}
        parameters.update(initial_parameters)

        super(TransportLookup, self).__init__(roles, parameters)

    def default_state(self):
        media_id = 'minimal_plus_amino_acids'
        make_media = Media()
        media = make_media.get_saved_media(media_id)

        environment_deltas = [key + self.exchange_key for key in media.keys()]  # TODO (Eran) -- pass in '_change' string

        # declare the states
        environment_state = media
        environment_state.update({key: 0 for key in environment_deltas})
        environment_state['volume'] = 10

        cell_state = {'volume': 1}

        return {
            'environment_deltas': environment_deltas,
            'external': environment_state,
            'internal': cell_state}

    def default_emitter_keys(self):
        keys = {
            'internal': [],
            'external': external_molecule_ids
        }
        return keys

    def next_update(self, timestep, states):

        volume = states['internal']['volume']

        # nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
        self.molar_to_counts = (self.nAvogadro * 1e-3) * (volume * 1e-15)

        # get transport fluxes
        self.transport_fluxes = self.look_up.look_up(
            self.lookup_type,
            self.media_id,
            self.transport_reaction_ids)

        # convert to counts
        delta_counts = self.flux_to_counts(self.transport_fluxes)  # TODO * timestep

        # Get the deltas for environmental molecules
        environment_deltas = {}
        for molecule_id in delta_counts.keys():
            if molecule_id in self.molecule_to_external_map:
                external_molecule_id = self.molecule_to_external_map[molecule_id]
                environment_deltas[external_molecule_id] = delta_counts[molecule_id]

        update = {
            'external': {mol_id + self.exchange_key: delta for mol_id, delta in environment_deltas.iteritems()}, # TODO (Eran) -- pass in this '_change' substring
        }

        return update

    # TODO (Eran) -- make this a util
    def flux_to_counts(self, fluxes):
        rxn_counts = {
            reaction_id: int(self.molar_to_counts * flux)
            for reaction_id, flux in fluxes.iteritems()}
        delta_counts = {}
        for reaction_id, rxn_count in rxn_counts.iteritems():
            stoichiometry = self.all_transport_reactions[reaction_id]['stoichiometry']
            substrate_counts = {
                substrate_id: coeff * rxn_count
                for substrate_id, coeff in stoichiometry.iteritems()}
            # add to delta_counts
            for substrate, delta in substrate_counts.iteritems():
                if substrate in delta_counts:
                    delta_counts[substrate] += delta
                else:
                    delta_counts[substrate] = delta
        return delta_counts

    def load_data(self):
        '''
        - Loads all reactions, including locations for enzymes.
        - Separates out the transport reactions as an class dictionary
        - Makes mappings from molecule ids with location tags to external molecules without location tags

        '''

        # use rate_law_utilities to get all_reactions
        all_reactions = load_reactions()

        # make dict of reactions in TRANSPORT_IDS_FILE
        self.all_transport_reactions = {}
        with open(TRANSPORT_IDS_FILE, 'rU') as tsvfile:
            reader = JsonReader(
                ifilter(lambda x: x.lstrip()[0] != '#', tsvfile), # Strip comments
                dialect = TSV_DIALECT)
            for row in reader:
                reaction_id = row['reaction id']
                stoichiometry = all_reactions[reaction_id]['stoichiometry']
                reversible = all_reactions[reaction_id]['is reversible']
                transporters_loc = all_reactions[reaction_id]['catalyzed by']

                self.all_transport_reactions[reaction_id] = {
                    'stoichiometry': stoichiometry,
                    'is reversible': reversible,
                    'catalyzed by': transporters_loc,
                }

        # Make map of external molecule_ids with a location tag (as used in reaction stoichiometry) to molecule_ids in the environment
        self.molecule_to_external_map = {}
        self.external_to_molecule_map = {}
        with open(EXTERNAL_MOLECULES_FILE, 'rU') as tsvfile:
            reader = JsonReader(
                ifilter(lambda x: x.lstrip()[0] != '#', tsvfile), # Strip comments
                dialect = TSV_DIALECT)
            for row in reader:
                molecule_id = row['molecule id']
                location = row['exchange molecule location']
                self.molecule_to_external_map[molecule_id + location] = molecule_id
                self.external_to_molecule_map[molecule_id] = molecule_id + location
