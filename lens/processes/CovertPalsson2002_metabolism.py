from __future__ import absolute_import, division, print_function

import os
from scipy import constants

from lens.actor.process import Process
from lens.data.spreadsheets import load_tsv
from lens.utils.units import units
from lens.utils.modular_fba import FluxBalanceAnalysis
from lens.environment.make_media import Media


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

def get_molecules_from_reactions(stoichiometry):
    # TODO -- merge with rate_law_utilities.get_molecules_from_reactions
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

        # parameters
        parameters = {
            'nAvogadro': constants.N_A,
            }

        self.load_data()
        all_molecule_ids = get_molecules_from_reactions(self.stoichiometry)
        self.internal_molecule_ids = [mol_id
            for mol_id in all_molecule_ids if mol_id not in self.external_molecule_ids + ['mass']]

        # initialize FBA
        self.fba = FluxBalanceAnalysis(
            reactionStoich=self.stoichiometry,
            externalExchangedMolecules=self.external_molecule_ids,
            objective=self.objective,
            objectiveType="standard",
            solver="glpk-linear",
        )

        roles = {
            'external': self.external_molecule_ids,
            'internal': self.internal_molecule_ids + ['volume']}
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
        internal = {'mass': 0.032}

        environment_deltas = [key + self.exchange_key for key in external.keys()]

        # declare the states
        external_molecules = merge_dicts(external,{key: 0 for key in environment_deltas})
        internal_molecules = merge_dicts(internal, {'volume': 1})  # fL


        import ipdb; ipdb.set_trace()

        return {
            'environment_deltas': environment_deltas,
            'external': external_molecules,
            'internal': internal_molecules}

    def next_update(self, timestep, states):

        volume = states['internal']['volume']
        external_state = states['external']

        # TODO -- set transport bounds
        exchange_molecules = self.fba.getExternalMoleculeIDs()
        external_concentrations = [external_state.get(molID, 0.0) for molID in exchange_molecules]
        self.fba.setExternalMoleculeLevels(external_concentrations)

        # TODO -- check units on exchange flux
        exchange_fluxes = self.fba.getExternalExchangeFluxes() # TODO * timestep

        # update state based on internal and external concentrations
        countsToMolar = 1 / (self.parameters['nAvogadro'] * volume * 1e-15)  # convert volume fL to L

        # Get the delta counts for environmental molecules
        delta_exchange_counts = ((1 / countsToMolar) * exchange_fluxes).astype(int)
        environment_deltas = dict(zip(self.external_molecule_ids, delta_exchange_counts))

        # TODO -- update internal state mass
        update = {
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

        self.stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_reactions']}
        transport_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_transport']}
        maintenance_stoichiometry = {reaction['Reaction']: reaction['Stoichiometry']
            for reaction in data['covert2002_maintenance_biomass_fluxes']}
        self.stoichiometry.update(transport_stoichiometry)
        self.stoichiometry.update(maintenance_stoichiometry)

        reverse_stoichiometry = get_reverse(data['covert2002_reactions'])
        reverse_transport_stoichiometry = get_reverse(data['covert2002_transport'])
        self.stoichiometry.update(reverse_stoichiometry)
        self.stoichiometry.update(reverse_transport_stoichiometry)

        # TODO -- remove growth exchange flux from file?
        self.external_molecule_ids = [reaction['Stoichiometry'].keys()[0]
            for reaction in data['covert2002_exchange_fluxes'] if reaction['Reaction'] != 'Growth']

        self.objective = {"mass": 1.0}

        self.transport_limits = {mol_id: 1.0 * (units.mmol / units.g / units.h)
            for mol_id in self.external_molecule_ids}

        flux_bounds = {flux['flux']: [flux['lower'], flux['upper']]
            for flux in data['covert2002_GLC_G6P_flux_bounds']}
        self.default_flux_bounds = flux_bounds['default']

if __name__ == '__main__':
    Metabolism()
