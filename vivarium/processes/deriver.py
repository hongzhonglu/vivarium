from __future__ import absolute_import, division, print_function

import os
import random
import copy

import numpy as np
from scipy import constants

from vivarium.actor.process import Process, convert_to_timeseries, plot_simulation_output
from vivarium.utils.units import units

AVOGADRO = constants.N_A * 1 / units.mol
GLOBALS = ['mass', 'volume', 'growth_rate', 'prior_mass']

class Deriver(Process):
    """
    Process for deriving volume and molecule concentrations
    from the cell state (mass, molecule counts)
    """
    def __init__(self, initial_parameters={}):

        # get counted molecules
        counted_molecules = initial_parameters.get('counted_molecules', [])

        roles = {
            'state': counted_molecules,
            'global': GLOBALS,
            'counts': counted_molecules}

        parameters = {
            'density': 1100 * units.g / units.L,
            'nAvogadro': constants.N_A * 1 / units.mol}
        parameters.update(initial_parameters)

        super(Deriver, self).__init__(roles, parameters)

    def default_settings(self):
        # default state
        mass = 1339 * units.fg  # wet mass in fg
        density = self.parameters['density']
        volume = mass/density
        global_state = {
            'growth_rate': 0.0,
            'mass': mass.magnitude,
            'prior_mass': mass.magnitude,
            'volume': volume.to('fL').magnitude}

        molecule_ids = self.roles['counts']
        counted_molecules = {mol_id: 0 for mol_id in molecule_ids}

        default_state = {
            'state': counted_molecules,
            'global': global_state,
            'counts': counted_molecules}

        # default emitter keys
        default_emitter_keys = {
            'state': molecule_ids,
            'global': ['volume', 'growth_rate']}

        # default updaters
        default_updaters = {
            'state': {state_id: 'set' for state_id in molecule_ids},
            'global': {state_id: 'set' for state_id in GLOBALS}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):

        # parameters
        density = self.parameters['density'] # units.g / units.L
        nAvogadro = self.parameters['nAvogadro']

        # states
        prior_mass = states['global']['prior_mass'] * units.fg
        mass = states['global']['mass'] * units.fg
        counts = states['counts']

        # update volume and growth rate
        volume =  mass / density
        growth_rate = (mass - prior_mass) / timestep / mass
        deriver_update = {
            'volume': volume.to('fL').magnitude,
            'growth_rate': growth_rate.magnitude,
            'prior_mass': mass.magnitude}

        # concentration update
        mmol_to_count = nAvogadro.to('1/mmol') * volume.to('L')
        concentration_update = {mol_id: (count / mmol_to_count).magnitude
            for mol_id, count in counts.items()}

        for mol_id, conc in concentration_update.items():
            assert conc >= 0, 'derived {} concentration < 0'.format(mol_id)

        return {
            'global': deriver_update,
            'state': concentration_update}



def test_deriver(total_time=10):

    expression_rates = {'protein': 0.05}
    growth_rate = 6e-3

    # configure process
    initial_parameters = {
        'counted_molecules': list(expression_rates.keys())}
    deriver = Deriver(initial_parameters)

    # get initial state and parameters
    settings = deriver.default_settings()
    state = settings['state']

    # initialize saved data
    saved_state = {}

    ## simulation
    timeline = [(total_time, {})]
    time = 0
    timestep = 1  # sec
    saved_state[time] = state
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.items():
                    state[key].update(change)

        # set mass, counts
        mass = state['global']['mass']
        state['global']['mass'] = mass * np.exp(growth_rate * timestep)
        for mol_id, rate in expression_rates.items():
            if random.random() < rate:
                state['counts'][mol_id] += 1

        # get update
        update = deriver.next_update(timestep, state)

        # set derived state
        state['state'].update(update['state'])

        # save state
        saved_state[time] = copy.deepcopy(state)

        state['global']['prior_mass'] = update['global']['prior_mass']

    return saved_state

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'global_derive')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}
    saved_data = test_deriver(100)
    del saved_data[0] # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
