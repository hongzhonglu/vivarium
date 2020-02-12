from __future__ import absolute_import, division, print_function

import os
import random
import copy

from scipy import constants
import numpy as np

from vivarium.actor.process import Process
from vivarium.utils.units import units


AVOGADRO = constants.N_A * 1 / units.mol


class DeriveGlobals(Process):
    """
    Process for deriving volume and molecule concentrations
    from the cell state (mass, molecule counts)
    """
    def __init__(self, initial_parameters={}):

        roles = {
            'global': ['mass', 'volume', 'growth_rate', 'prior_mass']}

        parameters = {
            'density': 1100 * units.g / units.L}
        parameters.update(initial_parameters)

        super(DeriveGlobals, self).__init__(roles, parameters)

    def default_settings(self):
        # default state
        mass = 1339 * units.fg  # wet mass in fg
        density = self.parameters['density']
        volume = mass/density
        global_state = {
            'growth_rate': 0.0,
            'mass': mass.magnitude,
            'volume': volume.to('fL').magnitude,
            'prior_mass': mass.magnitude}

        default_state = {
            'global': global_state}

        # default emitter keys
        default_emitter_keys = {
            'global': ['volume', 'growth_rate']}

        # default updaters
        set_states = ['volume', 'growth_rate', 'prior_mass']
        default_updaters = {
            'global': {state_id: 'set' for state_id in set_states}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):

        # parameters
        density = self.parameters['density'] # units.g / units.L

        # states
        prior_mass = states['global']['prior_mass'] * units.fg
        mass = states['global']['mass'] * units.fg

        # update volume and growth rate
        volume =  mass / density
        growth_rate = (mass - prior_mass) / timestep / mass
        deriver_update = {
            'volume': volume.to('fL').magnitude,
            'growth_rate': growth_rate.magnitude,
            'prior_mass': mass.magnitude}

        # combine updates
        update = {
            'global': deriver_update}

        return update



def test_deriver(total_time=10):

    growth_rate = 6e-3

    # configure process
    deriver = DeriveGlobals({})

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

        # get update
        update = deriver.next_update(timestep, state)

        # set derived state
        state['global'].update(update['global'])

        # save state
        saved_state[time] = copy.deepcopy(state)


    return saved_state

if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'derive_global')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {}
    saved_data = test_deriver(100)
