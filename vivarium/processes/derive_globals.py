from __future__ import absolute_import, division, print_function

import copy

from scipy import constants
import numpy as np

from vivarium.actor.process import Process
from vivarium.utils.units import units



AVOGADRO = constants.N_A * 1 / units.mol



class DeriveGlobals(Process):
    """
    Process for deriving volume, mmol_to_counts factor, and growth rate
    from the cell mass
    """
    def __init__(self, initial_parameters={}):

        roles = {
            'global': [
                'mass',
                'volume',
                'growth_rate',
                'prior_mass',
                'mmol_to_counts',
                'density']}

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveGlobals, self).__init__(roles, parameters)

    def default_settings(self):
        # default state
        mass = 1339 * units.fg  # wet mass in fg
        density = 1100 * units.g / units.L
        volume = mass/density
        mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
        global_state = {
            'growth_rate': 0.0,
            'mass': mass.magnitude,
            'volume': volume.to('fL').magnitude,
            'mmol_to_counts': mmol_to_counts.magnitude,
            'prior_mass': mass.magnitude,
            'density': density.magnitude}

        default_state = {
            'global': global_state}

        # default emitter keys
        default_emitter_keys = {
            'global': ['volume', 'growth_rate']}

        # schema
        set_states = ['volume', 'growth_rate', 'prior_mass', 'mmol_to_counts']
        schema = {
            'global': {
                state_id : {
                    'updater': 'set'}
                for state_id in set_states}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):

        # states
        density = states['global']['density'] * units.g / units.L
        prior_mass = states['global']['prior_mass'] * units.fg
        mass = states['global']['mass'] * units.fg

        # update volume and growth rate
        volume =  mass / density
        mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
        growth_rate = (mass - prior_mass) / timestep / mass
        deriver_update = {
            'volume': volume.to('fL').magnitude,
            'mmol_to_counts': mmol_to_counts.magnitude,
            'growth_rate': growth_rate.magnitude,
            'prior_mass': mass.magnitude}

        return {
            'global': deriver_update}



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
    saved_data = test_deriver(100)
    print(saved_data)
