from __future__ import absolute_import, division, print_function

from vivarium.actor.process import Process
from vivarium.utils.units import units


class DeriveVolume(Process):
    def __init__(self, initial_parameters={}):
        parameters = {'density': 1100}  # * units.g / units.L
        roles = {
            'internal': ['mass', 'volume'],
        }
        parameters.update(initial_parameters)

        super(DeriveVolume, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        mass = 1339 * units.fg  # 1.339e-12 g  # 1339 (wet mass in fg)
        density = self.parameters['density'] * units.g / units.L
        volume = mass/density
        internal = {
            'mass': mass.magnitude,
            'volume': volume.to('fL').magnitude}
        default_state = {'internal': internal}

        # default emitter keys
        default_emitter_keys = {'internal': ['mass', 'volume']}

        # default updaters
        default_updaters = {'internal': {'volume': 'set'}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):
        mass = states['internal']['mass'] * units.fg
        density = self.parameters['density'] * units.g / units.L
        volume =  mass / density
        update = {'internal': {'volume': volume.to('fL').magnitude}}

        return update