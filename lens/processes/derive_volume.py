from __future__ import absolute_import, division, print_function

from lens.actor.process import Process

class DeriveVolume(Process):
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['mass', 'density', 'volume']}
        parameters = {}

        super(DeriveVolume, self).__init__(roles, parameters, deriver=True)

    def default_state(self):

        mass = 1.339e-12  # g  # 1339 (wet mass in fg)
        # dry_mass = 403.0 * units.fg
        density = 1100  # g/L
        volume = mass/density * 1e15  # fL

        default_state = {
            'mass': mass,
            'density': density,
            'volume': volume,
        }

        return {
            'internal': default_state
        }

    def default_emitter_keys(self):
        keys = {
            'internal': ['mass', 'density', 'volume'],
        }
        return keys

    def next_update(self, timestep, states):
        mass = states['internal']['mass']
        density = states['internal']['density']

        volume =  mass / density * 1e15  # fL
        update = {
            'internal': {'volume': volume}}

        return update