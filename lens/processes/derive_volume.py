from __future__ import absolute_import, division, print_function

from lens.actor.process import Process
from lens.utils.units import units


class DeriveVolume(Process):
    def __init__(self, initial_parameters={}):
        roles = {
            'internal': ['mass', 'density', 'volume'],
        }
        parameters = {}

        super(DeriveVolume, self).__init__(roles, parameters, deriver=True)

    def default_state(self):

        mass = 1339 * units.fg  # 1.339e-12 g  # 1339 (wet mass in fg)
        density = 1100  * units.g/units.L
        volume = mass/density

        default_state = {
            'mass': mass.magnitude,
            'density': density.magnitude,
            'volume': volume.to('fL').magnitude,
        }

        return {
            'internal': default_state,
        }

    def default_emitter_keys(self):
        keys = {
            'internal': ['mass', 'density', 'volume'],
        }
        return keys

    def next_update(self, timestep, states):
        mass = states['internal']['mass'] * units.fg
        density = states['internal']['density'] * units.g/units.L
        volume =  mass / density
        update = {'internal': {'volume': volume.to('fL').magnitude}}

        return update