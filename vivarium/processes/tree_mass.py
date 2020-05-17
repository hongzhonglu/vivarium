from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Deriver
from vivarium.utils.units import units

from vivarium.processes.derive_globals import AVOGADRO

class TreeMass(Deriver):
    """
    Derives and sets total mass from individual molecular counts
    that have a mass schema in their stores .

    """

    def __init__(self, initial_parameters={}):
        ports = {
            'global': [
                'mass']}

        super(DeriveMass, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_state = {
            'global': {
                'mass': 0}}

        # emitter keys
        default_emitter_keys = {
            'global': ['mass', 'dark_mass']}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

    def ports_schema(self):
        return {
            'global': {
                'mass': {
                    '_default': 0,
                    '_updater': 'set'}}}

    def next_update(self, timestep, states):
        def calculate_mass(value, path, node):
            if 'mass' in node.properties:
                unit_mass = node.properties['mass']
                count = node.value
                mw = unit_mass * (units.g / units.mol)
                mol = count / AVOGADRO
                added_mass = mw * mol
                mass = added_mass.to('fg')
                return value + mass.magnitude
            else:
                return value

        return {
            'global': {
                'mass': {
                    '_reduce': {
                        'reducer': calculate_mass,
                        'from': ('..',),
                        'initial': 0.0}}}}
