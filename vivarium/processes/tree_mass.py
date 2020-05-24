from __future__ import absolute_import, division, print_function

from vivarium.core.process import Deriver
from vivarium.utils.units import units

from vivarium.processes.derive_globals import AVOGADRO


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


class TreeMass(Deriver):
    """
    Derives and sets total mass from individual molecular counts
    that have a mass schema in their stores .

    """

    defaults = {
        'from_path': ('..', '..')}

    def __init__(self, initial_parameters={}):
        self.from_path = self.or_default(initial_parameters, 'from_path')

        ports = {
            'global': [
                'mass']}

        super(TreeMass, self).__init__(ports, initial_parameters)

    def ports_schema(self):
        return {
            'global': {
                'mass': {
                    '_default': 0,
                    '_updater': 'set'}}}

    def next_update(self, timestep, states):
        return {
            'global': {
                'mass': {
                    '_reduce': {
                        'reducer': calculate_mass,
                        'from': self.from_path,
                        'initial': 0.0}}}}
