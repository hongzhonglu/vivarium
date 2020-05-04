from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Process
from vivarium.utils.units import units

from vivarium.processes.derive_globals import AVOGADRO

class DeriveMass(Process):
    """
    Derives and sets total mass from individual molecular counts
    that have a mass schema in their stores .

    """

    defaults = {
        'dark_mass': 0,  # fg of dark mass
    }

    def __init__(self, initial_parameters={}):

        self.source_ports = initial_parameters.get('source_ports', {})
        target_ports = initial_parameters.get('target_ports', {'global': []})
        if target_ports:
            assert len(target_ports) == 1, 'DeriveMass too many target ports'
            assert list(target_ports.keys())[0] == 'global', 'DeriveMass requires target port named global'

        self.dark_mass = initial_parameters.get(
            'dark_mass', self.defaults['dark_mass']) * units.fg

        if self.dark_mass > 0:
            self.source_ports['global'] = ['dark_mass']

        # make the ports
        ports = self.source_ports.copy()

        # add global variables to global port
        global_variables = ['mass', 'volume', 'mmol_to_counts']
        if 'global' in ports:
            variable_set = set(ports['global'])
            variable_set.update(global_variables)
            ports['global'] = list(variable_set)
        else:
            ports['global'] = global_variables

        super(DeriveMass, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_state = {
            'global': {'dark_mass': (self.dark_mass.to('g') * AVOGADRO).magnitude}
        }

        # emitter keys
        default_emitter_keys = {
            'global': ['mass', 'dark_mass']}

        # schema
        schema = {
            'global': {
                'mass': {
                    'updater': 'set'},
                'dark_mass': {
                    'mass': 1  # 1 g/mol of dark mass
                }
            }
        }

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

    def next_update(self, timestep, states):
        mass_schema = self.schema_properties(self.source_ports, 'mass')

        # get mass from all molecules in source_ports
        mass = 0
        for port, molecules in self.source_ports.items():
            for mol_id in molecules:
                count = states[port][mol_id]
                mw = mass_schema[port].get(mol_id, 0.0) * (units.g / units.mol)
                mol = count / AVOGADRO
                added_mass = mw * mol
                mass += added_mass.to('fg')

        return {
            'global': {'mass': mass.magnitude}}
