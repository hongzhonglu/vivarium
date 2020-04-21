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
        'dark_mass': 0,
    }

    def __init__(self, initial_parameters={}):

        self.dark_mass = initial_parameters.get('dark_mass', self.defaults['dark_mass'])
        # TODO -- add dark mass to schema?
        import ipdb; ipdb.set_trace()

        source_ports = initial_parameters['source_ports']
        target_ports = initial_parameters['target_ports']

        if target_ports:
            assert len(target_ports) == 1, 'DeriveMass too many target ports'
            assert list(target_ports.keys())[0] == 'global', 'DeriveMass requires target port named global'

        ports = {
            'global': ['mass', 'volume', 'mmol_to_counts']}
        ports.update(source_ports)

        self.in_ports = {
            port_id: keys
            for port_id, keys in ports.items()
            if port_id not in ['global']}

        super(DeriveMass, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_state = {}

        # emitter keys
        default_emitter_keys = {'global': ['mass']}

        # schema
        schema = {
            'global': {
                'mass': {
                    'updater': 'set'}}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

    def next_update(self, timestep, states):
        mass_schema = self.schema_properties(self.in_ports, 'mass')

        states_with_mass = {
            port: molecules
            for port, molecules in states.items()
            if port in self.in_ports}

        mass = 0
        for port, molecules in states_with_mass.items():
            for mol_id, count in molecules.items():
                mw = mass_schema[port].get(mol_id, 0.0) * (units.g / units.mol)
                mol = count / AVOGADRO
                added_mass = mw * mol
                mass += added_mass.to('fg')

        return {
            'global': {'mass': mass.magnitude}}
