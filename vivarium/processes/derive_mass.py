from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Process
from vivarium.utils.units import units

from vivarium.processes.derive_globals import AVOGADRO

class DeriveMass(Process):
    """"""
    def __init__(self, initial_parameters={}):

        ports = initial_parameters['ports']
        ports.update({
            'global': ['mass', 'volume', 'mmol_to_counts']})

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
        # schema = {
        #     'global': {
        #         'mass': {
        #             'updater': 'set'}}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            # 'schema': schema
        }

    def next_update(self, timestep, states):
        # volume = states['global']['volume'] * units.fL
        # mmol_to_counts = states['global']['mmol_to_counts'] * units.L / units.mmol
        mass_schema = self.get_schema(self.in_ports, 'mass')

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

        # print('mass deriver mass: {}'.format(mass))
        # print('mass deriver states: {}'.format(states_with_mass['internal']))
        # import ipdb; ipdb.set_trace()
        # TODO -- deriver is run 2X for every metabolism update????




        # return {
        #     'global': {'mass': mass}}

        return {}
