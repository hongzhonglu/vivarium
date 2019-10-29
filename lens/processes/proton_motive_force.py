from __future__ import absolute_import, division, print_function

from lens.actor.process import Process
from lens.utils.units import units


class ProtonMotiveForce(Process):
    '''
    Need to add a boot method for this process to lens/environment/boot.py for it to run on its own
    '''
    def __init__(self, initial_parameters={}):

        parameters = {
            'F': 96485.33289, # (charge / mol)  Faraday constant
        }


        roles = {
            'internal': ['Q_in'],
            'membrane': ['V_m'],
            'external': [],
        }
        parameters.update(initial_parameters)

        super(ProtonMotiveForce, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
        default_state = {
            'external': states (dict) -- external states ids with default initial values
            'internal': states (dict) -- internal states ids with default initial values
        '''

        return {
            'internal': {},
            'membrane': {},
            'external': {},
        }

    def default_emitter_keys(self):
        '''
        returns dictionary with:
        keys = {
            'internal': states (list), # a list of states to emit from internal
            'external': states (list), # a list of states to emit from external
        }
        '''
        keys = {
            'internal': [],
            'membrane': [],
            'external': [],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta,
        which is accumulated and passed to the environment at every exchange step'''
        keys = {
            'membrane': {
                'V_m': 'accumulate'}}
        return keys

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        Q_in = internal_state['Q_in']

        # V_m = electric potential difference
        # Q_in is the intracellular charge (in mole)
        # C is the membrane capacitance
        # F is the Faraday constant
        V_m = self.parameters['F'] * Q_in / C


        update = {
            # 'internal': {},
            'membrane': {'V_m': V_m},
            # 'external': {},
        }
        return update