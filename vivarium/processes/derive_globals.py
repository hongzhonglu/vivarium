from __future__ import absolute_import, division, print_function

import copy
import math

from scipy import constants
import numpy as np

from vivarium.compartment.process import Process
from vivarium.utils.units import units
from vivarium.utils.dict_utils import deep_merge


PI = math.pi
AVOGADRO = constants.N_A * 1 / units.mol



def length_from_volume(volume, width):
    '''
    get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
    a=length of cylinder without rounded caps, l=total length:

    V = (4/3)*PI*r^3 + PI*r^2*a
    l = a + 2*r
    '''
    radius = width / 2
    cylinder_length = (volume - (4/3) * PI * radius**3) / (PI * radius**2)
    total_length = cylinder_length + 2 * radius
    return total_length

def volume_from_length(length, width):
    '''
    inverse of length_from_volume
    '''
    radius = width / 2
    cylinder_length = length - width
    volume = cylinder_length * (PI * radius**2) + (4 / 3) * PI * radius**3
    return volume

def surface_area_from_length(length, width):
    '''
    SA = 3*PI*r^2 + 2*PI*r*a
    '''
    radius = width / 2
    cylinder_length = length - width
    surface_area = 3 * PI * radius**2 + 2 * PI * radius * cylinder_length
    return surface_area



class DeriveGlobals(Process):
    """
    Process for deriving volume, mmol_to_counts, and shape from the cell mass
    """

    defaults = {
        'width': 1,  # um
    }

    def __init__(self, initial_parameters={}):

        self.width = initial_parameters.get('width', self.defaults['width'])
        source_ports = initial_parameters.get('source_ports')
        target_ports = initial_parameters.get('target_ports')

        if source_ports:
            assert len(source_ports) == 1, 'DeriveGlobals too many source ports'
            assert list(source_ports.keys())[0] == 'global', 'DeriveGlobals requires source port named global'
        if target_ports:
            assert len(target_ports) == 1, 'DeriveGlobals too many target ports'
            assert list(target_ports.keys())[0] == 'global', 'DeriveGlobals requires target port named global'

        ports = {
            'global': [
                'mass',
                'volume',
                'mmol_to_counts',
                'density',
                'width',
                'length',
                'surface_area']}

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveGlobals, self).__init__(ports, parameters)

    def default_settings(self):
        # default state
        mass = 1339 * units.fg  # wet mass in fg
        density = 1100 * units.g / units.L
        volume = mass/density
        mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
        length = length_from_volume(volume.magnitude, self.width)
        surface_area = surface_area_from_length(length, self.width)

        global_state = {
            'mass': mass.magnitude,
            'volume': volume.to('fL').magnitude,
            'mmol_to_counts': mmol_to_counts.magnitude,
            'density': density.magnitude,
            'width': self.width,
            'length': length,
            'surface_area': surface_area,
        }

        default_state = {
            'global': global_state}

        # default emitter keys
        default_emitter_keys = {
            'global': ['volume', 'width', 'length', 'surface_area']}

        # schema
        set_states = ['volume', 'mmol_to_counts', 'length', 'surface_area']
        set_divide = ['density']
        schema = {
            'global': {
                state_id : {
                    'updater': 'set'}
                for state_id in set_states}}
        divide_schema = {
            'global': {
                state_id : {
                    'divide': 'set'}
                for state_id in set_divide}}
        schema = deep_merge(schema, divide_schema)

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        # states
        density = states['global']['density'] * units.g / units.L
        mass = states['global']['mass'] * units.fg

        # get volume from mass, and more variables from volume
        volume = mass / density
        mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
        length = length_from_volume(volume.magnitude, self.width)
        surface_area = surface_area_from_length(length, self.width)

        return {
            'global': {
                'volume': volume.to('fL').magnitude,
                'mmol_to_counts': mmol_to_counts.magnitude,
                'length': length,
                'surface_area': surface_area}}



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
