"""
Lattice

A two-dimensional lattice environmental model

## Physics
# Diffusion constant of glucose in 0.5 and 1.5 percent agarose gel = ~6 * 10^-10 m^2/s (Weng et al. 2005. Transport of glucose and poly(ethylene glycol)s in agarose gels).
# Conversion to micrometers: 6 * 10^-10 m^2/s = 600 micrometers^2/s.

## Cell biophysics
# rotational diffusion in liquid medium with viscosity = 1 mPa.s: Dr = 3.5+/-0.3 rad^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
# translational diffusion in liquid medium with viscosity = 1 mPa.s: Dt=100 micrometers^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)

"""

from __future__ import absolute_import, division, print_function

import os
import math
import random
import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from vivarium.actor.outer import EnvironmentSimulation
from vivarium.utils.multicell_physics import MultiCellPhysics
from vivarium.environment.make_media import Media

# Constants
N_AVOGADRO = constants.N_A
PI = math.pi

CELL_DENSITY = 1100

# Lattice parameters
TRANSLATION_JITTER = 0.1
ROTATION_JITTER = 0.05

DEFAULT_TIMESTEP = 1.0

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


def gaussian(deviation, distance):
    return np.exp(-np.power(distance, 2.) / (2 * np.power(deviation, 2.)))

class EnvironmentSpatialLattice(EnvironmentSimulation):
    def __init__(self, config):

        ## Simulation Parameters
        self._time = 0
        self._max_time = 10e6
        self.output_dir = config.get('output_dir')
        self.run_for = config.get('run_for', DEFAULT_TIMESTEP)  # simulation time step

        ## Cell Parameters
        self.cell_density = CELL_DENSITY  # TODO -- get density from cell sim, or just get mass
        self.cell_radius = config.get('cell_radius', 0.5)  # TODO -- this should be a property of cells.
        self.translation_jitter = config.get('translation_jitter', TRANSLATION_JITTER)
        self.rotation_jitter = config.get('rotation_jitter', ROTATION_JITTER)
        self.cell_placement = config.get('cell_placement')

        ## Lattice Parameters
        # if edge_length_y not provided, assume it is equal to edge_length_x.
        # adjust provided edge_length_y to make dx=dy for easier diffusion.
        self.edge_length_x = config.get('edge_length_x', 10.0)
        edge_length_y = config.get('edge_length_y', self.edge_length_x)
        self.patches_per_edge_x = config.get('patches_per_edge_x', 10)
        self.patches_per_edge_y = int(edge_length_y * self.patches_per_edge_x / self.edge_length_x)
        self.edge_length_y = self.patches_per_edge_y * self.edge_length_x / self.patches_per_edge_x
        self.depth = config.get('depth', 3000.0)  # um
        self.total_volume = (self.depth * self.edge_length_x * self.edge_length_y) * (10 ** -15) # (L)
        self.patch_volume = self.total_volume / (self.edge_length_x * self.edge_length_y)
        self.dx = self.edge_length_x / self.patches_per_edge_x
        self.dy = self.edge_length_y / self.patches_per_edge_y
        self.dx2 = self.dx * self.dy

        ## Media and Timeline Parameters
        self.make_media = Media()
        self.media_id = config.get('media_id', 'minimal')
        self.timeline = ()
        timeline_str = config.get('timeline_str')
        new_media = config.get('new_media', {})
        get_concentrations = config.get('concentrations',{})
        self.end_timeline = False

        # add a new media
        for new_media_id, media_dict in new_media.items():
            media_id = self.make_media.add_media(media_dict, new_media_id)

        # make the timeline and media
        if timeline_str:
            self.timeline = self.make_media.make_timeline(timeline_str)
            self._times = [t[0] for t in self.timeline]
            self.media_id = self.timeline[0][1]
            media = self.make_media.get_saved_media(self.media_id)
        elif get_concentrations:
            # make media and fill lattice patches with media concentrations
            media = get_concentrations
        else:
            media = self.make_media.get_saved_media(self.media_id)

        self._molecule_ids = list(media.keys())
        self.concentrations = media.values()
        self.molecule_index = {molecule: index for index, molecule in enumerate(self._molecule_ids)}
        self.fill_lattice(media)

        ## Concentration and Gradient Parameters
        self.static_concentrations = config.get('static_concentrations', False)
        self.diffusion = config.get('diffusion', 1e3)
        self.gradient = config.get('gradient', {'type': False})

        # upper limit on the time scale. (using 50% of the theoretical upper limit)
        if self.diffusion:
            self.diffusion_dt = 0.5 * self.dx**2 * self.dy**2 / (2 * self.diffusion * (self.dx**2 + self.dy**2))
            if self.diffusion_dt > self.run_for:
                # if large diffusion_dt, it still needs to be within self.run_for
                self.diffusion_dt = self.run_for
        else:
            self.diffusion_dt = self.run_for

        # add a gradient
        if self.gradient.get('type') == 'gaussian':
            # gaussian gradient multiplies the basal concentration of the given molecule
            # by a gaussian function of distance from center and deviation
            #
            # 'gradient': {
            #     'type': 'gradient',
            #     'molecules': {
            #         'mol_id1':{
            #             'center': [0.25, 0.5],
            #             'deviation': 30},
            #         'mol_id2': {
            #             'center': [0.75, 0.5],
            #             'deviation': 30}
            #     }},

            for molecule_id, specs in self.gradient['molecules'].items():
                mol_index = self._molecule_ids.index(molecule_id)
                center = [specs['center'][0] * self.edge_length_x,
                          specs['center'][1] * self.edge_length_y]
                deviation = specs['deviation']

                for x_patch in range(self.patches_per_edge_x):
                    for y_patch in range(self.patches_per_edge_y):
                        # distance from middle of patch to center coordinates
                        dx = (x_patch + 0.5) * self.edge_length_x / self.patches_per_edge_x - center[0]
                        dy = (y_patch + 0.5) * self.edge_length_y / self.patches_per_edge_y - center[1]
                        distance = np.sqrt(dx ** 2 + dy ** 2)
                        scale = gaussian(deviation, distance)
                        # multiply gradient by scale
                        self.lattice[mol_index][x_patch][y_patch] *= scale

        elif self.gradient.get('type') == 'linear':
            # linear gradient adds to the basal concentration of the given molecule
            # as a function of distance from center and slope.
            #
            # 'gradient': {
            #     'type': 'linear',
            #     'molecules': {
            #         'mol_id1':{
            #             'center': [0.0, 0.0],
            #             'slope': -10},
            #         'mol_id2': {
            #             'center': [1.0, 1.0],
            #             'slope': -5}
            #     }},

            for molecule_id, specs in self.gradient['molecules'].items():
                mol_index = self._molecule_ids.index(molecule_id)
                center = [specs['center'][0] * self.edge_length_x,
                          specs['center'][1] * self.edge_length_y]
                slope = specs['slope']

                for x_patch in range(self.patches_per_edge_x):
                    for y_patch in range(self.patches_per_edge_y):
                        # distance from middle of patch to center coordinates
                        dx = (x_patch + 0.5) * self.edge_length_x / self.patches_per_edge_x - center[0]
                        dy = (y_patch + 0.5) * self.edge_length_y / self.patches_per_edge_y - center[1]
                        distance = np.sqrt(dx ** 2 + dy ** 2)
                        added = distance * slope
                        # add gradient to basal concentration
                        self.lattice[mol_index][x_patch][y_patch] += added

                self.lattice[mol_index][self.lattice[mol_index] <= 0.0] = 0.0

        ## Initialize dictionaries
        self.simulations = {}       # map of agent_id to simulation state
        self.locations = {}         # map of agent_id to center location and orientation
        self.corner_locations = {}  # map of agent_id to corner location, for Lens visualization and multi-cell physics engine
        self.motile_forces = {}	    # map of agent_id to motile force, with magnitude and relative orientation

        # make physics object by passing in bounds and jitter
        bounds = [self.edge_length_x, self.edge_length_y]
        self.multicell_physics = MultiCellPhysics(
            bounds,
            self.translation_jitter,
            self.rotation_jitter)

        # configure emitter and emit lattice configuration
        self.emitter = config['emitter'].get('object')
        self.emit_fields = config.get('emit_fields', [])
        self.emit_indices = [self._molecule_ids.index(mol_id) for mol_id in self.emit_fields]
        self.emit_configuration()


    def evolve(self, timestep):
        ''' Evolve environment '''
        self.update_locations(timestep)
        self.update_media()

        if not self.static_concentrations:
            self.run_diffusion(timestep)

        # make sure all patches have concentrations of 0 or higher
        self.lattice[self.lattice < 0.0] = 0.0

        # emit current state
        self.emit_data()

    def update_locations(self, timestep):
        ''' Update the location of all agents '''

        # update all agents in multicell physics
        for agent_id, location in self.locations.items():

            # shape
            agent_state = self.simulations[agent_id]['state']
            volume = agent_state['volume']
            radius = self.cell_radius
            width = 2 * radius
            length = self.volume_to_length(volume, radius)
            mass = agent_state.get('mass', 1.0)  # TODO -- pass mass through message.

            # update length, width in multicell_physics
            agent_state['length'] = length
            agent_state['width'] = width
            self.multicell_physics.update_cell(agent_id, length, width, mass)

            # update motile forces in multicell_physics
            force = self.motile_forces[agent_id][0]
            torque = self.motile_forces[agent_id][1]
            self.multicell_physics.apply_motile_force(agent_id, force, torque)

        # run multicell physics
        self.multicell_physics.run_incremental(timestep)

        # set new agent location
        for agent_id, location in self.locations.items():
            # set location
            self.locations[agent_id] = self.multicell_physics.get_center(agent_id)
            self.corner_locations[agent_id] = self.multicell_physics.get_corner(agent_id)

            # enforce boundaries # TODO (Eran) -- make pymunk handle this
            self.locations[agent_id][0:2][self.locations[agent_id][0:2] < 0] = 0.0
            if self.locations[agent_id][0] > self.edge_length_x:
                self.locations[agent_id][0] = self.edge_length_x - self.dx / 2
            if self.locations[agent_id][1] > self.edge_length_y:
                self.locations[agent_id][1] = self.edge_length_y - self.dy / 2

    def add_cell_to_physics(self, agent_id, position, angle):
        ''' Add body to multi-cell physics simulation'''

        volume = self.simulations[agent_id]['state']['volume']
        width = self.cell_radius * 2
        length = self.volume_to_length(volume, self.cell_radius)
        mass = volume * self.cell_density  # TODO -- get units to work

        # add length, width to state, for use by visualization
        self.simulations[agent_id]['state']['length'] = length
        self.simulations[agent_id]['state']['width'] = width

        self.multicell_physics.add_cell_from_center(
            cell_id = agent_id,
            width = width,
            length = length,
            mass = mass,
            center_position = position,
            angle = angle,
        )

        # add to lattice
        self.locations[agent_id] = self.multicell_physics.get_center(agent_id)
        self.corner_locations[agent_id] = self.multicell_physics.get_corner(agent_id)

    def update_media(self):
        if self.timeline:
            current_index = [i for i, t in enumerate(self._times) if self.time() >= t][-1]
            current_event = self.timeline[current_index][1]
            if current_event == 'end':
                self.end_timeline = True
            elif current_event != self.media_id:
                self.media_id = self.timeline[current_index][1]
                new_media = self.make_media.get_saved_media(self.media_id)
                self.fill_lattice(new_media)

                print('Media condition: ' + str(self.media_id))

    def fill_lattice(self, media):
        # Create lattice and fill each site with concentrations dictionary
        # Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
        self.lattice = np.empty([len(self._molecule_ids)] + [self.patches_per_edge_x, self.patches_per_edge_y], dtype=np.float64)
        for index, molecule_id in enumerate(self._molecule_ids):
            self.lattice[index].fill(media[molecule_id])

    def run_diffusion(self, timestep):
        for index in range(len(self.lattice)):
            molecule = self.lattice[index]
            # run diffusion if molecule field is not uniform
            if len(set(molecule.flatten())) != 1:
                t=0.0
                while t<timestep:
                    molecule += self.diffusion_timestep(molecule, self.diffusion_dt)
                    t += self.diffusion_dt

    def diffusion_timestep(self, field, dt):
        ''' calculate concentration changes cause by diffusion'''
        change_field = self.diffusion * dt * convolve(field, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_field


    ## Conversion functions
    def volume_to_length(self, volume, radius):
        '''
        get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
        a=length of cylinder without rounded caps, l=total length:

        V = (4/3)*PI*r^3 + PI*r^2*a
        l = a + 2*r
        '''
        cylinder_length = (volume - (4/3) * PI * radius**3) / (PI * radius**2)
        total_length = cylinder_length + 2 * radius

        return total_length

    def count_to_concentration(self, count):
        ''' Convert count to concentrations '''
        return count / (self.patch_volume * N_AVOGADRO)


    # functions for getting values
    def simulation_parameters(self, agent_id):
        latest = max([
            simulation['time']
            for agent_id, simulation
            in self.simulations.items()])
        time = max(self._time, latest)

        return {
            'time': time}

    def simulation_state(self, agent_id):
        """
        Return the state the environment is tracking about the simulation given by `agent_id`.
        For use by a daughter cell.
        """
        return dict(
            self.simulations[agent_id],
            location=self.locations[agent_id])

    def get_molecule_ids(self):
        ''' Return the ids of all molecule species in the environment '''
        return self._molecule_ids

    def time(self):
        return self._time

    def run_for_time(self):
        return self.run_for

    def max_time(self):
        return self._max_time


    # Agent interface
    def run_incremental(self, run_until):
        ''' Simulate until run_until '''
        while self._time < run_until:
            self.evolve(self.run_for)
            self._time += self.run_for

    def add_simulation(self, agent_id, simulation):
        self.simulations.setdefault(agent_id, {}).update(simulation)

        if agent_id not in self.multicell_physics.cells:

            if simulation.get('parent_id', ''):
                index = simulation['index']
                parent_location = simulation['location'][0:2]
                orientation = simulation['location'][2]
                parent_volume = self.simulations[agent_id]['state']['volume'] * 2  # TODO -- get parent volume from other source
                parent_length = self.volume_to_length(parent_volume, self.cell_radius)

                daughter_locations = self.daughter_locations(parent_location, parent_length, orientation)
                location = daughter_locations[index]
            elif self.cell_placement:
                placement = np.array([
                    self.cell_placement[0] * self.edge_length_x,
                    self.cell_placement[1] * self.edge_length_y])
                location = placement
                orientation = np.random.uniform(0, 2 * PI)
            else:
                # Place cell at either the provided or a random initial location
                random_location = np.array([
                    np.random.uniform(0, self.edge_length_x),
                    np.random.uniform(0, self.edge_length_y)])
                location = simulation['agent_config'].get(
                    'location', random_location)
                orientation = simulation['agent_config'].get(
                    'orientation', np.random.uniform(0, 2 * PI))

            self.add_cell_to_physics(agent_id, location, orientation)

        # add to motile forces
        if agent_id not in self.motile_forces:
            self.motile_forces[agent_id] = [0.0, 0.0]

    def apply_inner_update(self, update, now):
        '''
        Use change counts from all the inner simulations, convert them to concentrations,
        and add to the environmental concentrations of each molecule at each simulation's location
        '''
        self.simulations.update(update)

        for agent_id, simulation in self.simulations.items():
            # only apply changes if we have reached this simulation's time point.
            if simulation['time'] <= now:
                # print('=== simulation update: {}'.format(simulation))
                state = simulation['state']

                if 'motile_force' in state:
                    self.motile_forces[agent_id] = state['motile_force']

                if not self.static_concentrations:
                    location = np.array([
                        self.locations[agent_id][0] * self.patches_per_edge_x / self.edge_length_x,
                        self.locations[agent_id][1] * self.patches_per_edge_y / self.edge_length_y
                    ])

                    patch_site = tuple(np.floor(location).astype(int))

                    for molecule, count in state['environment_change'].items():
                        concentration = self.count_to_concentration(count)
                        index = self.molecule_index[molecule]
                        self.lattice[index, patch_site[0], patch_site[1]] += concentration

    def generate_outer_update(self, now):
        '''Return a dict with {molecule_id: conc} for each sim at its current location'''
        update = {}
        for agent_id, simulation in self.simulations.items():
            # only provide concentrations if we have reached this simulation's time point.
            if simulation['time'] <= now:
                # get concentration from cell's given bin
                location = np.array([
                    self.locations[agent_id][0] * self.patches_per_edge_x / self.edge_length_x,
                    self.locations[agent_id][1] * self.patches_per_edge_y / self.edge_length_y
                ])
                patch_site = tuple(np.floor(location).astype(int))

                assert (0 <= patch_site[0] < self.patches_per_edge_x)
                assert (0 <= patch_site[1] < self.patches_per_edge_y)

                update[agent_id] = {}
                update[agent_id]['concentrations'] = dict(zip(
                    self._molecule_ids,
                    self.lattice[:, patch_site[0], patch_site[1]]))

                update[agent_id]['media_id'] = self.media_id

        return update

    def shutdown_environment(self, now):
        return self.end_timeline

    def daughter_locations(self, parent_location, parent_length, parent_angle):
        pos_ratios = [-0.25, 0.25]
        daughter_locations = []
        for daughter in range(2):
            dx = parent_length * pos_ratios[daughter] * math.cos(parent_angle)
            dy = parent_length * pos_ratios[daughter] * math.sin(parent_angle)
            location = parent_location + [dx, dy]

            daughter_locations.append(location)

        return daughter_locations

    def remove_simulation(self, agent_id):
        self.simulations.pop(agent_id, {})
        self.locations.pop(agent_id, {})
        self.multicell_physics.remove_cell(agent_id)


    # Emitters
    def emit_data(self):

        # emit lattice data
        emit_fields = {}
        for index, molecule_id in zip(self.emit_indices, self.emit_fields):
            emit_fields[molecule_id] = self.lattice[index].tolist()
        data = {
            'type': 'lattice-field',
            'fields': emit_fields,
            'time': self.time()}
        emit_config = {
            'table': 'history',
            'data': data}
        self.emitter.emit(emit_config)

        # emit data for each agent
        for agent_id, simulation in self.simulations.items():
            agent_location = self.locations[agent_id].tolist()  # [x, y, theta]
            agent_state = self.simulations[agent_id]['state']
            data = {
                'type': 'lattice',
                'agent_id': agent_id,
                'location': agent_location,
                'volume': agent_state['volume'],
                'width': agent_state['width'],
                'length': agent_state['length'],
                'time': self.time()}

            emit_config = {
                'table': 'history',
                'data': data}

            self.emitter.emit(emit_config)

    def emit_configuration(self):
        data = {
            'type': 'lattice',
            'cell_radius': self.cell_radius,
            'edge_length_x': self.edge_length_x,
            'edge_length_y': self.edge_length_y,
            'patches_per_edge_x': self.patches_per_edge_x,
            'patches_per_edge_y': self.patches_per_edge_y,
            'total_volume': self.total_volume,
            'timeline': self.timeline}

        emit_config = {
            'table': 'configuration',
            'data': data}

        self.emitter.emit(emit_config)


# tests
def tumble():
    tumble_jitter = 0.4  # (radians)
    force = 1.0
    torque = random.normalvariate(0, tumble_jitter)
    return [force, torque]

def run():
    force = 2.1
    torque = 0.0
    return [force, torque]


def test_diffusion(total_time=1):
    from vivarium.actor.emitter import get_emitter

    timestep = 0.01  # sec
    test_diffusion = 'GLC'
    edge_length = 100.0

    # get media
    media_id = 'GLC_G6P'
    make_media = Media()
    media = make_media.get_saved_media(media_id)

    # get emitter
    emitter = get_emitter({})  # TODO -- is an emitter really necessary?

    boot_config = {
        'concentrations': media,
        'run_for': timestep,
        'depth': 0.01,  # 3000 um is default
        'edge_length_x': edge_length,
        'patches_per_edge_x': int(edge_length/2),
        'diffusion': 1e3,
        'gradient': {
            'type': 'linear',
            'molecules': {
                test_diffusion: {
                    'center': [0.5, 0.5],
                    'slope': -20.0 / edge_length},
            }},
        'emitter': emitter
    }

    # configure lattice
    lattice = EnvironmentSpatialLattice(boot_config)

    # get test_diffusion field index
    test_diffusion_idx = lattice._molecule_ids.index(test_diffusion)
    center_patch = [int(boot_config['gradient']['molecules'][test_diffusion]['center'][0] * boot_config['patches_per_edge_x']),
              int(boot_config['gradient']['molecules'][test_diffusion]['center'][1] * boot_config['patches_per_edge_x'])]

    # run simulation
    time = 0
    previous_field = lattice.lattice[test_diffusion_idx].copy()
    saved_state = {
        'field': [previous_field],
        'time': [time]}
    while time < total_time:
        time += timestep

        # run lattice and get new locations
        lattice.run_incremental(time)

        # get field
        field = lattice.lattice[test_diffusion_idx].copy()

        # center should go down
        assert field[center_patch[0]][center_patch[1]] < previous_field[center_patch[0]][center_patch[1]]
        # corner should go up
        assert field[0][0] > previous_field[0][0]

        # save current field for comparison
        previous_field = field

        # save state for plotting
        saved_state['time'].append(time)
        saved_state['field'].append(field.tolist())

    data = {
        'saved_state': saved_state,
        'timestep': timestep,
        'x_length': edge_length,
        'y_length': edge_length,
    }
    return data


def test_lattice(total_time=10):
    from vivarium.actor.emitter import get_emitter

    timestep = 0.01  # sec
    edge_length = 100.0

    # time of motor behavior without chemotaxis
    run_time = 0.42  # s (Berg)
    tumble_time = 0.14  # s (Berg)

    # get emitter
    emitter = get_emitter({})  # TODO -- is an emitter really necessary?

    boot_config = {
        'translation_jitter': 0.0,
        'rotation_jitter': 0.0,
        'run_for': timestep,
        'depth': 0.01,  # 3000 um is default
        'static_concentrations': True,
        'edge_length_x': edge_length,
        'patches_per_edge_x': int(edge_length/2),
        'cell_placement': [0.5, 0.5],  # place cells at center of lattice
        'emitter': emitter
    }

    # configure lattice
    lattice = EnvironmentSpatialLattice(boot_config)

    # get simulations
    simulations = lattice.simulations

    # add an agent
    agent_id = '1'
    simulation = simulations.setdefault(agent_id, {})
    agent_state = {'volume': 1.0}
    agent_config = {
        'location': np.array([0.0, 0.0]),
        'orientation': np.array([0.0]),
    }
    simulation.update({
        'time': lattice.time(),
        'state': agent_state,
        'agent_config': agent_config,
        })
    lattice.add_simulation(agent_id, simulation)

    # run simulation
    saved_state = {
        'location': [],
        'time': [],
        'motor_state': []}

    # test run/tumble
    time = 0
    motor_state = 1 # 0 for run, 1 for tumble
    time_in_motor_state = 0
    while time < total_time:

        # get forces
        if motor_state == 1: # tumble
            if time_in_motor_state < tumble_time:
                motile_force = tumble()
                time_in_motor_state += timestep
            else:
                # switch
                motile_force = run()
                motor_state = 0
                time_in_motor_state = 0

        elif motor_state == 0:  # run
            if time_in_motor_state < run_time:
                motile_force = run()
                time_in_motor_state += timestep
            else:
                # switch
                motile_force = tumble()
                motor_state = 1
                time_in_motor_state = 0

        # apply forces
        lattice.motile_forces[agent_id] = motile_force

        # run lattice and get new locations
        lattice.run_incremental(time)
        locations = lattice.locations  # new location

        # save
        saved_state['time'].append(time)
        saved_state['location'].append(list(locations[agent_id]))
        saved_state['motor_state'].append(motor_state)

        # update time
        time += timestep

    # TODO -- assert division
    # TODO -- assert low reynolds number check

    data = {
        'saved_state': saved_state,
        'timestep': timestep,
        'x_length': edge_length,
        'y_length': edge_length,
    }
    return data


def plot_motility(data, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    expected_speed = 14.2  # um/s (Berg)
    expected_angle_between_runs = 68 # degrees (Berg)

    saved_state = data['saved_state']
    # timestep = data['timestep']
    locations = saved_state['location']
    motor_state = saved_state['motor_state'] # 0 for run, 1 for tumble
    times = saved_state['time']

    # get speed, angle between runs
    speed_vec = [0]
    previous_time = times[0]
    previous_loc = locations[0]
    previous_motor_state = motor_state[0]
    run_angle = 0
    angle_between_runs = []
    for time, location, m_state in zip(times[1:], locations[1:], motor_state[1:]):

        # get speed
        dt = time - previous_time
        distance = ((location[0] - previous_loc[0]) ** 2 + (location[1] - previous_loc[1]) ** 2) ** 0.5
        speed_vec.append(distance / dt)  # um/sec

        # get change in angles between motor states
        if m_state == 0 and previous_motor_state == 1:
            angle_change = abs(location[2] - run_angle) / np.pi * 180 % 360 # convert to absolute degrees
            angle_between_runs.append(angle_change)
        elif m_state == 0:
            run_angle = location[2]

        # update previous states
        previous_time = time
        previous_loc = location
        previous_motor_state = m_state

    avg_speed = sum(speed_vec) / len(speed_vec)
    avg_angle_between_runs = sum(angle_between_runs) / len(angle_between_runs)

    # plot results
    cols = 1
    rows = 3
    plt.figure(figsize=(6 * cols, 1.5 * rows))
    plt.rcParams.update({'font.size': 12})

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)

    ax1.plot(times, speed_vec)
    ax1.axhline(y=avg_speed, color='b', linestyle='dashed', label='mean')
    ax1.axhline(y=expected_speed, color='r', linestyle='dashed', label='expected mean')
    ax1.set_ylabel(u'speed (\u03bcm/sec)')
    # ax1.set_xlabel('time')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plot change in direction between runs
    ax2.plot(angle_between_runs)
    ax2.axhline(y=avg_angle_between_runs, color='b', linestyle='dashed', label='mean angle between runs')
    ax2.axhline(y=expected_angle_between_runs, color='r', linestyle='dashed', label='exp. angle between runs')
    ax2.set_ylabel('angle between runs')
    ax2.set_xlabel('run #')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # motile force
    ax3.plot(times, motor_state)
    ax3.set_xlabel('time (s)')
    ax3.set_yticks([0.0, 1.0])
    ax3.set_yticklabels(["run", "tumble"])

    # save figure
    fig_path = os.path.join(out_dir, 'motility')
    plt.subplots_adjust(wspace=0.7, hspace=0.9)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


def plot_trajectory(data, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    x_length = data['x_length']
    y_length = data['y_length']
    y_ratio = y_length/x_length
    saved_state = data['saved_state']
    locations = saved_state['location']
    times = np.array(saved_state['time'])

    # plot trajectory
    fig = plt.figure(figsize=(8, 8*y_ratio))

    # get locations and convert to 2D array
    locations_array = np.array(locations)
    x_coord = locations_array[:, 0]
    y_coord = locations_array[:, 1]

    # make multi-colored trajectory
    points = np.array([x_coord, y_coord]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=plt.get_cmap('cool'))
    lc.set_array(times)
    lc.set_linewidth(3)

    # plot line
    line = plt.gca().add_collection(lc)

    # color bar
    cbar = plt.colorbar(line)
    cbar.set_label('time', rotation=270)

    # plt.plot(x_coord, y_coord, 'b-')  # trajectory
    plt.plot(x_coord[0], y_coord[0], color=(0.0, 0.8, 0.0), marker='*')  # starting point
    plt.plot(x_coord[-1], y_coord[-1], color='r', marker='*')  # ending point

    plt.xlim((0, x_length))
    plt.ylim((0, y_length))

    fig_path = os.path.join(out_dir, 'trajectory')
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')
    plt.close(fig)


def plot_field(data, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    x_length = data['x_length']
    y_length = data['y_length']
    saved_state = data['saved_state']
    fields = saved_state['field']
    times = saved_state['time']

     # plot fields
    time_vec = times
    n_snapshots = 6
    n_fields = 1
    plot_steps = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int).tolist()
    snapshot_times = [time_vec[i] for i in plot_steps]

    # make figure
    fig = plt.figure(figsize=(20 * n_snapshots, 10))
    grid = plt.GridSpec(n_fields, n_snapshots, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})
    for index, time_index in enumerate(plot_steps, 0):

        field_data = fields[time_index]

        ax = fig.add_subplot(grid[0, index])  # grid is (row, column)
        # ax.title.set_text('time: {:.4f} hr | field: {}'.format(float(time) / 60. / 60., field_id))
        ax.set(xlim=[0, x_length], ylim=[0, y_length], aspect=1)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

        plt.imshow(field_data,
                   vmin=0,
                   vmax=15.0,
                   origin='lower',
                   extent=[0, x_length, 0, y_length],
                   interpolation='nearest',
                   cmap='YlGn')

    plt.colorbar()
    fig_path = os.path.join(out_dir, 'field')
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output1 = test_lattice(50)
    plot_motility(output1, out_dir)
    plot_trajectory(output1, out_dir)

    # output2 = test_diffusion(2)
    # plot_field(output2, out_dir)
