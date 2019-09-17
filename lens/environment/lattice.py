"""
Lattice

A two-dimensional lattice environmental model

## Physics
# Diffusion constant of glucose in 0.5 and 1.5 percent agarose gel = ~6 * 10^-10 m^2/s (Weng et al. 2005. Transport of glucose and poly(ethylene glycol)s in agarose gels).
# Conversion to micrometers: 6 * 10^-10 m^2/s = 600 micrometers^2/s.

## Cell biophysics
# rotational diffusion in liquid medium with viscosity = 1 mPa.s: Dr = 3.5+/-0.3 rad^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
# translational diffusion in liquid medium with viscosity = 1 mPa.s: Dt=100 micrometers^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)


@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import os
import math
import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from lens.environment.make_media import Media
from lens.actor.outer import EnvironmentSimulation
from lens.utils.multicell_physics import MultiCellPhysics

# Constants
N_AVOGADRO = constants.N_A
PI = math.pi

CELL_DENSITY = 1100

# Lattice parameters
N_DIMS = 2

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


def gaussian(deviation, distance):
    return np.exp(-np.power(distance, 2.) / (2 * np.power(deviation, 2.)))

class EnvironmentSpatialLattice(EnvironmentSimulation):
    def __init__(self, config):
        self._time = 0
        self._timestep = 1.0
        self._max_time = 10e6

        # constants
        self.cell_density = CELL_DENSITY  # TODO (Eran) -- get mass from cell sim rather than using density to convert volume to mass

        # configured parameters
        self.output_dir = config.get('output_dir')
        self.run_for = config.get('run_for', 5.0)
        self.edge_length = config.get('edge_length', 10.0)
        self.patches_per_edge = config.get('patches_per_edge', 10)
        self.cell_radius = config.get('cell_radius', 0.5)
        self.static_concentrations = config.get('static_concentrations', False)
        self.diffusion = config.get('diffusion', 0.1)
        self.gradient = {
            'seed': False,
            'molecules': {
                'GLC':{
                    'center': [0.5, 0.5],
                    'deviation': 10.0},
            }}
        self.gradient.update(config.get('gradient', {}))
        self.translation_jitter = config.get('translation_jitter', 0.5)
        self.rotation_jitter = config.get('rotation_jitter', 0.005)
        self.depth = config.get('depth', 3000.0)  # um
        self.timeline = config.get('timeline')
        self.media_id = config.get('media_id', 'minimal')
        if self.timeline:
            self._times = [t[0] for t in self.timeline]

        # derived parameters
        self.total_volume = (self.depth * self.edge_length ** 2) * (10 ** -15) # (L)
        self.patch_volume = self.total_volume / (self.patches_per_edge ** 2)
        # intervals in x- directions (assume y- direction equivalent)
        self.dx = self.edge_length / self.patches_per_edge
        self.dx2 = self.dx * self.dx
        # upper limit on the time scale (go with at least 50% of this)
        self.dt = 0.5 * self.dx2 * self.dx2 / (2 * self.diffusion * (self.dx2 + self.dx2)) if self.diffusion else 0

        self.simulations = {}       # map of agent_id to simulation state
        self.locations = {}         # map of agent_id to center location and orientation
        self.corner_locations = {}  # map of agent_id to corner location, for Lens visualization and multi-cell physics engine
        self.motile_forces = {}	    # map of agent_id to motile force, with magnitude and relative orientation

        # make physics object by passing in bounds and jitter
        bounds = [self.edge_length, self.edge_length]
        self.multicell_physics = MultiCellPhysics(
            bounds,
            self.translation_jitter,
            self.rotation_jitter)

        # make media object for making new media
        self.make_media = Media()  # TODO (Eran) -- pass in make_media through config to include media from all timelines!

        # make media and fill lattice patches with media concentrations
        media = config['concentrations']
        self._molecule_ids = config['concentrations'].keys()
        self.concentrations = config['concentrations'].values()
        self.molecule_index = {molecule: index for index, molecule in enumerate(self._molecule_ids)}
        self.fill_lattice(media)

        # Add gradient
        if self.gradient['seed']:
            for molecule_id, specs in self.gradient['molecules'].iteritems():
                center = [x * self.edge_length for x in specs['center']]
                deviation = specs['deviation']
                for x_patch in xrange(self.patches_per_edge):
                    for y_patch in xrange(self.patches_per_edge):
                        # distance from middle of patch to center coordinates
                        dx = (x_patch + 0.5) * self.edge_length / self.patches_per_edge - center[0]
                        dy = (y_patch + 0.5) * self.edge_length / self.patches_per_edge - center[1]
                        distance = np.sqrt(dx ** 2 + dy ** 2)
                        scale = gaussian(deviation, distance)
                        # multiply glucose gradient by scale
                        self.lattice[self._molecule_ids.index(molecule_id)][x_patch][y_patch] *= scale

        # # Output
        # self.create_lattice_table()
        # # Track agent tables
        # self.agent_tables = {}

    def evolve(self):
        ''' Evolve environment '''
        self.update_locations()
        self.update_media()

        if not self.static_concentrations:
            self.run_diffusion()

        # make sure all patches have concentrations of 0 or higher
        self.lattice[self.lattice < 0.0] = 0.0

        # self.append_agent_tables()  # TODO (Eran) -- use emitter

    def update_locations(self):
        ''' Update location for all agent_ids '''
        for agent_id, location in self.locations.iteritems():

            # shape
            agent_state = self.simulations[agent_id]['state']
            volume = agent_state['volume']
            radius = self.cell_radius
            width = 2 * radius
            length = self.volume_to_length(volume, radius)
            mass = agent_state.get('mass', 1.0)  # TODO -- pass mass through state message

            # update length, width, update in multicell_physics
            agent_state['length'] = length
            agent_state['width'] = width
            self.multicell_physics.update_cell(agent_id, length, width, mass)

            # Motile forces
            force = self.motile_forces[agent_id][0]
            torque = self.motile_forces[agent_id][1]
            self.multicell_physics.apply_motile_force(agent_id, force, torque)

        self.multicell_physics.run_incremental(self.run_for)

        for agent_id, location in self.locations.iteritems():
            # update location
            self.locations[agent_id] = self.multicell_physics.get_center(agent_id)
            self.corner_locations[agent_id] = self.multicell_physics.get_corner(agent_id)

            # enforce boundaries # TODO (Eran) -- make pymunk handle this
            self.locations[agent_id][0:2][self.locations[agent_id][0:2] > self.edge_length] = self.edge_length - self.dx / 2
            self.locations[agent_id][0:2][self.locations[agent_id][0:2] < 0] = 0.0

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
            if self.media_id != self.timeline[current_index][1]:
                self.media_id = self.timeline[current_index][1]
                new_media = self.make_media.get_saved_media(self.media_id)
                self.fill_lattice(new_media)

                print('Media condition: ' + str(self.media_id))

    def fill_lattice(self, media):
        # Create lattice and fill each site with concentrations dictionary
        # Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
        self.lattice = np.empty([len(self._molecule_ids)] + [self.patches_per_edge for dim in xrange(N_DIMS)], dtype=np.float64)
        for index, molecule_id in enumerate(self._molecule_ids):
            self.lattice[index].fill(media[molecule_id])

    def run_diffusion(self):
        change_lattice = np.zeros(self.lattice.shape)
        for index in xrange(len(self.lattice)):
            molecule = self.lattice[index]

            # run diffusion if molecule field is not uniform
            if len(set(molecule.flatten())) != 1:
                change_lattice[index] = self.diffusion_timestep(molecule)

        self.lattice += change_lattice

    def diffusion_timestep(self, lattice):
        ''' calculate concentration changes cause by diffusion'''
        change_lattice = self.diffusion * self._timestep * convolve(lattice, LAPLACIAN_2D, mode='reflect') / self.dx2
        return change_lattice


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
            in self.simulations.iteritems()])
        time = max(self._time, latest)

        return {
            'time': time,
        }

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
            self._time += self._timestep
            self.evolve()

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
            else:
                # Place cell at either the provided or a random initial location
                location = simulation['agent_config'].get(
                    'location', np.random.uniform(0, self.edge_length, N_DIMS))
                orientation = simulation['agent_config'].get(
                    'orientation', np.random.uniform(0, 2 * PI))

            self.add_cell_to_physics(agent_id, location, orientation)

        # add to motile forces
        if agent_id not in self.motile_forces:
            self.motile_forces[agent_id] = [0.0, 0.0]

        # create output file for each cell to log location data
        # self.create_agent_table(agent_id)  # TODO -- get tableWriter path

    def apply_inner_update(self, update, now):
        '''
        Use change counts from all the inner simulations, convert them to concentrations,
        and add to the environmental concentrations of each molecule at each simulation's location
        '''
        self.simulations.update(update)

        for agent_id, simulation in self.simulations.iteritems():
            # only apply changes if we have reached this simulation's time point.
            if simulation['time'] <= now:
                # print('=== simulation update: {}'.format(simulation))
                state = simulation['state']

                if 'motile_force' in state:
                    self.motile_forces[agent_id] = state['motile_force']

                if not self.static_concentrations:
                    location = self.locations[agent_id][0:2] * self.patches_per_edge / self.edge_length
                    patch_site = tuple(np.floor(location).astype(int))

                    for molecule, count in state['environment_change'].iteritems():
                        concentration = self.count_to_concentration(count)
                        index = self.molecule_index[molecule]
                        self.lattice[index, patch_site[0], patch_site[1]] += concentration

    def generate_outer_update(self, now):
        '''Return a dict with {molecule_id: conc} for each sim at its current location'''
        update = {}
        for agent_id, simulation in self.simulations.iteritems():
            # only provide concentrations if we have reached this simulation's time point.
            if simulation['time'] <= now:
                # get concentration from cell's given bin
                location = self.locations[agent_id][0:2] * self.patches_per_edge / self.edge_length
                patch_site = tuple(np.floor(location).astype(int))

                assert (0 <= patch_site[0] < self.patches_per_edge)
                assert (0 <= patch_site[1] < self.patches_per_edge)

                update[agent_id] = {}
                update[agent_id]['concentrations'] = dict(zip(
                    self._molecule_ids,
                    self.lattice[:, patch_site[0], patch_site[1]]))

                update[agent_id]['media_id'] = self.media_id

        return update

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
        # self.agent_tables.pop(agent_id, {})
        self.multicell_physics.remove_cell(agent_id)

    # TODO (Eran) -- do all logging through emitters
    # # TableWriter functions
    # def create_lattice_table(self):
    #     table = TableWriter(self.output_dir)
    #     table.writeAttributes(
    #         edge_length=self.edge_length,
    #     )
    #
    # def create_agent_table(self, agent_id):
    #     # TODO (Eran) -- why is this called more than once for every agent? It should only be called once at initialization
    #     if agent_id not in self.agent_tables:
    #         table_writer_file = os.path.join(self.output_dir, agent_id)
    #         table = TableWriter(table_writer_file)
    #         table.writeAttributes(
    #             start_time = self.time(),
    #         )
    #         self.agent_tables[agent_id] = table
    #
    # def append_agent_tables(self):
    #     for agent_id, table in self.agent_tables.iteritems():
    #         agent_location = self.locations[agent_id]
    #
    #         agent_state = self.simulations[agent_id]['state']
    #         agent_volume = agent_state['volume']
    #         agent_width = agent_state['width']
    #         agent_length = agent_state['length']
    #
    #         table.append(
    #             location=agent_location,
    #             volume=agent_volume,
    #             width=agent_width,
    #             length=agent_length
    #         )
