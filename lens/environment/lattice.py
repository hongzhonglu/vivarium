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
import numpy as np
from scipy import constants
from scipy.ndimage import convolve

from lens.actor.outer import EnvironmentSimulation
from lens.utils.multicell_physics import MultiCellPhysics
from lens.environment.make_media import Media

# Constants
N_AVOGADRO = constants.N_A
PI = math.pi

CELL_DENSITY = 1100

# Lattice parameters
TRANSLATION_JITTER = 0.1
ROTATION_JITTER = 0.05

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])


def gaussian(deviation, distance):
    return np.exp(-np.power(distance, 2.) / (2 * np.power(deviation, 2.)))

class EnvironmentSpatialLattice(EnvironmentSimulation):
    def __init__(self, config):

        ## Simulation Parameters
        self._time = 0
        self._timestep = 1.0
        self._max_time = 10e6
        self.output_dir = config.get('output_dir')
        self.run_for = config.get('run_for', 5.0)

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
        for new_media_id, media_dict in new_media.iteritems():
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

        self._molecule_ids = media.keys()
        self.concentrations = media.values()
        self.molecule_index = {molecule: index for index, molecule in enumerate(self._molecule_ids)}
        self.fill_lattice(media)

        ## Concentration and Gradient Parameters
        self.static_concentrations = config.get('static_concentrations', False)
        self.diffusion = config.get('diffusion', 0.1)
        self.gradient = config.get('gradient', {'type': False})

        # upper limit on the time scale. dy is assumed to be the same as dx. (using 50% of the theoretical upper limit)
        self.dt = 0.5 * self.dx2 * self.dx2 / (2 * self.diffusion * (self.dx2 + self.dx2)) if self.diffusion else 0

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

            for molecule_id, specs in self.gradient['molecules'].iteritems():
                mol_index = self._molecule_ids.index(molecule_id)
                center = [specs['center'][0] * self.edge_length_x,
                          specs['center'][1] * self.edge_length_y]
                deviation = specs['deviation']

                for x_patch in xrange(self.patches_per_edge_x):
                    for y_patch in xrange(self.patches_per_edge_y):
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

            for molecule_id, specs in self.gradient['molecules'].iteritems():
                mol_index = self._molecule_ids.index(molecule_id)
                center = [specs['center'][0] * self.edge_length_x,
                          specs['center'][1] * self.edge_length_y]
                slope = specs['slope']

                for x_patch in xrange(self.patches_per_edge_x):
                    for y_patch in xrange(self.patches_per_edge_y):
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


    def evolve(self):
        ''' Evolve environment '''
        self.update_locations()
        self.update_media()

        if not self.static_concentrations:
            self.run_diffusion()

        # make sure all patches have concentrations of 0 or higher
        self.lattice[self.lattice < 0.0] = 0.0

        # emit current state
        self.emit_data()

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
            elif self.cell_placement:
                placement = np.array([
                    self.cell_placement[0] * self.edge_length_x,
                    self.cell_placement[1] * self.edge_length_y])
                location = simulation['agent_config'].get(
                    'location', placement)
                orientation = simulation['agent_config'].get(
                    'orientation', np.random.uniform(0, 2 * PI))
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

        for agent_id, simulation in self.simulations.iteritems():
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
        for agent_id, simulation in self.simulations.iteritems():
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



def test_lattice():
    from lens.actor.emitter import get_emitter

    # get media
    media_id = 'GLC_G6P'
    make_media = Media()
    media = make_media.get_saved_media(media_id)

    # get emitter
    emitter = get_emitter({})

    boot_config = {
        'concentrations': media,
        'run_for': 10.0,
        'depth': 0.0001,  # 3000 um is default
        'edge_length': 10.0,
        'patches_per_edge': 1,
        'emitter': emitter
    }

    # configure lattice
    lattice = EnvironmentSpatialLattice(boot_config)

    # get simulations
    simulations = lattice.simulations

    # add cell simulation
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

    # apply forces
    run_for = 1.0
    force = [1.0, 1.0]
    lattice.motile_forces[agent_id] = force

    # run lattice and get new locations
    lattice.run_incremental(run_for)
    locations = lattice.locations  # new location

    import ipdb; ipdb.set_trace()

    # diffusion

    # division
    return {}

def plot_lattice(output, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    # plot
    plt.figure()

    fig_path = os.path.join(out_dir, 'lattice')
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    output = test_lattice()
    plot_lattice(output, out_dir)
