from __future__ import absolute_import, division, print_function

import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import random
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb
import pymunk

# vivarium imports
from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process,
    convert_to_timeseries)


# constants
PI = math.pi
ELASTICITY = 0.95
DAMPING = 0.05  # simulates viscous forces to reduce velocity at low Reynolds number (1 = no damping, 0 = full damping)
ANGULAR_DAMPING = 0.7  # less damping for angular velocity seems to improve behavior
FRICTION = 0.9  # TODO -- does this do anything?
PHYSICS_TS = 0.005
FORCE_SCALING = 15  # scales from pN
JITTER_FORCE = 1e-3  # pN

DEFAULT_BOUNDS = [10, 10]

# colors for phylogeny initial bodies
HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]

def random_body_position(body):
    ''' pick a random point along the boundary'''
    width, length = body.dimensions
    if random.randint(0, 1) == 0:
        # force along ends
        if random.randint(0, 1) == 0:
            # force on the left end
            location = (random.uniform(0, width), 0)
        else:
            # force on the right end
            location = (random.uniform(0, width), length)
    else:
        # force along length
        if random.randint(0, 1) == 0:
            # force on the bottom end
            location = (0, random.uniform(0, length))
        else:
            # force on the top end
            location = (width, random.uniform(0, length))
    return location



class Multibody(Process):
    ''''''
    def __init__(self, initial_parameters={}):

        # hardcoded parameters
        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.damping = DAMPING
        self.angular_damping = ANGULAR_DAMPING
        self.force_scaling = FORCE_SCALING
        self.physics_dt = PHYSICS_TS

        # configured parameters
        self.jitter_force = initial_parameters.get('jitter_force', JITTER_FORCE)
        bounds = initial_parameters.get('bounds', DEFAULT_BOUNDS)

        # initialize pymunk space
        self.space = pymunk.Space()

        # add static barriers
        # TODO -- mother machine configuration
        self.add_barriers(bounds)

        # initialize bodies
        self.bodies = {}
        bodies = initial_parameters.get('bodies')
        for body_id, specs in bodies.items():
            self.add_body_from_center(body_id, specs)

        # make ports
        # TODO -- need option to add/remove ports throughout simulation
        ports = {'bodies': list(self.bodies.keys())}

        parameters = {}
        parameters.update(initial_parameters)

        super(Multibody, self).__init__(ports, parameters)

    def default_settings(self):
        bodies = {body_id: self.get_body_specs(body_id)
            for body_id in self.bodies.keys()}

        schema = {'bodies': {body_id : {'updater': 'merge'}
            for body_id, body in bodies.items()}}

        return {
            'state': {'bodies': bodies},
            'schema': schema,}

    def next_update(self, timestep, states):
        bodies = states['bodies']

        # TODO -- check if any body has been removed?
        for body_id, specs in bodies.items():
            if body_id in self.bodies:
                self.update_body(body_id, specs)
            else:
                self.add_body_from_center(body_id, specs)

        # run simulation
        self.run(timestep)

        # get new bodies specs
        new_bodies = {
            body_id: self.get_body_specs(body_id)
                for body_id in self.bodies.keys()}

        return {'bodies': new_bodies}

    def run(self, timestep):
        assert self.physics_dt < timestep

        time = 0
        while time < timestep:
            time += self.physics_dt

            # apply forces
            for body in self.space.bodies:
                self.apply_jitter_force(body)
                self.apply_motile_force(body)
                self.apply_viscous_force(body)

            # run for a physics timestep
            self.space.step(self.physics_dt)

    def apply_motile_force(self, body):
        width, length = body.dimensions

        # motile forces
        motile_location = (width / 2, 0)  # apply force at back end of body
        motile_force = [0.0, 0.0]
        if hasattr(body, 'motile_force'):
            thrust, torque = body.motile_force
            motile_force = [thrust, 0.0]

            # add directly to angular velocity
            body.angular_velocity += torque
            ## force-based torque
            # if torque != 0.0:
            #     motile_force = get_force_with_angle(thrust, torque)

        scaled_motile_force = [thrust * self.force_scaling for thrust in motile_force]
        body.apply_force_at_local_point(scaled_motile_force, motile_location)

    def apply_jitter_force(self, body):
        jitter_location = random_body_position(body)
        jitter_force = [
            random.normalvariate(0, self.jitter_force),
            random.normalvariate(0, self.jitter_force)]
        scaled_jitter_force = [force * self.force_scaling for force in jitter_force]
        body.apply_force_at_local_point(scaled_jitter_force, jitter_location)

    def apply_viscous_force(self, body):
        # dampen the velocity
        body.velocity = body.velocity * self.damping + (body.force / body.mass) * self.physics_dt
        body.angular_velocity = body.angular_velocity * self.angular_damping + body.torque / body.moment * self.physics_dt

    def add_barriers(self, bounds):
        """ Create static barriers """
        x_bound = bounds[0]
        y_bound = bounds[1]

        static_body = self.space.static_body
        static_lines = [
            pymunk.Segment(static_body, (0.0, 0.0), (x_bound, 0.0), 0.0),
            pymunk.Segment(static_body, (x_bound, 0.0), (x_bound, y_bound), 0.0),
            pymunk.Segment(static_body, (x_bound, y_bound), (0.0, y_bound), 0.0),
            pymunk.Segment(static_body, (0.0, y_bound), (0.0, 0.0), 0.0),
        ]
        for line in static_lines:
            line.elasticity = 0.0  # no bounce
            line.friction = 0.9
        self.space.add(static_lines)

    def add_body_from_center(self, body_id, body):
        width = body['width']
        length = body['length']
        mass = body['mass']
        center_position = body['location']
        angle = body['angle']
        angular_velocity = body.get('angular_velocity', 0.0)

        half_length = length / 2
        half_width = width / 2

        shape = pymunk.Poly(None, (
            (-half_length, -half_width),
            (half_length, -half_width),
            (half_length, half_width),
            (-half_length, half_width)))

        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = (center_position[0], center_position[1])
        body.angle = angle
        body.dimensions = (width, length)
        body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add body to bodies dictionary
        self.bodies[body_id] = (body, shape)

    def update_body(self, body_id, specs):

        length = specs.get('length')
        width = specs.get('width')
        mass = specs.get('mass')

        body, shape = self.bodies[body_id]
        position = body.position
        angle = body.angle

        # make shape, moment of inertia, and add a body
        half_length = length/2
        half_width = width/2
        new_shape = pymunk.Poly(None, (
            (-half_length, -half_width),
            (half_length, -half_width),
            (half_length, half_width),
            (-half_length, half_width)))

        inertia = pymunk.moment_for_poly(mass, new_shape.get_vertices())
        new_body = pymunk.Body(mass, inertia)
        new_shape.body = new_body

        new_body.position = position
        new_body.angle = angle
        new_body.angular_velocity = body.angular_velocity
        new_body.dimensions = (width, length)

        new_shape.elasticity = shape.elasticity
        new_shape.friction = shape.friction

        # swap bodies
        self.space.add(new_body, new_shape)
        self.space.remove(body, shape)

        # update body
        self.bodies[body_id] = (new_body, new_shape)

    def get_body_specs(self, body_id):
        body, shape = self.bodies[body_id]
        width, length = body.dimensions
        position = body.position

        # TODO -- velocity? angular velocity? forces?
        return {
            'location': [position[0], position[1]],
            'angle': body.angle,
            'length': length,
            'width': width,
            'mass': body.mass}



# test functions
def n_body_config(n_bodies=10, bounds=[10, 10]):
    bodies = {
        body_id: {
            'location': [
                np.random.uniform(0, bounds[0]),
                np.random.uniform(0, bounds[1])],
            'angle': np.random.uniform(0, 2 * PI),
            'volume': 1,
            'length': 1.0,
            'width': 0.5,
            'mass': 1,
            'forces': [0, 0]}
        for body_id in range(n_bodies)}

    return {
        'bodies': bodies,
        'bounds': bounds,
        'jitter_force': 1e1}


def plot_body(ax, data, color):

    # location, orientation, length
    x_center = data['location'][0]
    y_center = data['location'][1]
    theta = data['angle'] / PI * 180 + 90 # rotate 90 degrees to match field
    length = data['length']
    width = data['width']

    # get bottom left position
    x_offset = (width / 2)
    y_offset = (length / 2)
    theta_rad = math.radians(theta)
    dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
    dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

    x = x_center - dx
    y = y_center - dy

    # get color, convert to rgb
    rgb = hsv_to_rgb(color)

    # Create a rectangle
    rect = patches.Rectangle(
        (x, y), width, length, angle=theta, linewidth=2, edgecolor='w', facecolor=rgb)

    ax.add_patch(rect)


def plot_snapshots(data, config, out_dir='out', filename='multibody'):
    n_snapshots = 6
    bounds = config.get('bounds', DEFAULT_BOUNDS)

    # time steps that will be used
    time_vec = data['time']
    time_indices = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
    snapshot_times = [time_vec[i] for i in time_indices]

    # get bodies
    bodies = data['bodies']

    body_colors = {}
    for body_id in bodies:
        hue = random.choice(HUES)  # select random initial hue
        color = [hue] + DEFAULT_SV
        body_colors[body_id] = color

    # make the figure
    n_rows = 1
    n_cols = n_snapshots + 1  # one column for the colorbar
    fig = plt.figure(figsize=(12 * n_cols, 12 * n_rows))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    # plot snapshot data in each subsequent column
    for col_idx, (time_idx, time) in enumerate(zip(time_indices, snapshot_times), 1):
        row_idx = 0
        ax = init_axes(fig, bounds[0], bounds[1], grid, row_idx, col_idx, time)
        for body_id, series in bodies.items():
            body_data = {
                'location': series[time_idx]['location'],
                'angle': series[time_idx]['angle'],
                'length': series[time_idx]['length'],
                'width': series[time_idx]['width']}
            color = body_colors[body_id]
            plot_body(ax, body_data, color)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def init_axes(fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.4f} hr'.format(float(time) / 60. / 60.)
        plt.title(plot_title, y=1.08)
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax

def test_multibody(config=n_body_config(), time=1):
    multibody = Multibody(config)
    settings = {
        'total_time': time,
        'environment_port': 'external',
        'environment_volume': 1e-2}
    return simulate_process(multibody, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'multibody')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = n_body_config(10)
    saved_data = test_multibody(config, 20)
    timeseries = convert_to_timeseries(saved_data)
    plot_snapshots(timeseries, config, out_dir, 'bodies')
