from __future__ import absolute_import, division, print_function

import os

os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import pygame
from pygame.locals import *
from pygame.color import *

# Python imports
import random
import math
import numpy as np

# pymunk imports
import pymunkoptions
pymunkoptions.options["debug"] = False
import pymunk
import pymunk.pygame_util

PI = math.pi

ELASTICITY = 0.95
DAMPING = 0.05  # simulates viscous forces to reduce velocity at low Reynolds number (1 = no damping, 0 = full damping)
ANGULAR_DAMPING = 0.7  # less damping for angular velocity seems to improve behavior
FRICTION = 0.9  # TODO -- does this do anything?
PHYSICS_TS = 0.005
FORCE_SCALING = 15  # scales from pN


def get_force_with_angle(force, angle):
    x = force * math.cos(angle)
    y = force * math.sin(angle)
    return [x, y]

def front_from_corner(width, length, corner_position, angle):
    half_width = width/2
    dx = length * math.cos(angle) + half_width * math.cos(angle + PI/2)  # PI/2 gives a half-rotation for the width component
    dy = length * math.sin(angle) + half_width * math.sin(angle + PI/2)
    front_position = [corner_position[0] + dx, corner_position[1] + dy]
    return np.array([front_position[0], front_position[1], angle])

def corner_from_center(width, length, center_position, angle):
    half_length = length/2
    half_width = width/2
    dx = half_length * math.cos(angle) + half_width * math.cos(angle + PI/2)
    dy = half_length * math.sin(angle) + half_width * math.sin(angle + PI/2)
    corner_position = [center_position[0] - dx, center_position[1] - dy]

    return np.array([corner_position[0], corner_position[1], angle])

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



class MultiCellPhysics(object):
    ''''''
    def __init__(self, bounds, jitter_force, debug=False, debug_scale=20):

        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.damping = DAMPING
        self.angular_damping = ANGULAR_DAMPING
        self.force_scaling = FORCE_SCALING
        self.jitter_force = jitter_force

        # Space
        self.space = pymunk.Space()

        # Physics
        self.physics_dt = PHYSICS_TS

        # Debugging with pygame
        self.pygame_viz = debug
        self.pygame_scale = 1
        if self.pygame_viz:
            self.pygame_scale = debug_scale
            pygame.init()
            self._screen = pygame.display.set_mode((
                int(bounds[0]*self.pygame_scale), int(bounds[1]*self.pygame_scale)))
            self._clock = pygame.time.Clock()
            self._draw_options = pymunk.pygame_util.DrawOptions(self._screen)

        # Static barriers
        self.add_barriers(bounds)

        # Cells
        self.cells = {}

    # pygame functions
    def _process_events(self):
        for event in pygame.event.get():
            if event.type == QUIT:
                self._running = False
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                self._running = False

    def _clear_screen(self):
        self._screen.fill(THECOLORS["white"])

    def _draw_objects(self):
        self.space.debug_draw(self._draw_options)

    def _update_screen(self):
        self._process_events()
        self._clear_screen()
        self._draw_objects()
        pygame.display.flip()
        # Delay fixed time between frames
        self._clock.tick(5)

    def update_motile_force(self, cell_id, force, torque):
        body, shape = self.cells[cell_id]
        body.motile_force = (force, torque)

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
        # jitter forces
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

    def run_incremental(self, run_for):
        assert self.physics_dt < run_for

        time = 0
        while time < run_for:
            time += self.physics_dt

            # apply forces
            for body in self.space.bodies:
                self.apply_jitter_force(body)
                self.apply_motile_force(body)
                self.apply_viscous_force(body)

            self.space.step(self.physics_dt)

        if self.pygame_viz:
            self._update_screen()

    def update_cell(self, cell_id, length, width, mass):
        ''' create a new body and new shape at the cell's same center position and angle '''

        body, shape = self.cells[cell_id]
        position = body.position
        angle = body.angle

        # make shape, moment of inertia, and add a body
        half_length = length/2 * self.pygame_scale
        half_width = width/2 * self.pygame_scale
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

        # update cell
        self.cells[cell_id] = (new_body, new_shape)

    def add_cell_from_center(self, cell_id, width, length, mass, center_position, angle, angular_velocity=0.0):
        '''
        add a cell to the physics engine by specifying a dictionary with values:
            - cell_id (str): the cell's id
            - width (float)
            - length (float)
            - mass (float)
            - center_position (float)
            - angle (float)
            - angular_velocity (float) -- this value is optional, if it is not given, it is set to 0.
        '''

        half_length = length/2 * self.pygame_scale
        half_width = width/2 * self.pygame_scale

        shape = pymunk.Poly(None, (
            (-half_length, -half_width),
            (half_length, -half_width),
            (half_length, half_width),
            (-half_length, half_width)))

        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = (center_position[0] * self.pygame_scale, center_position[1] * self.pygame_scale)
        body.angle = angle
        body.dimensions = (width, length)
        body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add cell
        self.cells[cell_id] = (body, shape)

    def remove_cell(self, cell_id):
        body, shape = self.cells[cell_id]
        self.space.remove(body, shape)
        del self.cells[cell_id]

    def get_front(self, cell_id):
        body, shape = self.cells[cell_id]
        width, length = body.dimensions
        corner_position = body.position
        angle = body.angle
        front_position = front_from_corner(
            width,
            length,
            corner_position,
            angle)
        return np.array([front_position[0] / self.pygame_scale, front_position[1] / self.pygame_scale, angle])

    def get_center(self, cell_id):
        body, shape = self.cells[cell_id]
        center_position = body.position
        angle = body.angle
        return np.array([center_position[0] / self.pygame_scale, center_position[1] / self.pygame_scale, angle])

    def get_corner(self, cell_id):
        body, shape = self.cells[cell_id]
        width, length = body.dimensions
        center_position = body.position
        angle = body.angle
        corner_position = corner_from_center(
            width,
            length,
            center_position,
            angle)
        return np.array([corner_position[0] / self.pygame_scale, corner_position[1] / self.pygame_scale, angle])

    def add_barriers(self, bounds):
        """ Create static barriers """
        x_bound = bounds[0] * self.pygame_scale
        y_bound = bounds[1] * self.pygame_scale

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



# testing functions
if __name__ == '__main__':

    total_time = 100
    time_step = 1

    bounds = [10.0, 10.0]
    jitter_force = 5.0e0

    # make physics instance
    physics = MultiCellPhysics(
        bounds,
        jitter_force,
        True)

    # initialize n agents
    n_agents = 3
    agent_ids = list(range(n_agents))

    # agent initial conditions
    growth = 0.1
    volume = 1.0
    width = 0.5
    length = 2.0
    cell_density = 1100
    mass = volume * cell_density

    # add cells
    for agent_id in agent_ids:
        position = np.array([
            np.random.uniform(0, bounds[0]),
            np.random.uniform(0, bounds[1])])
        angle = np.random.uniform(0, 2 * PI)

        physics.add_cell_from_center(
            cell_id=agent_id,
            width=width,
            length=length,
            mass=mass,
            center_position=position,
            angle=angle)

    # run simulation
    time = 0
    while time < total_time:
        time += time_step

        for agent_id in agent_ids:
            (body, shape) = physics.cells[agent_id]
            (width, length) = body.dimensions
            length += growth
            physics.update_cell(agent_id, length, width, mass)
            physics.run_incremental(time_step)
