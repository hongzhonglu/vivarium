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
FRICTION = 0.9


class MultiCellPhysics(object):
    ''''''
    def __init__(self, bounds, translation_jitter, rotation_jitter, debug=False):
        self.pygame_scale = 700 / bounds[0]
        self.pygame_viz = debug
        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.translation_jitter = translation_jitter
        self.rotation_jitter = rotation_jitter

        # Space
        self.space = pymunk.Space()

        # Physics
        self.timestep = 1
        self.physics_steps_per_frame = 10
        self.physics_dt = self.timestep / self.physics_steps_per_frame

        if self.pygame_viz:
            pygame.init()
            self._screen = pygame.display.set_mode((710, 710))
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

    def random_body_position(self, body):
        ''' pick a random point along the boundary'''
        width, length = body.dimensions
        if random.randint(0, 1) == 0:
            # force along ends
            if random.randint(0, 1) == 0:
                # force on the left end
                location = (0, random.uniform(0, width))
            else:
                # force on the right end
                location = (length, random.uniform(0, width))
        else:
            # force along length
            if random.randint(0, 1) == 0:
                # force on the bottom end
                location = (random.uniform(0, length), 0)
            else:
                # force on the top end
                location = (random.uniform(0, length), width)

        return location

    def apply_motile_force(self, cell_id, force, torque):
        body, shape = self.cells[cell_id]
        body.motile_force = (force, torque)

    def run_incremental(self, run_for):
        time = 0
        while time < run_for:
            time += self.timestep

            # Progress time forward
            for x in range(self.physics_steps_per_frame * self.timestep):
                for body in self.space.bodies:
                    width, length = body.dimensions

                    # random jitter
                    jitter_torque = random.normalvariate(0, self.rotation_jitter)
                    jitter_force = [
                        random.normalvariate(0, self.translation_jitter),
                        random.normalvariate(0, self.translation_jitter)]
                    location = (length/2, width/2)  #self.random_body_position(body)

                    # motile forces
                    motile_torque = 0.0
                    motile_force = [0.0, 0.0]
                    if hasattr(body, 'motile_force'):
                        force, motile_torque = body.motile_force
                        motile_force = [force, 0.0]  # force is applied in the positive x-direction (forward)

                    body.angular_velocity = (jitter_torque + motile_torque)  # TODO (eran) add to angular velocity rather than replace it. Needs better damping first
                    total_force = [a + b for a, b in zip(jitter_force, motile_force)]
                    body.apply_force_at_local_point(total_force, location)

                self.space.step(self.physics_dt)

            # Disable momentum at low Reynolds number (the ratio of inertial and viscous forces)
            for body in self.space.bodies:
                body.velocity -= (0.5 * body.velocity * self.timestep)
                body.angular_velocity -= (0.5 * body.angular_velocity * self.timestep)  # TODO (Eran) this should be function of viscosity

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
        front_position = self.front_from_corner(width*self.pygame_scale, length*self.pygame_scale, corner_position, angle)
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
        corner_position = self.corner_from_center(
            width * self.pygame_scale,
            length * self.pygame_scale,
            center_position,
            angle)
        return np.array([corner_position[0] / self.pygame_scale, corner_position[1] / self.pygame_scale, angle])

    def front_from_corner(self, width, length, corner_position, angle):
        half_width = width/2
        dx = length * math.cos(angle) + half_width * math.cos(angle + PI/2)  # PI/2 gives a half-rotation for the width component
        dy = length * math.sin(angle) + half_width * math.sin(angle + PI/2)
        front_position = [corner_position[0] + dx, corner_position[1] + dy]
        return np.array([front_position[0], front_position[1], angle])

    def corner_from_center(self, width, length, center_position, angle):
        half_length = length/2
        half_width = width/2
        dx = half_length * math.cos(angle) + half_width * math.cos(angle + PI/2)
        dy = half_length * math.sin(angle) + half_width * math.sin(angle + PI/2)
        corner_position = [center_position[0] - dx, center_position[1] - dy]

        return np.array([corner_position[0], corner_position[1], angle])

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
def set_motile_force(physics, agent_id, object_id):
    body, shape = physics.cells[agent_id]
    obj_body, obj_shape = physics.cells[object_id]

    obj_position = physics.get_center(object_id)  # obj_body.position
    position = physics.get_front(agent_id)
    angle = body.angle

    obj_distance = obj_position - position
    obj_angle = math.atan2(obj_distance[1],obj_distance[0])
    obj_relative_angle = obj_angle - angle

    # run/tumble
    if abs(obj_relative_angle) < PI/4:
        # run
        force = 15000.0
        torque = 0.0
        print('RUN!')
    else:
        # tumble
        force = 5000.0
        torque = random.normalvariate(0, 0.02)  #random.uniform(-0.005, 0.005)  #random.uniform(-50000.0, 50000.0)  # random.normalvariate(0, 5000.0) #5000.0 #random.uniform(0.0, 2*PI)
        print('TUMBLE!')

    physics.apply_motile_force(agent_id, force, torque)


# For testing with pygame
if __name__ == '__main__':

    bounds = [20.0, 20.0]

    agent_id = 1
    volume = 1.0
    width = 0.5
    length = 2.0
    cell_density = 1100
    mass = volume * cell_density
    translation_jitter = 0.5
    rotation_jitter = 0.005

    position = (2.0, 2.0)
    angle = PI/2

    # make physics instance
    physics = MultiCellPhysics(
        bounds,
        translation_jitter,
        rotation_jitter,
        True)

    # add cell
    physics.add_cell_from_center(
        cell_id = agent_id,
        width = width,
        length = length,
        mass = mass,
        position = position,
        angle = angle,
    )

    # add object
    object_id = 999
    physics.add_cell_from_center(
        cell_id=object_id,    # agent_id,
        width=0.5,          # width
        length=0.5,          # length
        mass=100000.0,       # mass
        position=(14.0, 14.0),     # position
        angle=0.0,          # angle
    )

    running = True
    growth = 0.1  # 0.02
    while running:
        length += growth
        physics.update_cell(agent_id, length, width, mass)
        # set_motile_force(physics, agent_id, object_id)
        physics.run_incremental(5)
