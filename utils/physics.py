from __future__ import absolute_import, division, print_function

__version__ = "$Id:$"
__docformat__ = "reStructuredText"

# Library imports
import pygame
from pygame.locals import *
from pygame.color import *

# Python imports
import random
import math
import numpy as np

# pymunk imports
import pymunk
import pymunk.pygame_util

pymunk_scale = 60

PI = math.pi

ELASTICITY = 0.95
FRICTION = 0.9


class MultiCellPhysics(object):
    ''''''
    def __init__(self, bounds, translation_jitter, rotation_jitter):

        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.translation_jitter = translation_jitter * pymunk_scale
        self.rotation_jitter = rotation_jitter

        # Space
        self.space = pymunk.Space()

        # Physics
        self.timestep = 1
        self.physics_steps_per_frame = 60
        self.physics_dt = self.timestep / self.physics_steps_per_frame

        # pygame
        pygame.init()
        self._screen = pygame.display.set_mode((700, 700))
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
        # pygame
        self._process_events()
        self._clear_screen()
        self._draw_objects()
        pygame.display.flip()
        # Delay fixed time between frames
        self._clock.tick(5)


    def run_incremental(self, run_for):
        time = 0
        while time < run_for:
            time += self.timestep

            # Progress time forward
            for x in range(self.physics_steps_per_frame * self.timestep):
                for body in self.space.bodies:
                    # Add jitter to cells
                    force = (
                        random.normalvariate(0, self.rotation_jitter) * self.physics_dt,
                        random.normalvariate(0, self.rotation_jitter) * self.physics_dt)
                    location = (
                        random.normalvariate(0, self.translation_jitter) * self.physics_dt,
                        random.normalvariate(0, self.translation_jitter) * self.physics_dt)
                    body.apply_force_at_local_point(force, location)

                self.space.step(self.physics_dt)

            self._update_screen()


    def add_cell_from_corner(self, cell_id, width, length, mass, corner_position, angle, angular_velocity=None):

        shape = pymunk.Poly(None, (
            (0, 0),
            (length * pymunk_scale, 0),
            (length * pymunk_scale, width * pymunk_scale),
            (0, width * pymunk_scale)))

        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = (corner_position[0] * pymunk_scale, corner_position[1] * pymunk_scale)
        body.angle = angle
        body.dimensions = (width, length)
        if angular_velocity:
            body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add cell
        self.cells[cell_id] = (body, shape)

    def add_cell_from_center(self, cell_id, width, length, mass, center_position, angle, angular_velocity=None):
        corner_position = self.corner_from_center(width, length, center_position, angle)
        self.add_cell_from_corner(
            cell_id,
            width,
            length,
            mass,
            corner_position,
            angle,
            angular_velocity,
        )

    def center_from_corner(self, width, length, corner_position, angle):
        half_length = length/2
        half_width = width/2
        dx = half_length * math.cos(angle) + half_width * math.cos(angle + PI/2)
        dy = half_length * math.sin(angle) + half_width * math.sin(angle + PI/2)
        center_position = [corner_position[0] + dx, corner_position[1] + dy]

        return center_position

    def corner_from_center(self, width, length, center_position, angle):
        half_length = length/2
        half_width = width/2
        dx = half_length * math.cos(angle) + half_width * math.cos(angle + PI/2)
        dy = half_length * math.sin(angle) + half_width * math.sin(angle + PI/2)
        corner_position = [center_position[0] - dx, center_position[1] - dy]

        return corner_position

    def update_cell(self, cell_id, length, width, mass):

        body, shape = self.cells[cell_id]

        width_0, length_0 = body.dimensions
        d_width = width - width_0
        d_length = length - length_0

        # make shape, moment of inertia, and add a body
        new_shape = pymunk.Poly(None, (
            (0, 0),
            (length * pymunk_scale, 0),
            (length * pymunk_scale, width * pymunk_scale),
            (0, width * pymunk_scale)))

        inertia = pymunk.moment_for_poly(mass, new_shape.get_vertices())
        new_body = pymunk.Body(mass, inertia)
        new_shape.body = new_body

        # reposition on center
        dx = -1 * d_width/2 * math.cos(body.angle)
        dy = -1 * d_length/2 * math.sin(body.angle)

        new_body.position = body.position + [dx * pymunk_scale, dy * pymunk_scale]
        new_body.angle = body.angle
        new_body.angular_velocity = body.angular_velocity
        new_body.dimensions = (width, length)

        new_shape.elasticity = shape.elasticity
        new_shape.friction = shape.friction

        # swap bodies
        self.space.add(new_body, new_shape)
        self.space.remove(body, shape)

        # update cell
        self.cells[cell_id] = (new_body, new_shape)


    def divide(self, cell_id, daughter_ids):

        body, shape = self.cells[cell_id]
        width, length = body.dimensions
        mass = body.mass # TODO -- divide mass?
        new_length = length / 2

        pos_ratios = [0, 0.5]
        for index, daughter_id in enumerate(daughter_ids):
            dx = length * pos_ratios[index] * math.cos(body.angle)
            dy = length * pos_ratios[index] * math.sin(body.angle)

            position = body.position/pymunk_scale + [dx, dy]

            self.add_cell_from_corner(
                daughter_id,
                width,
                new_length,
                mass,
                position,
                body.angle,
                body.angular_velocity,
            )

            self._update_screen()

        self.remove_cell(cell_id)
        # TODO -- return positions? new cell_ids?

    def remove_cell(self, cell_id):
        body, shape = self.cells[cell_id]
        self.space.remove(body, shape)
        del self.cells[cell_id]

    def get_center(self, cell_id):
        body, shape = self.cells[cell_id]
        width, length = body.dimensions
        corner_position = body.position
        angle = body.angle

        # get center
        center_position = self.center_from_corner(width, length, corner_position, angle)

        return np.array([center_position[0] / pymunk_scale, center_position[1] / pymunk_scale, angle])

    def get_corner(self, cell_id):
        body, shape = self.cells[cell_id]
        corner_position = body.position
        angle = body.angle

        return np.array([corner_position[0] / pymunk_scale, corner_position[1] / pymunk_scale, angle])

    def add_barriers(self, bounds):
        """ Create static barriers """

        x_bound = bounds[0] * pymunk_scale
        y_bound = bounds[1] * pymunk_scale

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










# For testing. Not used in object.

def volume_to_length(volume, radius):
    '''
    get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
    a=length of cylinder without rounded caps, l=total length:

    V = (4/3)*PI*r^3 + PI*r^2*a
    l = a + 2*r
    '''

    PI = math.pi

    cylinder_length = (volume - (4 / 3) * PI * radius ** 3) / (PI * radius ** 2)
    total_length = cylinder_length + 2 * radius

    return total_length

if __name__ == '__main__':
    cell_density = 1100

    bounds = [10.0, 10.0]
    translation_jitter = 5.0  # 1.0
    rotation_jitter = 50000.0  # 20.0
    physics = MultiCellPhysics(bounds, translation_jitter, rotation_jitter)

    volume = 1
    radius = 0.5
    width = radius * 2
    length = volume_to_length(volume, radius)
    mass = volume * cell_density  # TODO -- get units to work

    # add cells
    n_cells = 5
    for cell_id in xrange(n_cells):
        position = (random.uniform(length, bounds[0] - length), random.uniform(length, bounds[1] - length))
        angle = random.uniform(0, 2 * PI)

        physics.add_cell_from_corner(
            cell_id,
            width,
            length,
            mass,
            position,
            angle,
        )

    growth = 0.1
    division_length = length * 2

    running = True
    while running:

        for cell_id in physics.cells.keys():
            body, shape = physics.cells[cell_id]
            width, length = body.dimensions
            mass = body.mass  # TODO -- update mass

            # grow
            length += growth

            if length >= division_length:
                max_cell_id = max(physics.cells.keys())
                daughter_ids = [max_cell_id + 1, max_cell_id + 2]
                physics.divide(cell_id, daughter_ids)
            else:
                physics.update_cell(cell_id, length, width, mass)

        physics.run_incremental(5)
