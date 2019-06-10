from __future__ import absolute_import, division, print_function

__version__ = "$Id:$"
__docformat__ = "reStructuredText"

# Library imports
import pygame
from pygame.key import *
from pygame.locals import *
from pygame.color import *

# Python imports
import random
import math
import numpy as np

# pymunk imports
import pymunk
import pymunk.pygame_util

INITIAL_MASS = 10
RADIUS = 13
INITIAL_LENGTH = 50
DIVISION_LENGTH = 3  #3
DIVISION_RANGE = [-0.25, 0.25]

ELASTICITY = 0.95
FRICTION = 0.9

PI = math.pi

pymunk_scale = 60

class MultiCellPhysics(object):
    """
    This class implements a simple scene in which there is a static platform (made up of a couple of lines)
    that don't move. Balls appear occasionally and drop onto the platform. They bounce around.
    """
    def __init__(self, bounds, translation_jitter, rotation_jitter):

        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.translation_jitter = translation_jitter * pymunk_scale
        self.rotation_jitter = rotation_jitter

        # Space
        self.space = pymunk.Space()

        # Physics
        # Time step
        self.timestep = 1
        self.physics_steps_per_frame = 60
        self.physics_dt = self.timestep / self.physics_steps_per_frame

        # self._dt = 1.0 / 60.0
        # # Number of physics steps per screen frame
        # self._steps_per_frame = 5

        # pygame
        pygame.init()
        self._screen = pygame.display.set_mode((700, 700))
        self._clock = pygame.time.Clock()
        self._draw_options = pymunk.pygame_util.DrawOptions(self._screen)

        # Static barrier walls (lines) that the balls bounce off of
        self.add_barriers(bounds)

        # objects that exist in the world
        self.cells = {}

        # Execution control and time until the next ball spawns
        self._running = True
        self._ticks_to_next_grow = 2


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

    def add_cell(self, cell_id, radius, length, mass, position, angle, angular_velocity=None):
        # shape = pymunk.Poly(None, ((0, 0), (2*radius, 0), (2*radius, length), (0, length)))
        shape = pymunk.Poly(None, (
            (0, 0),
            (length * pymunk_scale, 0),
            (length * pymunk_scale, radius * 2 * pymunk_scale),
            (0, radius * 2 * pymunk_scale)))

        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = (position[0] * pymunk_scale, position[1] * pymunk_scale)
        body.angle = angle
        body.dimensions = (radius, length)
        if angular_velocity:
            body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add cell
        self.cells[cell_id] = (body, shape)

        # self._update_screen()

    def update_cell(self, cell_id, length, radius, mass):

        body, shape = self.cells[cell_id]

        # make shape, moment of inertia, and add a body
        new_shape = pymunk.Poly(None, (
            (0, 0),
            (length * pymunk_scale, 0),
            (length * pymunk_scale, radius * 2 * pymunk_scale),
            (0, radius * 2 * pymunk_scale)))

        inertia = pymunk.moment_for_poly(mass, new_shape.get_vertices())
        new_body = pymunk.Body(mass, inertia)
        new_shape.body = new_body

        # TODO - reposition on center?
        new_body.position = body.position
        new_body.angle = body.angle
        new_body.angular_velocity = body.angular_velocity
        new_body.dimensions = (radius, length)

        new_shape.elasticity = shape.elasticity
        new_shape.friction = shape.friction

        # swap bodies
        self.space.add(new_body, new_shape)
        self.space.remove(body, shape)

        # update cell
        self.cells[cell_id] = (new_body, new_shape)

    def remove_cell(self, cell_id):
        # get body and shape from cell_id, remove from space and from cells
        body, shape = self.cells[cell_id]

        self.space.remove(body, shape)
        del self.cells[cell_id]

    def get_position(self, cell_id):
        body, shape = self.cells[cell_id]
        position = body.position
        angle = body.angle

        return np.array([position[0] / pymunk_scale, position[1] / pymunk_scale, angle])

    def add_barriers(self, bounds):
        """
        Create the static barriers.
        """

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




# For testing.

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
    length = volume_to_length(volume, radius)
    mass = volume * cell_density  # TODO -- get units to work

    # add cells
    n_cells = 5
    for cell_id in xrange(n_cells):
        position = (random.uniform(length, bounds[0] - length), random.uniform(length, bounds[1] - length))
        angle = random.uniform(0, 2 * PI)

        physics.add_cell(
            cell_id,
            radius,
            length,
            mass,
            position,
            angle,
        )

    growth = 0.1
    division_length = length * 2

    running = True
    while running:

        for cell_id, cell in physics.cells.iteritems():
            body, shape = cell
            radius, length = body.dimensions
            mass = body.mass  # TODO -- update mass

            # grow
            length += growth

            physics.update_cell(cell_id, length, radius, mass)

            # if length >= division_length:


        physics.run_incremental(5)
