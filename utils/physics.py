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

    def update_screen(self):
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
                
            self.update_screen()

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

        # self.update_screen()



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

    bounds = [10.0, 10.0]
    translation_jitter = 5.0  # 1.0
    rotation_jitter = 10000.0  # 20.0
    physics = MultiCellPhysics(bounds, translation_jitter, rotation_jitter)

    cell_density = 1100
    volume = 1
    radius = 0.5
    length = volume_to_length(volume, radius)
    mass = volume * cell_density  # TODO -- get units to work

    position = (random.uniform(length, bounds[0] - length), random.uniform(length, bounds[1] - length))
    angle = random.uniform(0, 2 * PI)

    cell_id = 1
    # create first cell
    physics.add_cell(
        cell_id,
        radius,
        length,
        mass,
        position,
        angle,
    )

    running = True
    while running:
        physics.run_incremental(5)

    # physics.run()