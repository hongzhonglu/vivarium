
__version__ = "$Id:$"
__docformat__ = "reStructuredText"

# Python imports
import random
import math

# Library imports
import pygame
from pygame.key import *
from pygame.locals import *
from pygame.color import *

# pymunk imports
import pymunk
import pymunk.pygame_util

jitter_force_sigma = 500
jitter_location_sigma = 10

INITIAL_MASS = 10
RADIUS = 25
INITIAL_LENGTH = 50
DIVISION_LENGTH = 100
ELASTICITY = 0.95
FRICTION = 0.9

class MultiCell(object):
    """
    This class implements a simple scene in which there is a static platform (made up of a couple of lines)
    that don't move. Balls appear occasionally and drop onto the platform. They bounce around.
    """
    def __init__(self):
        # Space
        self._space = pymunk.Space()
        # self._space.gravity = (0.0, -900.0)

        # Physics
        # Time step
        self._dt = 1.0 / 60.0
        # Number of physics steps per screen frame
        self._physics_steps_per_frame = 1

        # pygame
        pygame.init()
        self._screen = pygame.display.set_mode((600, 600))
        self._clock = pygame.time.Clock()

        self._draw_options = pymunk.pygame_util.DrawOptions(self._screen)

        # Static barrier walls (lines) that the balls bounce off of
        self._add_static_scenery()

        # Balls that exist in the world
        self._balls = []

        # Execution control and time until the next ball spawns
        self._running = True
        self._ticks_to_next_ball = 10
        self._ticks_to_next_grow = 10

    def run(self):
        """
        The main loop of the game.
        :return: None
        """

        # create first cell
        self._add_cell(
            RADIUS,
            INITIAL_LENGTH,
            INITIAL_MASS,
            (250, 250),
            0,
            ELASTICITY,
            FRICTION,
        )

        # Main loop
        while self._running:
            # Progress time forward
            for x in range(self._physics_steps_per_frame):
                for body in self._space.bodies:
                    # Add jitter to cells
                    force = (random.normalvariate(0, jitter_force_sigma), random.normalvariate(0, jitter_force_sigma))
                    location = (
                    random.normalvariate(0, jitter_location_sigma), random.normalvariate(0, jitter_location_sigma))
                    body.apply_force_at_local_point(force, location)

                self._space.step(self._dt)

            self._process_events()
            # self._update_cells()
            self._grow_cells()
            self._divide_cells()
            self._clear_screen()
            self._draw_objects()
            pygame.display.flip()
            # Delay fixed time between frames
            self._clock.tick(50)
            pygame.display.set_caption("fps: " + str(self._clock.get_fps()))


    def _grow_cells(self):
        self._ticks_to_next_grow -= 1
        if self._ticks_to_next_grow <= 0:
            for body in self._space.bodies:
                shape = list(body.shapes)[0]  # assumes only one shape in each body
                self.grow(body, shape)

            self._ticks_to_next_grow = 100

    def _divide_cells(self):

        for body in self._space.bodies:
            radius, length = body.dimensions
            if length > DIVISION_LENGTH:
                shape = list(body.shapes)[0]  # assumes only one shape in each body
                self.divide(body, shape)

    def divide(self, body, shape):

        radius, length = body.dimensions
        mass = body.mass

        new_length = length / 2

        pos_ratios = [0, 0.5]
        for daughter in range(2):

            dy = length * pos_ratios[daughter] * math.sin(body.angle + math.pi/2) # add rotation to correc
            dx = length * pos_ratios[daughter] * math.cos(body.angle + math.pi/2)
            position = body.position + [dx, dy]
            
            self._add_cell(
                radius,
                new_length,
                mass,
                position,
                body.angle,
                shape.elasticity,
                shape.friction,
                body.angular_velocity,
            )

        self._space.remove(body, shape)
        self._balls.remove(shape)

    def _add_cell(self, radius, length, mass, position, angle, elasticity, friction, angular_velocity=None):
        shape = pymunk.Poly(None, ((0, 0), (radius, 0), (radius, length), (0, length)))
        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = position
        body.angle = angle
        body.dimensions = (radius, length)
        if angular_velocity:
            body.angular_velocity = angular_velocity

        shape.elasticity = elasticity
        shape.friction = friction

        # add body and shape to space
        self._space.add(body, shape)
        self._balls.append(shape)


    def grow(self, body, shape):

        radius, length = body.dimensions

        mass = body.mass  # TODO -- need to update mass as well!
        length += 10  # Growth. TODO -- Pass this in.

        # make shape, moment of inertia, and add a body
        new_shape = pymunk.Poly(None, ((0, 0), (radius, 0), (radius, length), (0, length)))
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
        self._space.add(new_body, new_shape)
        self._space.remove(body, shape)
        self._balls.append(new_shape)
        self._balls.remove(shape)

    def _add_static_scenery(self):
        """
        Create the static bodies.
        :return: None
        """
        static_body = self._space.static_body
        static_lines = [
            pymunk.Segment(static_body, (10.0, 10.0), (590.0, 10.0), 0.0),
            pymunk.Segment(static_body, (590.0, 10.0), (590.0, 590.0), 0.0),
            pymunk.Segment(static_body, (590.0, 590.0), (10.0, 590.0), 0.0),
            pymunk.Segment(static_body, (10.0, 590.0), (10.0, 10.0), 0.0),
            ]
        for line in static_lines:
            line.elasticity = 0.95
            line.friction = 0.9
        self._space.add(static_lines)

    def _process_events(self):
        """
        Handle game and events like keyboard input. Call once per frame only.
        :return: None
        """
        for event in pygame.event.get():
            if event.type == QUIT:
                self._running = False
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                self._running = False
            elif event.type == KEYDOWN and event.key == K_p:
                pygame.image.save(self._screen, "bouncing_balls.png")

    def _clear_screen(self):
        """
        Clears the screen.
        :return: None
        """
        self._screen.fill(THECOLORS["white"])

    def _draw_objects(self):
        """
        Draw the objects.
        :return: None
        """
        self._space.debug_draw(self._draw_options)


if __name__ == '__main__':
    game = MultiCell()
    game.run()