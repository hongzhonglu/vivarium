"""This example spawns (bouncing) balls randomly on a L-shape constructed of
two segment shapes. Not interactive.
"""

__version__ = "$Id:$"
__docformat__ = "reStructuredText"

# Python imports
import random

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

class BouncyBalls(object):
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
            self._update_cells()
            self._grow_cells()
            self._divide_cells()
            self._clear_screen()
            self._draw_objects()
            pygame.display.flip()
            # Delay fixed time between frames
            self._clock.tick(50)
            pygame.display.set_caption("fps: " + str(self._clock.get_fps()))

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

    def _update_cells(self):
        """
        Create/remove balls as necessary. Call once per frame only.
        :return: None
        """
        self._ticks_to_next_ball -= 1
        if self._ticks_to_next_ball <= 0:
            self._create_cell()
            self._ticks_to_next_ball = 100

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

        # import ipdb; ipdb.set_trace()
        length = length / 2

        for daughter in range(2):

            # TODO -- place new cells in correct positions
            new_body, new_shape = self._add_cell(
                radius,
                length,
                mass,
                body.position,
                body.angle,
                body.angular_velocity,
                shape.elasticity,
                shape.friction,
            )

            # swap bodies
            self._space.add(new_body, new_shape)
            self._balls.append(new_shape)


        self._space.remove(body, shape)
        self._balls.remove(shape)



    def _add_cell(self, radius, length, mass, position, angle, angular_velocity, elasticity, friction):
        shape = pymunk.Poly(None, ((0, 0), (radius, 0), (radius, length), (0, length)))
        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = position
        body.angle = angle
        body.angular_velocity = angular_velocity
        body.dimensions = (radius, length)

        shape.elasticity = elasticity
        shape.friction = friction

        return body, shape


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


    # TODO us _add_cell instead
    def _create_cell(self):
        """
        Create a cell.
        """
        mass = INITIAL_MASS
        radius = RADIUS
        length = INITIAL_LENGTH

        # make shape, moment of inertia, and add a body
        shape = pymunk.Poly(None, ((0, 0), (radius, 0), (radius, length), (0, length)))
        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        x = random.randint(115, 350)
        body.position = x, 400
        body.dimensions = (radius, length)

        shape.elasticity = 0.95
        shape.friction = 0.9
        self._space.add(body, shape)
        self._balls.append(shape)


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
    game = BouncyBalls()
    game.run()