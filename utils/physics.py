from __future__ import absolute_import, division, print_function

# Python imports
import random
import math
import numpy as np

# pymunk imports
import pymunk
# import pymunk.pygame_util  # TODO -- is this needed?

ELASTICITY = 0.95
FRICTION = 0.9


class MultiCellPhysics(object):
    ''''''
    def __init__(self, bounds, translation_jitter, rotation_jitter):

        self.elasticity = ELASTICITY
        self.friction = FRICTION
        self.translation_jitter = translation_jitter
        self.rotation_jitter = rotation_jitter
        
        # Space
        self.space = pymunk.Space()

        # Time step
        self.timestep = 1
        self.physics_steps_per_frame = 60
        self.physics_dt = self.timestep / self.physics_steps_per_frame

        # Static barriers
        self.add_barriers(bounds)

        # Objects that exist in the world
        self.cells = {}


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


    def add_cell(self, cell_id, radius, length, mass, position, angle, angular_velocity=None):
        shape = pymunk.Poly(None, ((0, 0), (2*radius, 0), (2*radius, length), (0, length)))
        inertia = pymunk.moment_for_poly(mass, shape.get_vertices())
        body = pymunk.Body(mass, inertia)
        shape.body = body

        body.position = position
        body.angle = angle
        body.dimensions = (2*radius, length)
        if angular_velocity:
            body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add cell
        self.cells[cell_id] = (body, shape)

    def get_position(self, cell_id):
        body, shape = self.cells[cell_id]
        position = body.position
        angle = body.angle

        return np.array([position[0], position[1], angle])

    def update_cell(self, cell_id, length, radius, mass):

        body, shape = self.cells[cell_id]

        # make shape, moment of inertia, and add a body
        new_shape = pymunk.Poly(None, ((0, 0), (radius*2, 0), (radius*2, length), (0, length)))
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


    def add_barriers(self, bounds):
        """
        Create the static barriers.
        """

        x_bound = bounds[0]
        y_bound = bounds[0]

        static_body = self.space.static_body
        static_lines = [
            pymunk.Segment(static_body, (0.0, 0.0), (x_bound, 0.0), 0.0),
            pymunk.Segment(static_body, (x_bound, 0.0), (x_bound, y_bound), 0.0),
            pymunk.Segment(static_body, (x_bound, y_bound), (0.0, y_bound), 0.0),
            pymunk.Segment(static_body, (0.0, y_bound), (0.0, 0.0), 0.0),
            ]
        for line in static_lines:
            line.elasticity = 0.95
            line.friction = 0.9
        self.space.add(static_lines)
