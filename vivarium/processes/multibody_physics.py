from __future__ import absolute_import, division, print_function

import os
import argparse
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import random
import math

import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb
from matplotlib.collections import LineCollection

# pymunk imports
import pymunkoptions
pymunkoptions.options["debug"] = False
import pymunk
import pymunk.pygame_util

# pygame for debugging
import pygame
from pygame.key import *
from pygame.locals import *
from pygame.color import *

# vivarium imports
from vivarium.compartment.emitter import timeseries_from_data
from vivarium.compartment.process import (
    Process,
    COMPARTMENT_STATE)
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_process,
    simulate_compartment)
from vivarium.processes.Vladimirov2008_motor import run, tumble
from vivarium.processes.derive_globals import (
    volume_from_length)



DEBUG_SIZE = 600  # size of the pygame debug screen
DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi

# colors for phylogeny initial agents
HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_HUE = HUES[0]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]

# agent port keys
AGENT_KEYS = ['location', 'angle', 'volume', 'length', 'width', 'mass', 'forces']
NON_AGENT_KEYS = ['fields', 'time', 'global', COMPARTMENT_STATE]





def random_body_position(body):
    # pick a random point along the boundary
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

def daughter_locations(parent_location, parent_length, parent_angle):
    pos_ratios = [-0.25, 0.25]
    daughter_locations = []
    for daughter in range(2):
        dx = parent_length * pos_ratios[daughter] * math.cos(parent_angle)
        dy = parent_length * pos_ratios[daughter] * math.sin(parent_angle)
        location = [parent_location[0] + dx, parent_location[1] + dy]
        daughter_locations.append(location)

    return daughter_locations



class Multibody(Process):
    """
    A multi-body physics process using pymunk

    Notes:
    - rotational diffusion in liquid medium with viscosity = 1 mPa.s: Dr = 3.5+/-0.3 rad^2/s
        (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
    - translational diffusion in liquid medium with viscosity = 1 mPa.s: Dt=100 micrometers^2/s
        (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)

    """

    defaults = {
        'initial_agents': {},
        'elasticity': 0.9,
        'damping': 0.05, # simulates viscous forces to reduce velocity at low Reynolds number (1 = no damping, 0 = full damping)
        'angular_damping': 0.7,  # less damping for angular velocity seems to improve behavior
        'friction': 0.9,  # TODO -- does this do anything?
        'physics_dt': 0.005,
        'force_scaling': 100,  # scales from pN
        'jitter_force': 1e-3,  # pN
        'bounds': DEFAULT_BOUNDS,
        'mother_machine': False,
        'animate': False,
        'debug': False,
    }


    def __init__(self, initial_parameters={}):

        # hardcoded parameters
        self.elasticity = self.defaults['elasticity']
        self.friction = self.defaults['friction']
        self.damping = self.defaults['damping']
        self.angular_damping = self.defaults['angular_damping']
        self.force_scaling = self.defaults['force_scaling']
        self.physics_dt = self.defaults['physics_dt']

        # configured parameters
        self.jitter_force = initial_parameters.get('jitter_force', self.defaults['jitter_force'])
        self.bounds = initial_parameters.get('bounds', self.defaults['bounds'])

        # initialize pymunk space
        self.space = pymunk.Space()

        # debug screen with pygame
        self.pygame_viz = initial_parameters.get('debug', self.defaults['debug'])
        self.pygame_scale = 1  # pygame_scale scales the debug screen
        if self.pygame_viz:
            max_bound = max(self.bounds)
            self.pygame_scale = DEBUG_SIZE / max_bound
            self.force_scaling *= self.pygame_scale
            pygame.init()
            self._screen = pygame.display.set_mode((
                int(self.bounds[0]*self.pygame_scale),
                int(self.bounds[1]*self.pygame_scale)), RESIZABLE)
            self._clock = pygame.time.Clock()
            self._draw_options = pymunk.pygame_util.DrawOptions(self._screen)

        # add static barriers
        # TODO -- mother machine configuration
        self.mother_machine = initial_parameters.get('mother_machine', self.defaults['mother_machine'])
        self.add_barriers(self.bounds)

        # initialize agents
        self.agent_bodies = {}
        self.initial_agents = initial_parameters.get('agents', self.defaults['initial_agents'])
        for agent_id, specs in self.initial_agents.items():
            self.add_body_from_center(agent_id, specs)

        # interactive plot for visualization
        self.animate = initial_parameters.get('animate', self.defaults['animate'])
        if self.animate:
            plt.ion()
            self.ax = plt.gca()
            self.ax.set_aspect('equal')
            self.animate_frame(self.initial_agents)

        # all initial agents get a key under a single port
        ports = {'agents': ['agents']}

        parameters = {}
        parameters.update(initial_parameters)

        super(Multibody, self).__init__(ports, parameters)

    def default_settings(self):

        state = {'agents': {'agents': self.initial_agents}}

        schema = {'agents': {'agents': {'updater': 'merge'}}}

        default_emitter_keys = {
            port_id: keys for port_id, keys in self.ports.items()}

        return {
            'state': state,
            'schema': schema,
            'emitter_keys': default_emitter_keys,
            'time_step': 2
        }

    def next_update(self, timestep, states):
        agents = states['agents']['agents']

        # animate before update
        if self.animate:
            self.animate_frame(agents)

        # if an agent has been removed from the agents store,
        # remove it from space and agent_bodies
        removed_agents = [
            agent_id for agent_id in self.agent_bodies.keys()
            if agent_id not in agents.keys()]
        for agent_id in removed_agents:
            body, shape = self.agent_bodies[agent_id]
            self.space.remove(body, shape)
            del self.agent_bodies[agent_id]

        # update agents, add new agents
        for agent_id, specs in agents.items():
            if agent_id in self.agent_bodies:
                self.update_body(agent_id, specs)
            else:
                self.add_body_from_center(agent_id, specs)

        # run simulation
        self.run(timestep)

        # get new agent position
        agent_position = {
            agent_id: self.get_body_position(agent_id)
            for agent_id in self.agent_bodies.keys()}

        return {'agents': {'agents': agent_position}}

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

        if self.pygame_viz:
            self._update_screen()

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
            # force-based torque
            # if torque != 0.0:
            #     motile_force = get_force_with_angle(thrust, torque)

        scaled_motile_force = [thrust * self.force_scaling for thrust in motile_force]
        body.apply_force_at_local_point(scaled_motile_force, motile_location)

    def apply_jitter_force(self, body):
        jitter_location = random_body_position(body)
        jitter_force = [
            random.normalvariate(0, self.jitter_force),
            random.normalvariate(0, self.jitter_force)]
        scaled_jitter_force = [
            force * self.force_scaling
            for force in jitter_force]
        body.apply_force_at_local_point(
            scaled_jitter_force,
            jitter_location)

    def apply_viscous_force(self, body):
        # dampen the velocity
        body.velocity = body.velocity * self.damping + (body.force / body.mass) * self.physics_dt
        body.angular_velocity = body.angular_velocity * self.angular_damping + body.torque / body.moment * self.physics_dt

    def add_barriers(self, bounds):
        """ Create static barriers """
        thickness = 0.2

        x_bound = bounds[0] * self.pygame_scale
        y_bound = bounds[1] * self.pygame_scale

        static_body = self.space.static_body
        static_lines = [
            pymunk.Segment(static_body, (0.0, 0.0), (x_bound, 0.0), thickness),
            pymunk.Segment(static_body, (x_bound, 0.0), (x_bound, y_bound), thickness),
            pymunk.Segment(static_body, (x_bound, y_bound), (0.0, y_bound), thickness),
            pymunk.Segment(static_body, (0.0, y_bound), (0.0, 0.0), thickness),
        ]

        if self.mother_machine:
            channel_height = self.mother_machine.get('channel_height') * self.pygame_scale
            channel_space = self.mother_machine.get('channel_space') * self.pygame_scale

            n_lines = math.floor(x_bound/channel_space)

            machine_lines = [
                pymunk.Segment(
                    static_body,
                    (channel_space * line, 0),
                    (channel_space * line, channel_height), thickness)
                for line in range(n_lines)]
            static_lines += machine_lines

        for line in static_lines:
            line.elasticity = 0.0  # no bounce
            line.friction = 0.9
        self.space.add(static_lines)

    def add_body_from_center(self, body_id, body):
        width = body['width'] * self.pygame_scale
        length = body['length'] * self.pygame_scale
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

        body.position = (
            center_position[0] * self.pygame_scale,
            center_position[1] * self.pygame_scale)
        body.angle = angle
        body.dimensions = (width, length)
        body.angular_velocity = angular_velocity

        shape.elasticity = self.elasticity
        shape.friction = self.friction

        # add body and shape to space
        self.space.add(body, shape)

        # add body to agents dictionary
        self.agent_bodies[body_id] = (body, shape)

    def update_body(self, body_id, specs):

        length = specs['length'] * self.pygame_scale
        width = specs['width'] * self.pygame_scale
        mass = specs['mass']
        motile_force = specs.get('motile_force', [0, 0])

        body, shape = self.agent_bodies[body_id]
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
        new_body.motile_force = motile_force

        new_shape.elasticity = shape.elasticity
        new_shape.friction = shape.friction

        # swap bodies
        self.space.remove(body, shape)
        self.space.add(new_body, new_shape)

        # update body
        self.agent_bodies[body_id] = (new_body, new_shape)

    def get_body_position(self, agent_id):
        body, shape = self.agent_bodies[agent_id]
        position = body.position
        rescaled_position = [
            position[0] / self.pygame_scale,
            position[1] / self.pygame_scale]

        # enforce bounds
        rescaled_position = [
            0 if pos<0 else pos
            for idx, pos in enumerate(rescaled_position)]
        rescaled_position = [
            self.bounds[idx] if pos>self.bounds[idx] else pos
            for idx, pos in enumerate(rescaled_position)]

        return {
            'location': rescaled_position,
            'angle': body.angle,
        }


    ## matplotlib interactive plot
    def animate_frame(self, agents):
        plt.cla()
        for agent_id, data in agents.items():
            # location, orientation, length
            x_center = data['location'][0]
            y_center = data['location'][1]
            angle = data['angle'] / PI * 180 + 90  # rotate 90 degrees to match field
            length = data['length']
            width = data['width']

            # get bottom left position
            x_offset = (width / 2)
            y_offset = (length / 2)
            theta_rad = math.radians(angle)
            dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
            dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

            x = x_center - dx
            y = y_center - dy

            # Create a rectangle
            rect = patches.Rectangle((x, y), width, length, angle=angle, linewidth=1, edgecolor='b')
            self.ax.add_patch(rect)

        plt.xlim([0, self.bounds[0]])
        plt.ylim([0, self.bounds[1]])
        plt.draw()
        plt.pause(0.05)


    ## pygame visualization (for debugging)
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



# configs
def get_n_dummy_agents(n_agents):
    return {agent_id: None for agent_id in range(n_agents)}

def random_body_config(config):
    # cell dimensions
    width = 1
    length = 2
    volume = volume_from_length(length, width)

    n_agents = config['n_agents']
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    agents = get_n_dummy_agents(n_agents)
    agent_config = {
        agent_id: {
            'location': [
                np.random.uniform(0, bounds[0]),
                np.random.uniform(0, bounds[1])],
            'angle': np.random.uniform(0, 2 * PI),
            'volume': volume,
            'length': length,
            'width': width,
            'mass': 1,
            'forces': [0, 0]}
        for agent_id in agents.keys()}

    return {
        'agents': agent_config,
        'bounds': bounds}

def mother_machine_body_config(config):
    # cell dimensions
    width = 1
    length = 2
    volume = volume_from_length(length, width)

    n_agents = config['n_agents']
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    channel_space = config.get('channel_space', 1)

    # possible locations, shuffled for index-in
    n_spaces = math.floor(bounds[0]/channel_space)
    assert n_agents < n_spaces, 'more agents than mother machine spaces'

    possible_locations = [
        [x*channel_space - channel_space/2, 0.01]
        for x in range(1, n_spaces)]
    random.shuffle(possible_locations)

    agents = get_n_dummy_agents(n_agents)
    agent_config = {
        agent_id: {
            'location': possible_locations[index],
            'angle': PI/2,
            'volume': volume,
            'length': length,
            'width': width,
            'mass': 1,
            'forces': [0, 0]}
        for index, agent_id in enumerate(agents.keys())}

    return {
        'agents': agent_config,
        'bounds': bounds}

# tests and simulations
def test_multibody(config={'n_agents':1}, time=1):
    body_config = random_body_config(config)
    multibody = Multibody(body_config)

    # initialize agent's boundary state
    initial_agents_state = body_config['agents']
    compartment = process_in_compartment(multibody)
    compartment.send_updates({'agents': [{'agents': initial_agents_state}]})

    settings = {
        'total_time': time,
        'return_raw_data': True,
        'environment_port': 'external',
        'environment_volume': 1e-2}

    return simulate_compartment(compartment, settings)

def simulate_motility(config, settings):
    # time of motor behavior without chemotaxis
    run_time = 0.42  # s (Berg)
    tumble_time = 0.14  # s (Berg)

    total_time = settings['total_time']
    timestep = settings['timestep']
    initial_agents_state = config['agents']

    # make the process
    multibody = Multibody(config)
    compartment = process_in_compartment(multibody)

    # initialize agent boundary state
    compartment.send_updates({'agents': [{'agents': initial_agents_state}]})

    # get initial agent state
    agents_store = compartment.states['agents']
    agents_state = agents_store.state

    # initialize hidden agent motile states, and update agent motile_forces in agent store
    agent_motile_states = {}
    motile_forces = {}
    for agent_id, specs in agents_state['agents'].items():
        motile_force = run()
        agent_motile_states[agent_id] = {
            'motor_state': 1,  # 0 for run, 1 for tumble
            'time_in_motor_state': 0}
        motile_forces[agent_id] = {'motile_force': motile_force}
    compartment.send_updates({'agents': [{'agents': motile_forces}]})

    ## run simulation
    # test run/tumble
    time = 0
    while time < total_time:
        print('time: {}'.format(time))
        time += timestep

        # update motile force and apply to state
        motile_forces = {}
        for agent_id, motile_state in agent_motile_states.items():
            motor_state = motile_state['motor_state']
            time_in_motor_state = motile_state['time_in_motor_state']

            if motor_state == 1:  # tumble
                if time_in_motor_state < tumble_time:
                    motile_force = run()
                    time_in_motor_state += timestep
                else:
                    # switch
                    motile_force = run()
                    motor_state = 0
                    time_in_motor_state = 0

            elif motor_state == 0:  # run
                if time_in_motor_state < run_time:
                    motile_force = tumble()
                    time_in_motor_state += timestep
                else:
                    # switch
                    motile_force = tumble()
                    motor_state = 1
                    time_in_motor_state = 0

            agent_motile_states[agent_id] = {
                'motor_state': motor_state,  # 0 for run, 1 for tumble
                'time_in_motor_state': time_in_motor_state}
            motile_forces[agent_id] = {'motile_force': motile_force}

        compartment.send_updates({'agents': [{'agents': motile_forces}]})
        update = compartment.update(timestep)

    return compartment.emitter.get_data()

def run_mother_machine():
    bounds = [30, 30]
    channel_height = 0.7 * bounds[1]
    channel_space = 1.5

    settings = {
        'growth_rate': 0.02,
        'growth_rate_noise': 0.02,
        'division_volume': 2.6,
        'channel_height': channel_height,
        'total_time': 240}
    mm_config = {
        'animate': True,
        'mother_machine': {
            'channel_height': channel_height,
            'channel_space': channel_space},
        'jitter_force': 2e-2,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'channel_height': channel_height,
        'channel_space': channel_space,
        'n_agents': 5}
    mm_config.update(mother_machine_body_config(body_config))
    mm_data = simulate_growth_division(mm_config, settings)

    # make snapshot
    agents = {time: time_data['agents']['agents'] for time, time_data in mm_data.items()}
    fields = {}
    plot_snapshots(agents, fields, mm_config, out_dir, 'mother_machine_snapshots')


def run_motility(out_dir):
    # test motility
    bounds = [100, 100]
    motility_sim_settings = {
        'timestep': 0.05,
        'total_time': 5}
    motility_config = {
        'animate': True,
        'jitter_force': 0,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'n_agents': 6}
    motility_config.update(random_body_config(body_config))

    # run motility sim
    motility_data = simulate_motility(motility_config, motility_sim_settings)

    # make motility plot
    reduced_data = {time: data['agents'] for time, data in motility_data.items()}
    motility_timeseries = timeseries_from_data(reduced_data)
    plot_motility(motility_timeseries, out_dir)
    plot_trajectory(motility_timeseries, motility_config, out_dir)

    # make motility snapshot
    agents = {time: time_data['agents']['agents'] for time, time_data in motility_data.items()}
    fields = {}
    plot_snapshots(agents, fields, motility_config, out_dir, 'motility_snapshots')

def run_growth_division():
    bounds = [20, 20]
    settings = {
        'growth_rate': 0.02,
        'growth_rate_noise': 0.02,
        'division_volume': 2.6,
        'total_time': 300}

    gd_config = {
        'animate': True,
        'jitter_force': 1e-1,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'n_agents': 1}
    gd_config.update(random_body_config(body_config))
    gd_data = simulate_growth_division(gd_config, settings)

    # make snapshot
    agents = {time: time_data['agents']['agents'] for time, time_data in gd_data.items()}
    fields = {}
    plot_snapshots(agents, fields, gd_config, out_dir, 'growth_division_snapshots')


def simulate_growth_division(config, settings):

    # make the process
    multibody = Multibody(config)
    compartment = process_in_compartment(multibody)

    # get initial agent state
    agents_store = compartment.states['agents']

    ## run simulation
    # get simulation settings
    growth_rate = settings.get('growth_rate', 0.0006)
    growth_rate_noise = settings.get('growth_rate_noise', 0.0)
    division_volume = settings.get('division_volume', 0.4)
    channel_height = settings.get('channel_height')
    total_time = settings.get('total_time', 10)
    timestep = compartment.time_step

    time = 0
    while time < total_time:
        print('time: {}'.format(time))
        time += timestep

        agents_state = agents_store.state
        agent_updates = {}
        remove_agents = []
        add_agents = {}
        for agent_id, state in agents_state['agents'].items():
            location = state['location']
            angle = state['angle']
            length = state['length']
            width = state['width']
            mass = state['mass']

            # update
            growth_rate2 = (growth_rate + np.random.normal(0.0, growth_rate_noise)) * timestep
            new_mass = mass + mass * growth_rate2
            new_length = length + length * growth_rate2
            new_volume = volume_from_length(new_length, width)

            if channel_height and location[1] > channel_height:
                remove_agents.append(agent_id)
            elif new_volume > division_volume:
                daughter_ids = [str(agent_id) + '0', str(agent_id) + '1']

                # daughter state with updated values
                half_mass = new_mass / 2
                half_length = new_length / 2
                new_locations = daughter_locations(location, length, angle)

                daughter_states = {}
                for index, daughter_id in enumerate(daughter_ids):
                    daughter_state = state.copy()
                    daughter_state.update({
                        'location': new_locations[index],
                        'mass': half_mass,
                        'length': half_length})
                    daughter_states[daughter_id] = daughter_state

                # remove mother from store, add daughters
                remove_agents.append(agent_id)
                add_agents.update(daughter_states)
                agent_updates.update(daughter_states)  # TODO -- why is this not updating the store?

            else:
                agent_updates[agent_id] = {
                    'volume': new_volume,
                    'length': new_length,
                    'mass': new_mass}

        for agent_id in remove_agents:
            # remove from store
            del agents_state['agents'][agent_id]

        if add_agents:
            # add to store
            agents_state['agents'].update(add_agents)

        # update compartment
        compartment.send_updates({'agents': [{'agents': agent_updates}]})
        compartment.update(timestep)

    return compartment.emitter.get_data()


# plotting
def plot_agent(ax, data, color):
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

def plot_agents(ax, agents, agent_colors={}):
    '''
    - ax: the axis for plot
    - agents: a dict with {agent_id: agent_data} and
        agent_data a dict with keys location, angle, length, width
    - agent_colors: dict with {agent_id: hsv color}
    '''
    for agent_id, agent_data in agents.items():
        color = agent_colors.get(agent_id, [DEFAULT_HUE]+DEFAULT_SV)
        plot_agent(ax, agent_data, color)

def plot_snapshots(agents, fields, config, out_dir='out', filename='snapshots'):
    '''
        - agents (dict): with {time: agent_data}
        - fields TODO
        - config (dict): the environment config for the simulation
    '''
    n_snapshots = 6
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    edge_length_x = bounds[0]
    edge_length_y = bounds[1]

    # time steps that will be used
    time_vec = list(agents.keys())
    time_indices = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
    snapshot_times = [time_vec[i] for i in time_indices]

    # get fields id and range
    field_ids = []
    if fields:
        field_ids = list(fields[time_vec[0]].keys())
        field_range = {}
        for field_id in field_ids:
            field_min = min([field_data[field_id].min() for t, field_data in fields.items()])
            field_max = max([field_data[field_id].max() for t, field_data in fields.items()])
            field_range[field_id] = [field_min, field_max]

    # get agent ids
    agent_ids = set()
    for time, time_data in agents.items():
        current_agents = list(time_data.keys())
        agent_ids.update(current_agents)
    agent_ids = list(agent_ids)

    # set agent colors
    agent_colors = {}
    for agent_id in agent_ids:
        hue = random.choice(HUES)  # select random initial hue
        color = [hue] + DEFAULT_SV
        agent_colors[agent_id] = color

    # make the figure
    n_rows = max(len(field_ids), 1)
    n_cols = n_snapshots + 1  # one column for the colorbar
    fig = plt.figure(figsize=(12 * n_cols, 12 * n_rows))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    # plot snapshot data in each subsequent column
    for col_idx, (time_idx, time) in enumerate(zip(time_indices, snapshot_times), 1):
        agents_now = agents[time]
        if field_ids:
            for row_idx, field_id in enumerate(field_ids):
                ax = init_axes(fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time)
                # transpose field to align with agent
                field = np.transpose(np.array(fields[time][field_id])).tolist()
                vmin, vmax = field_range[field_id]
                im = plt.imshow(field,
                                origin='lower',
                                extent=[0, edge_length_x, 0, edge_length_y],
                                vmin=vmin,
                                vmax=vmax,
                                cmap='BuPu')

                plot_agents(ax, agents_now, agent_colors)
        else:
            row_idx = 0
            ax = init_axes(fig, bounds[0], bounds[1], grid, row_idx, col_idx, time)
            plot_agents(ax, agents_now, agent_colors)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def plot_trajectory(agent_timeseries, config, out_dir='out', filename='trajectory'):
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    x_length = bounds[0]
    y_length = bounds[1]
    y_ratio = y_length / x_length

    # get agents
    times = np.array(agent_timeseries['time'])
    agents = agent_timeseries['agents']

    # get each agent's trajectory
    trajectories = {}
    for agent_id, data in agents.items():
        trajectories[agent_id] = []
        for time_data in data:
            x, y = time_data['location']
            theta = time_data['angle']
            pos = [x, y, theta]
            trajectories[agent_id].append(pos)

    # make the figure
    fig = plt.figure(figsize=(8, 8*y_ratio))

    for agent_id, agent_trajectory in trajectories.items():
        # convert trajectory to 2D array
        locations_array = np.array(agent_trajectory)
        x_coord = locations_array[:, 0]
        y_coord = locations_array[:, 1]

        # make multi-colored trajectory
        points = np.array([x_coord, y_coord]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=plt.get_cmap('cool'))
        lc.set_array(times)
        lc.set_linewidth(3)

        # plot line
        line = plt.gca().add_collection(lc)
        plt.plot(x_coord[0], y_coord[0], color=(0.0, 0.8, 0.0), marker='*')  # starting point
        plt.plot(x_coord[-1], y_coord[-1], color='r', marker='*')  # ending point

    plt.xlim((0, x_length))
    plt.ylim((0, y_length))

    # color bar
    cbar = plt.colorbar(line, ticks=[times[0], times[-1]])  # TODO --adjust this for full timeline
    cbar.set_label('time (s)', rotation=270)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def plot_motility(timeseries, out_dir='out', filename='motility_analysis'):
    expected_speed = 14.2  # um/s (Berg)
    expected_angle_between_runs = 68 # degrees (Berg)

    times = timeseries['time']
    agents = timeseries['agents']

    motility_analysis = {
        agent_id: {
            'speed': [],
            'angle': [],
            'thrust': [],
            'torque': []}
        for agent_id in agents.keys()}

    for agent_id, agent_data in agents.items():
        previous_location = [0,0]
        previous_time = times[0]

        # go through each time point for this agent
        for time, time_data in zip(times, agent_data):
            angle = time_data['angle']
            thrust, torque = time_data['motile_force']
            location = time_data['location']

            # get speed since last time
            if time != times[0]:
                dt = time - previous_time
                distance = (
                    (location[0] - previous_location[0]) ** 2 +
                    (location[1] - previous_location[1]) ** 2
                        ) ** 0.5
                speed = distance / dt  # um/sec
            else:
                speed = 0.0

            # save data
            motility_analysis[agent_id]['speed'].append(speed)
            motility_analysis[agent_id]['angle'].append(angle)
            motility_analysis[agent_id]['thrust'].append(thrust)
            motility_analysis[agent_id]['torque'].append(torque)

            # save previous location and time
            previous_location = location
            previous_time = time

    # plot results
    cols = 1
    rows = 3
    fig = plt.figure(figsize=(6 * cols, 1.5 * rows))
    plt.rcParams.update({'font.size': 12})

    ax1 = plt.subplot(rows, cols, 1)
    for agent_id, analysis in motility_analysis.items():
        speed = analysis['speed']
        avg_speed = np.mean(speed)
        ax1.plot(times, speed, label=agent_id)
        # ax1.axhline(y=avg_speed, color='b', linestyle='dashed', label='mean')

    ax1.axhline(y=expected_speed, color='r', linestyle='dashed', label='expected mean')
    ax1.set_ylabel(u'speed \n (\u03bcm/sec)')
    ax1.set_xlabel('time')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2 = plt.subplot(rows, cols, 2)
    for agent_id, analysis in motility_analysis.items():
        thrust = analysis['thrust']
        ax2.plot(times, thrust, label=agent_id)
    ax2.set_ylabel('thrust')

    ax3 = plt.subplot(rows, cols, 3)
    for agent_id, analysis in motility_analysis.items():
        torque = analysis['torque']
        ax3.plot(times, torque, label=agent_id)
    ax3.set_ylabel('torque')

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def init_axes(fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.4f} s'.format(float(time))
        plt.title(plot_title, y=1.08)
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'multibody')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='multibody')
    parser.add_argument('--mother', '-m', action='store_true', default=False)
    parser.add_argument('--motility', '-o', action='store_true', default=False)
    parser.add_argument('--growth', '-g', action='store_true', default=False)
    args = parser.parse_args()

    if args.mother:
        run_mother_machine()
    elif args.motility:
        run_motility(out_dir)
    elif args.growth:
        data = run_growth_division()
