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
from vivarium.compartment.emitter import timeseries_from_data
from vivarium.compartment.process import (
    Process,
    COMPARTMENT_STATE)
from vivarium.compartment.composition import (
    process_in_compartment,
    simulate_process)
from vivarium.processes.Vladimirov2008_motor import run, tumble

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

# colors for phylogeny initial agents
HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_HUE = HUES[0]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]

# agent port keys
AGENT_KEYS = ['location', 'angle', 'volume', 'length', 'width', 'mass', 'forces']
NON_AGENT_KEYS = ['fields', 'time', 'global', COMPARTMENT_STATE]


def get_volume(length, width):
    '''
    V = (4/3)*PI*r^3 + PI*r^2*a
    l = a + 2*r
    '''
    radius = width / 2
    cylinder_length = length - width
    volume = cylinder_length * (PI * radius ** 2) + (4 / 3) * PI * radius ** 3

    return volume

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

        # initialize agents
        self.agents = {}
        agents = initial_parameters.get('agents', {})
        for agent_id, specs in agents.items():
            self.add_body_from_center(agent_id, specs)

        # all initial agents get a key under a single port
        ports = {'agents': ['agents']}

        parameters = {}
        parameters.update(initial_parameters)

        super(Multibody, self).__init__(ports, parameters)

    def default_settings(self):
        agents = {agent_id: self.get_body_specs(agent_id)
                for agent_id in self.agents.keys()}
        state = {'agents': {'agents': agents}}

        schema = {'agents': {'agents': {'updater': 'merge'}}}

        default_emitter_keys = {
            port_id: keys for port_id, keys in self.ports.items()}

        return {
            'state': state,
            'schema': schema,
            'emitter_keys': default_emitter_keys,
        }

    def next_update(self, timestep, states):
        agents = states['agents']['agents']

        # check if an agent has been removed
        removed_agents = [
            agent_id for agent_id in self.agents.keys() if agent_id not in agents.keys()]
        for agent_id in removed_agents:
            del self.agents[agent_id]

        # update agents, add new agents
        for agent_id, specs in agents.items():
            if agent_id in self.agents:
                self.update_body(agent_id, specs)
            else:
                self.add_body_from_center(agent_id, specs)

        # run simulation
        self.run(timestep)

        # get new agent specs
        new_agents = {
            agent_id: self.get_body_specs(agent_id)
            for agent_id in self.agents.keys()}

        return {'agents': {'agents': new_agents}}

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

            motile_force = [10, 0.0]  # TODO -- remove this

            # add directly to angular velocity
            body.angular_velocity += torque
            # force-based torque
            # if torque != 0.0:
            #     motile_force = get_force_with_angle(thrust, torque)

        scaled_motile_force = [thrust * self.force_scaling for thrust in motile_force]

        # import ipdb;
        # ipdb.set_trace()
        # TODO -- get motile force through
        # TODO -- need to set lattice size?


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

        # add body to agents dictionary
        self.agents[body_id] = (body, shape)

    def update_body(self, body_id, specs):

        length = specs.get('length')
        width = specs.get('width')
        mass = specs.get('mass')
        motile_force = specs.get('motile_force', [0, 0])

        body, shape = self.agents[body_id]
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
        self.space.add(new_body, new_shape)
        self.space.remove(body, shape)

        # update body
        self.agents[body_id] = (new_body, new_shape)

    def get_body_specs(self, agent_id):
        body, shape = self.agents[agent_id]
        width, length = body.dimensions
        position = body.position

        # TODO -- velocity? angular velocity? forces?
        return {
            'location': [position[0], position[1]],
            'angle': body.angle,
            'length': length,
            'width': width,
            'mass': body.mass,
            'motile_force': [0, 0]}


# test functions
def get_n_dummy_agents(n_agents):
    return {agent_id: None for agent_id in range(n_agents)}

def random_body_config(n_agents=10, bounds=[10, 10]):
    agents = get_n_dummy_agents(n_agents)
    agent_config = {
        agent_id: {
            'location': [
                np.random.uniform(0, bounds[0]),
                np.random.uniform(0, bounds[1])],
            'angle': np.random.uniform(0, 2 * PI),
            'volume': 1,
            'length': 1.0,
            'width': 0.5,
            'mass': 1,
            'forces': [0, 0]}
        for agent_id in agents.keys()}

    return {
        'agents': agent_config,
        'bounds': bounds,
        'jitter_force': 1e1}

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
                ax = init_axes(fig, bounds[0], bounds[1], grid, row_idx, col_idx, time)
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

def init_axes(fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.4f} hr'.format(float(time) / 60. / 60.)
        plt.title(plot_title, y=1.08)
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax

def test_multibody(config=random_body_config(), time=1):
    multibody = Multibody(config)
    settings = {
        'total_time': time,
        'return_raw_data': True,
        'environment_port': 'external',
        'environment_volume': 1e-2}
    return simulate_process(multibody, settings)


motility_test_settings = {
    'timestep': 1,
    'total_time': 50}

def test_motility(config, settings=motility_test_settings):
    total_time = settings['total_time']
    timestep = settings['timestep']
    # n_agents = len(config['agents'])

    # make the process
    multibody = Multibody(config)
    compartment = process_in_compartment(multibody)
    agents_store = compartment.states['agents']
    agents_state = agents_store.state

    # initialize hidden agent motile states, and update agent motile_forces in agent store
    agent_motile_states = {}
    motile_forces = {}
    for agent_id, specs in agents_state['agents'].items():
        motor_force = run()
        agent_motile_states[agent_id] = {
            'motor_state': 1,  # 0 for run, 1 for tumble
            'time_in_motor_state': 0}
        motile_forces[agent_id] = {'motile_force': motor_force}
    compartment.send_updates({'agents': [{'agents': motile_forces}]})

    ## run simulation
    # test run/tumble
    time = 0
    while time < total_time:
        time += timestep

        # TODO -- update motile force and apply to state

        update = compartment.update(timestep)

    return compartment.emitter.get_data()


def plot_motility(timeseries, out_dir='out', filename='motility'):
    # time of motor behavior without chemotaxis
    run_time = 0.42  # s (Berg)
    tumble_time = 0.14  # s (Berg)

    expected_speed = 14.2  # um/s (Berg)
    expected_angle_between_runs = 68 # degrees (Berg)

    times = timeseries['time']
    agents = timeseries['agents']

    motility_analysis = {}
    for agent_id, agent_data in agents.items():
        previous_location = []

        # go through each time point
        for time, time_data in zip(times, agent_data):
            if agent_id not in motility_analysis:
                previous_time = time
                previous_location = time_data['location']
                motility_analysis[agent_id] = {
                    'speed': [0],
                }
            else:
                location = time_data['location']

                # get speed since last time
                dt = time - previous_time
                distance = (
                    (location[0] - previous_location[0]) ** 2 +
                    (location[1] - previous_location[1]) ** 2
                        ) ** 0.5
                speed = distance / dt  # um/sec
                motility_analysis[agent_id]['speed'].append(speed)

                # previous_location
                previous_location = location

    # plot results
    cols = 1
    rows = 1
    fig = plt.figure(figsize=(6 * cols, 1.5 * rows))
    plt.rcParams.update({'font.size': 12})

    ax1 = plt.subplot(rows, cols, 1)
    for agent_id, analysis in motility_analysis.items():
        speed = analysis['speed']
        avg_speed = np.mean(speed)
        ax1.plot(times, speed, label=agent_id)

    ax1.axhline(y=avg_speed, color='b', linestyle='dashed', label='mean')
    ax1.axhline(y=expected_speed, color='r', linestyle='dashed', label='expected mean')
    ax1.set_ylabel(u'speed \n (\u03bcm/sec)')
    ax1.set_xlabel('time')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))


    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)




if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'multibody')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # config = random_body_config(10)
    # data = test_multibody(config, 20)
    #
    # # make snapshot
    # agents = {time: time_data['agents']['agents'] for time, time_data in data.items()}
    # fields = {}
    # plot_snapshots(agents, fields, config, out_dir, 'snapshots')

    # test motility
    motility_config = {
        'jitter_force': 0,
        'bounds': [100, 100],
    }
    motility_config.update(random_body_config(1))
    motility_data = test_motility(motility_config)

    # make motility plot
    reduced_data = {time: data['agents'] for time, data in motility_data.items()}
    motility_timeseries = timeseries_from_data(reduced_data)
    plot_motility(motility_timeseries, out_dir)

    # make motility snapshot
    agents = {time: time_data['agents']['agents'] for time, time_data in motility_data.items()}
    fields = {}
    plot_snapshots(agents, fields, motility_config, out_dir, 'motility_snapshots')
