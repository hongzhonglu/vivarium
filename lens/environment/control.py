from __future__ import absolute_import, division, print_function

import time
import uuid

from lens.actor.control import ActorControl, AgentCommand


class ShepherdControl(ActorControl):
    """Send messages to agents in the system to control execution."""

    def __init__(self, agent_config):
        super(ShepherdControl, self).__init__(str(uuid.uuid1()), agent_config)

    def add_cell(self, agent_type, agent_config):
        # TODO(jerry): Bring back the --variant choice?
        self.add_agent(
            str(uuid.uuid1()),
            agent_type,
            agent_config)

    def lattice_experiment(self, args):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id()
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        # make media
        timeline_str = args.get('timeline','')
        media_id = args.get('media', 'minimal')
        lattice_config = {
            'timeline_str': timeline_str,
            'media_id': media_id,
            'boot': 'lens.environment.boot',
            'run_for': 4.0,
            'media_id': media_id,
            'translation_jitter': 0.5,
            'rotation_jitter': 0.005}

        self.add_agent(experiment_id, 'lattice', lattice_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'lookup', {
                'boot': args.get('agent_boot', 'lens.environment.boot'),
                'boot_config': {},
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index})

    def large_lattice_experiment(self, args):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('large_lattice')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        timeline_str = args.get('timeline','')
        media_id = args.get('media')
        lattice_config = {
            'timeline_str': timeline_str,
            'media_id': media_id,
            'run_for': 2.0,
            'edge_length': 50.0,
            'patches_per_edge': 10,
        }

        self.add_agent(experiment_id, 'lattice', lattice_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'ecoli', {
                'boot': 'lens.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index})

    def glc_g6p_experiment(self, args):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('glc-g6p')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        experiment_config = {
            'run_for': 2.0,
            'media_id': 'GLC_G6P'}

        self.add_agent(experiment_id, 'lattice', experiment_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'metabolism', {
                'boot': 'lens.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index})

    def chemotaxis_experiment(self, args):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('chemotaxis')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        media_id = 'MeAsp'
        media = {'GLC': 20.0,
                 'MeAsp': 1.0}
        new_media = {media_id: media}
        timeline_str = '0 {}, 3600 end'.format(media_id)

        chemotaxis_config = {
            'timeline_str': timeline_str,
            'new_media': new_media,
            'run_for' : 2.0,
            'static_concentrations': True,
            'gradient': {
                'seed': True,
                'molecules': {
                    'GLC':{
                        'center': [0.5, 0.5],
                        'deviation': 30.0},
                    'MeAsp': {
                        'center': [0.5, 0.5],
                        'deviation': 30.0}
                }},
            'diffusion': 0.0,
            'translation_jitter': 0.5,
            'rotation_jitter': 0.005,
            'edge_length': 200.0,
            'patches_per_edge': 50}
        self.add_agent(experiment_id, 'lattice', chemotaxis_config)

        # give lattice time before adding the cells
        time.sleep(15)

        for index in range(num_cells):
            self.add_cell(args['type'] or 'chemotaxis', {
                'boot': 'lens.environment.boot',
                'outer_id': experiment_id,
                'seed': index})


class EnvironmentCommand(AgentCommand):
    """
    Extend `AgentCommand` with new commands related to the lattice and ecoli experiments
    """

    def __init__(self, choices=[], description=''):
        full_description = '''
    Run an agent for the environmental context simulation.
    
    The commands are:
    `add --id OUTER_ID [--type T] [--config C]` ask the Shepherd to add an agent of
        type T with JSON configuration C to the environment OUTER_ID,
    `remove --id AGENT_ID` ask all Shepherds to remove agent AGENT_ID,
    `remove --prefix ID` ask all Shepherds to remove agents "ID*",
    `run [--id OUTER_ID]` start or resume one or all simulations,
    `pause [--id OUTER_ID]` pause one or all simulations,
    `divide --id AGENT_ID` ask a cell agent to divide,
    `shutdown [--id OUTER_ID]` shut down one or all environment agents and their
         connected agents,
    `experiment [--number N] [--type T] [--working-dir D]` ask the Shepherd to run
        a lattice environment with N agents of type T,
    'glc-g6p-experiment [--number N] [--type T]` ask the Shepherd to run a
        chemotaxis environment with N agents of type T
    'chemotaxis-experiment [--number N] [--type T]` ask the Shepherd to run a
        chemotaxis environment with N agents of type T
    ''' + description

        full_choices = [
			'large-experiment',
            'chemotaxis-experiment',
            'glc-g6p-experiment'] + choices

        super(EnvironmentCommand, self).__init__(
			full_choices,
			full_description)

    def experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config})
        control.lattice_experiment(args)
        control.shutdown()

    def large_experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config})
        control.large_lattice_experiment(args)
        control.shutdown()

    def glc_g6p_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config})
        control.glc_g6p_experiment(args)
        control.shutdown()

    def chemotaxis_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config})
        control.chemotaxis_experiment(args)
        control.shutdown()

    def add_arguments(self, parser):
        parser = super(EnvironmentCommand, self).add_arguments(parser)

        parser.add_argument(
            '-e', '--experiment_id',
            type=str,
            help='The experiment id')

        parser.add_argument(
            '-m', '--media',
            type=str,
            default='minimal',
            help='The environment media')

        parser.add_argument(
            '-t', '--timeline',
            type=str,
            # default='0 minimal',
            help='The timeline')

        return parser

def run():
    command = EnvironmentCommand()
    command.execute()

if __name__ == '__main__':
    run()
