from __future__ import absolute_import, division, print_function

import copy
import time
import uuid

from vivarium.actor.control import ActorControl, AgentCommand
from vivarium.environment.make_media import Media
from vivarium.utils.units import units


class ShepherdControl(ActorControl):
    """Send messages to agents in the system to control execution."""

    def __init__(self, actor_config):
        super(ShepherdControl, self).__init__(str(uuid.uuid1()), actor_config)

    def add_cell(self, agent_type, actor_config):
        # TODO(jerry): Bring back the --variant choice?
        self.add_agent(
            str(uuid.uuid1()),
            agent_type,
            actor_config)

    def lattice_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id()
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        # make media
        media_id = args.get('media', 'minimal')
        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 {}, 7200 end'.format(media_id)

        lattice_config = copy.deepcopy(actor_config)
        lattice_config.update({
            'name': 'lattice_experiment',
            'timeline_str': timeline_str,
            'media_id': media_id,
            'emit_fields': ['GLC'],
            'boot': 'vivarium.environment.boot',
            'run_for': 4.0,
            'diffusion': 1000,
            'edge_length_x': 20,
            'patches_per_edge_x': 10,
            'translation_jitter': 0.1,
            'rotation_jitter': 0.01})

        self.add_agent(experiment_id, 'lattice', lattice_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        cell_config = dict(actor_config, **{
            'boot': args.get('agent_boot', 'vivarium.environment.boot'),
            'outer_id': experiment_id,
            'working_dir': args['working_dir']})

        for index in range(num_cells):
            self.add_cell(
                args['type'] or 'lookup',
                dict(cell_config, seed=index))

    def long_lattice_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id()
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        # make media
        media_id = args.get('media', 'minimal')
        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 {}, 3600 end'.format(media_id)

        emit_fields = ['GLC']

        lattice_config = {
            'name': 'long_lattice_experiment',
            'timeline_str': timeline_str,
            'media_id': media_id,
            'emit_fields': emit_fields,
            'boot': 'vivarium.environment.boot',
            'run_for': 4.0,
            'edge_length_x': 80.0,
            'edge_length_y': 20.0,
            'patches_per_edge_x': 8,
            'translation_jitter': 0.1,
            'rotation_jitter': 0.01}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'lookup', dict(actor_config, **{
                'boot': args.get('agent_boot', 'vivarium.environment.boot'),
                'boot_config': {},
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def large_lattice_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('large_lattice')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        media_id = args.get('media', 'minimal')
        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 {}, 7200 end'.format(media_id)

        emit_fields = ['GLC']

        lattice_config = {
            'name': 'large_lattice_experiment',
            'timeline_str': timeline_str,
            'media_id': media_id,
            'emit_fields': emit_fields,
            'run_for': 2.0,
            'edge_length_x': 50.0,
            'patches_per_edge_x': 10,
            'gradient': {
                'type': 'linear',
                'molecules': {
                    'GLC': {
                        'center': [0.0, 0.0],
                        'slope': -1.0 / 250.0},
                }},
        }

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'ecoli', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def small_lattice_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('small_lattice')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        media_id = args.get('media', 'minimal')
        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 {}, 3600 end'.format(media_id)

        emit_fields = ['GLC']

        lattice_config = {
            'name': 'small_lattice_experiment',
            'timeline_str': timeline_str,
            'run_for': 2.0,
            'emit_fields': emit_fields,
            'edge_length_x': 10.0,
            'patches_per_edge_x': 10}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'metabolism', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def glc_g6p_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('glc-g6p')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        media_id = 'GLC_G6P'
        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 {}, 3600 end'.format(media_id)

        emit_fields = ['GLC', 'G6P']

        lattice_config = {
            'name': 'glc_g6p_experiment',
            'timeline_str': timeline_str,
            'run_for': 2.0,
            'emit_fields': emit_fields}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        time.sleep(10)  # TODO(jerry): Wait for the Lattice to boot

        for index in range(num_cells):
            self.add_cell(args['type'] or 'metabolism', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def ecoli_core_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('glc-g6p')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        timeline_str = args.get('timeline')
        if not timeline_str:
            timeline_str = '0 ecoli_core_GLC 1.0 L + lac__D_e 2.0 mmol 0.1 L, 21600 end'

        lattice_config = {
            'name': 'ecoli_core_experiment',
            'timeline_str': timeline_str,
            'edge_length_x': 15.0,
            'patches_per_edge_x': 15,
            'run_for': 2.0,
            'diffusion': 1e-3,
            'depth': 1e-2,
            'translation_jitter': 1.0,
            'emit_fields': [
                'glc__D_e',
                'lac__D_e']}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        time.sleep(10)

        for index in range(num_cells):
            self.add_cell(args['type'] or 'kinetic_FBA', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def chemotaxis_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('chemotaxis')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        media_id = 'MeAsp'
        media = {'GLC': 20.0,
                 'MeAsp': 3600.0}
        new_media = {media_id: media}
        # timeline_str = '0 {}, 14400 end'.format(media_id)
        timeline_str = '0 {}, 1800 end'.format(media_id)

        lattice_config = {
            'name': 'chemotaxis_experiment',
            'description': 'a long environment with a static gradient of glucose and a-methyl-DL-aspartic acid (MeAsp) '
                           'for observing chemotactic cells in action. Optimal chemotaxis is observed in a narrow range '
                           'of CheA activity, where concentration of CheY-P falls into the operating range of flagellar motors.',
            'timeline_str': timeline_str,
            'new_media': new_media,


            'emit_fields': ['GLC','MeAsp'],
            'run_for': 0.1,  # high coupling between cell and env requires short exchange timestep
            'static_concentrations': True,
            'cell_placement': [0.05, 0.5],  # place cells at bottom of gradient
            'gradient': {
                'type': 'linear',
                'molecules': {
                    'GLC': {
                        'center': [1.0, 0.5],
                        'slope': -1e-2},
                    'MeAsp': {
                        'center': [1.0, 0.5],
                        'slope': -2e0}
                }},
            'edge_length_x': 2000.0,
            'edge_length_y': 400.0,
            'translation_jitter': 1.0,
            # 'rotation_jitter': 0.05,
            'patches_per_edge_x': 100}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        # give lattice time before adding the cells
        time.sleep(15)

        for index in range(num_cells):
            self.add_cell(args['type'] or 'chemotaxis', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'seed': index}))


    def swarm_experiment(self, args, actor_config):
        experiment_id = args['experiment_id']
        if not experiment_id:
            experiment_id = self.get_experiment_id('chemotaxis')
        num_cells = args['number']
        print('Creating lattice agent_id {} and {} cell agents\n'.format(
            experiment_id, num_cells))

        ## Make media: GLC_G6P with MeAsp
        # get GLC_G6P media
        make_media = Media()
        media_id = 'GLC_G6P'
        media1 = make_media.get_saved_media('GLC_G6P', True)

        # make MeAsp media
        ingredients = {
            'MeAsp': {
                'counts': 1.0 * units.mmol,
                'volume': 0.001 * units.L}}
        media2 = make_media.make_recipe(ingredients, True)

        # combine the medias
        media = make_media.combine_media(media1, 0.999 * units.L, media2, 0.001 * units.L)
        media_id = 'GLC_G6P_MeAsp'

        # make timeline with new media
        new_media = {media_id: media}
        timeline_str = '0 {}, 3600 end'.format(media_id)

        lattice_config = {
            'name': 'swarm_experiment',
            'description': 'a large experiment for running swarms of chemotactic cells',
            'cell_placement': [0.5, 0.5], # place cells at center of lattice
            'timeline_str': timeline_str,
            'new_media': new_media,
            'run_for' : 2.0,
            'emit_fields': ['MeAsp', 'GLC'],
            'static_concentrations': False,
            'diffusion': 0.001,
            'edge_length_x': 100.0,
            'edge_length_y': 100.0,
            'patches_per_edge_x': 50}

        actor_config['boot_config'].update(lattice_config)
        self.add_agent(experiment_id, 'lattice', actor_config)

        # give lattice time before adding the cells
        time.sleep(15)

        for index in range(num_cells):
            self.add_cell(args['type'] or 'chemotaxis', dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'seed': index}))

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
            'long-experiment',
            'large-experiment',
            'small-experiment',
            'chemotaxis-experiment',
            'swarm-experiment',
            'glc-g6p-experiment',
            'ecoli-core-experiment'] + choices

        super(EnvironmentCommand, self).__init__(
            full_choices,
            {}, # no additional args
            full_description)

    def experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.lattice_experiment(args, self.actor_config)
        control.shutdown()

    def long_experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.long_lattice_experiment(args, self.actor_config)
        control.shutdown()

    def large_experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.large_lattice_experiment(args, self.actor_config)
        control.shutdown()

    def small_experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.small_lattice_experiment(args, self.actor_config)
        control.shutdown()

    def glc_g6p_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.glc_g6p_experiment(args, self.actor_config)
        control.shutdown()

    def ecoli_core_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config})
        control.ecoli_core_experiment(args)
        control.shutdown()

    def chemotaxis_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.chemotaxis_experiment(args, self.actor_config)
        control.shutdown()

    def swarm_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.kafka_config()})
        control.swarm_experiment(args, self.actor_config)
        control.shutdown()

    def add_arguments(self, parser):
        parser = super(EnvironmentCommand, self).add_arguments(parser)

        parser.add_argument(
            '-e', '--experiment_id',
            type=str,
            help='The experiment id')

        parser.add_argument(
            '-f', '--emit_field',
            type=str,
            default='minimal',
            help='emitted media field')

        parser.add_argument(
            '-m', '--media',
            type=str,
            default='minimal',
            help='The environment media')

        parser.add_argument(
            '-t', '--timeline',
            type=str,
            default=None,
            help='The timeline')

        return parser

def run():
    command = EnvironmentCommand()
    command.execute()

if __name__ == '__main__':
    run()
