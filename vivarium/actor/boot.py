from __future__ import absolute_import, division, print_function

import copy
import json
import argparse

from vivarium.actor.control import DEFAULT_KAFKA_CONFIG, DEFAULT_EMITTER_CONFIG, DEFAULT_DESTINATIONS, distribute_arguments
from vivarium.actor.outer import Outer
from vivarium.actor.inner import Inner
from vivarium.actor.stub import SimulationStub, EnvironmentStub
from vivarium.utils.dict_utils import deep_merge

def boot_outer(agent_id, agent_type, actor_config):
    """
    Initialize the `EnvironmentStub`, pass it to the `Outer` agent and launch the process.

    This is a demonstration of how to initialize an Outer component. In the place of
    `EnvironmentStub` you would substitute your own environment class that meets the interface
    defined in `Outer`.
    """

    volume = 1
    concentrations = {
        'yellow': 5,
        'green': 11,
        'red': 44,
        'blue': 12}

    environment = EnvironmentStub(volume, concentrations)
    print('outer! {}'.format(actor_config))
    return Outer(agent_id, agent_type, actor_config, environment)

def boot_inner(agent_id, agent_type, actor_config):
    """
    Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

    This is a demonstration of how to initialize an Inner component. When creating your
    own simulation you would supply a class that meets the same interface as the `SimulationStub`
    that would be driven by the Inner agent in response to messages from its corresponding
    Outer agent.
    """
    if 'outer_id' not in actor_config:
        raise ValueError("--outer-id required")

    agent_id = agent_id
    outer_id = actor_config['outer_id']
    simulation = SimulationStub()

    print('inner... {}'.format(actor_config))
    return Inner(
        agent_id,
        outer_id,
        agent_type,
        actor_config,
        simulation)

class BootAgent(object):
    """
    Boot agents from the command line.
    """

    def __init__(self):
        self.description='Boot agents for the environmental context simulation'
        self.agent_types = {
            'outer': boot_outer,
            'inner': boot_inner}

    def add_arguments(self, parser):
        parser.add_argument(
            '-i', '--id',
            required=True,
            help='id of the new agent')

        parser.add_argument(
            '-o', '--outer-id',
            help="ID of the new agent's outer environment agent")

        parser.add_argument(
            '-t', '--type',
            required=True,
            choices=self.agent_types,
            help='type of the new agent')

        parser.add_argument(
            '--config',
            default='{}',
            help='''JSON configuration dictionary for the new agent.''')

        parser.add_argument(
            '--kafka-host',
            default=None,
            help='address for Kafka server')

        parser.add_argument(
            '--mongo-host',
            default=None,
            help='address for Mongo server')

        parser.add_argument(
            '--emitter',
            default='database',
            help='which emitter to use')

        return parser

    def execute(self):
        parser = argparse.ArgumentParser(description=self.description)
        parser = self.add_arguments(parser)
        parse_args = parser.parse_args()

        args = vars(parse_args)
        actor_config = dict(json.loads(args['config']))

        kafka_config = copy.deepcopy(DEFAULT_KAFKA_CONFIG)
        emitter_config = copy.deepcopy(DEFAULT_EMITTER_CONFIG[args['emitter']])
        default_config = {
            'kafka_config': kafka_config,
            'boot_config': {
                'emitter': emitter_config}}

        actor_config = deep_merge(default_config, actor_config)
        actor_config = distribute_arguments(
            DEFAULT_DESTINATIONS,
            args,
            actor_config)

        if args['outer_id']:
            actor_config.setdefault('outer_id', args['outer_id'])

        boot = self.agent_types[args['type']]
        agent = boot(args['id'], args['type'], actor_config)
        agent.start()

def run():
    boot = BootAgent()
    boot.execute()

if __name__ == '__main__':
    run()
