from __future__ import absolute_import, division, print_function

import uuid
import time
import copy

import vivarium.actor.event as event
from vivarium.actor.actor import Actor


class Simulation(object):
    """Interface for the Inner agent's cell simulation."""

    def time(self):
        """Return the current time according to this Simulation."""

    def apply_outer_update(self, update):
        """Apply the update received from the environment to this simulation."""

    def run_incremental(self, run_until):
        """Run this Simulation until the given time."""

    def generate_inner_update(self):
        """
        Generate the update that will be sent to the environment based on changes calculated
        by the Simulation during `run_incremental(run_until)`.

        If the dictionary returned by this function contains a `division` key it will trigger
        preparations for cell division in the environment. The value for this key is a pair of
        dictionaries, of which the only required key is `id` containing the new daughter id.
        Additional keys are specific to the particular implementation of `EnvironmentSimulation`
        that receives this update (see `outer.py`).
        """

    def divide(self):
        """Perform cell division on the simulation and return information about the daughters."""
        return []

    def finalize(self):
        """Release any resources and perform any final cleanup."""


class Inner(Actor):

    """
    Inner: an independent cell simulation in a larger environmental context.

    This class wraps an instance of Simulation into an Agent and mediates
    the message passing communication with the coordinating Outer agent running
    an environmental simulation.
    """

    def __init__(self, agent_id, agent_type, actor_config, sim_initialize):
        """
        Construct the agent.

        Args:
            agent_id (str): Unique identifier for this agent.
              This agent will only respond to messages addressed to its inner agent_id.
            agent_type (str): The type of this agent, for coordination with the agent shepherd.
            actor_config (dict): A dictionary containing any information needed to run this
              inner agent. The only required keys are `outer_id`, `boot` and `kafka_config`
                * `outer_id`: the id of the outer this inner is embedded in.
                * `declare`: dict of values to pass to outer in CELL_DECLARE message.
                * `boot_config`: dict of options to pass in to the simulation during initialization.
                * `kafka_config`: anything this agent needs to know about kafka, including:
                    * `host`: the Kafka server host address.
                    * `topics`: a dictionary mapping topic roles to specific topics used by the agent
                        to communicate with other agents. The relevant ones to this agent are:
                        * `cell_receive`: The topic this agent will receive messages on from the
                            environment or relevant control processes.
                        * `environment_receive`: The topic this agent will send messages to its
                            associated outer agent (given by `outer_id`) and environmental simulation
                        * `shepherd_receive`: The topic this agent will send messages on for
                            adding agents to and removing agents from the environment.
            sim_initialize: the function for initializing a simulation. actor_config['boot_config']
              and the response from the environment accompanying the CELL_SYNCHRONIZE message will
              be passed in.
        """

        self.outer_id = actor_config['outer_id']
        self.declare = actor_config.get('declare', {})
        self.boot_config = actor_config.get('boot_config', {})
        self.generation = actor_config.get('generation', 0)

        # mutating in place
        kafka_config = actor_config['kafka_config']
        kafka_config['subscribe'].append(
            kafka_config['topics']['cell_receive'])

        self.sim_initialize = sim_initialize

        super(Inner, self).__init__(agent_id, agent_type, actor_config)

    def preinitialize(self):
        time.sleep(1.0)

        self.send(self.kafka_config['topics']['environment_receive'], {
            'event': event.CELL_DECLARE,
            'agent_id': self.outer_id,
            'inner_id': self.agent_id,
            'agent_config': self.actor_config,
            'state': self.declare})

    def initialize_simulation(self, message):

        self.boot_config.update(message['state'])
        self.boot_config.update({
            'simulation_id': message['inner_id'],
            'experiment_id': message['outer_id']})
        self.simulation = self.sim_initialize(self.boot_config)
        self.send_initialize()

    def send_initialize(self):
        """
        Initialization: Register this inner agent with the outer agent.
        """
        now = self.simulation.time()
        state = self.simulation.generate_inner_update()

        self.send(self.topics['environment_receive'], {
            'time': now,
            'event': event.CELL_INITIALIZE,
            'outer_id': self.outer_id,
            'inner_id': self.agent_id,
            'agent_config': self.actor_config,
            'state': state})

    def cell_exchange(self, message):
        """
        Notify the environment simulation about what changes this cell simulation has produced.
        Also, handle the case of the cell dividing during this update by adding the information
        to its environment update and calling `divide_cell`.
        """

        stop = self.simulation.time()
        update = self.simulation.generate_inner_update()
        division = update.get('division', [])

        for index, daughter in enumerate(division):
            print('daughter {}: {}'.format(index+1, daughter))

        for daughter in division:
            assert daughter.get('id')

        self.send(self.topics['environment_receive'], {
            'event': event.CELL_EXCHANGE,
            'time': stop,
            'outer_id': self.outer_id,
            'inner_id': self.agent_id,
            'message_id': message['message_id'],
            'state': update})

        if division:
            self.divide_cell({}, division)

    def environment_update(self, message):
        """
        Apply the update from the environment to the cell simulation and respond with changes.
        """

        self.simulation.apply_outer_update(message['state'])
        self.simulation.run_incremental(message['run_until'])
        self.cell_exchange(message)

    def divide_cell(self, message, division):
        """
        Perform agent cell division.

        The generation count is increased and added to the daughter cells' actor_config.

        This sends three messages to the agent shepherd: one `ADD_AGENT` for each new daughter cell,
        and finally a `REMOVE_AGENT` for itself. These new agents will initialize and notify the
        outer agent, inheriting properties from their parent cell.
        """

        generation = self.generation + 1
        for index, daughter in enumerate(division):
            daughter_id = daughter['id']
            daughter_type = daughter.get(
                'type',
                message.get(
                    'daughter_type',
                    self.agent_type))

            daughter_config = copy.deepcopy(self.actor_config)
            daughter_config.update(daughter)

            actor_config = dict(
                daughter_config,
                index=index,
                parent_id=self.agent_id,
                outer_id=self.outer_id,
                generation=generation)

            # Send the inherited state data as a blob instead of a file path.
            inherited_state_path = actor_config.pop('inherited_state_path', None)
            add_agent_message = {
                'event': event.ADD_AGENT,
                'agent_id': daughter_id,
                'agent_type': daughter_type,
                'agent_config': actor_config}

            if inherited_state_path:
                with open(inherited_state_path, 'rb') as f:
                    add_agent_message['blobs'] = [f.read()]
                # TODO(jerry): Delete the file?

            self.send(self.topics['shepherd_receive'], add_agent_message)

        self.send(self.topics['shepherd_receive'], {
            'event': event.REMOVE_AGENT,
            'agent_id': self.agent_id})

    def cell_shutdown(self, message):
        self.send(self.topics['environment_receive'], {
            'event': event.CELL_SHUTDOWN,
            'outer_id': self.outer_id,
            'inner_id': self.agent_id})

        self.shutdown()

    def finalize(self):
        """Do any cleanup the simulation needs to perform before exiting."""

        self.simulation.finalize()

    def receive(self, topic, message):
        """
        Respond to messages from the environment.

        The inner agent responds to four messages:

        * ENVIRONMENT_SYNCHRONIZE: Receive any relevant information from the environment before
            the main cycle of mutual updates begins.
        * ENVIRONMENT_UPDATE: Receive the latest state from the environment simulation. The
            relevant keys in this update are:

            * `state`: a dictionary containing the updated state from the environment.
            * `run_until`: how long to run the cell simulation until before reporting back the new
                environmental changes.
            * `message_id`: the id of the message as provided by the outer agent,
                to be returned as an acknowledgement that the message was processed along with
                the updated environmental changes.

        * DIVIDE_CELL: Perform cell division immediately, whether the simulation is ready or not.
        * SHUTDOWN_AGENT: Shutdown this agent and terminate the process.
        """

        if message.get('inner_id', message.get('agent_id')) == self.agent_id:
            self.print_message(topic, message)
            message_event = message['event']

            if message_event == event.ENVIRONMENT_SYNCHRONIZE:
                self.initialize_simulation(message)

            elif message_event == event.ENVIRONMENT_UPDATE:
                self.environment_update(message)

            elif message_event == event.DIVIDE_CELL:
                self.simulation.divide()

            elif message_event == event.PAUSE_AGENT:
                pass

            elif message_event == event.SHUTDOWN_AGENT:
                self.cell_shutdown(message)

            else:
                print('unexpected event {}: {}'.format(message_event, message))
