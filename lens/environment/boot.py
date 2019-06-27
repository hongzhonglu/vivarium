from __future__ import absolute_import, division, print_function

from lens.agent import Outer
from lens.agent import Inner
from lens.agent import BootAgent

from environment.lattice import EnvironmentSpatialLattice
from environment.surrogates.metabolism import Metabolism
from environment.surrogates.chemotaxis_minimal import ChemotaxisMinimal
from environment.surrogates.chemotaxis_MWC_sensors import Chemotaxis
from environment.surrogates.endocrine import Endocrine
from environment.surrogates.transport_kinetics import TransportKinetics
from environment.surrogates.transport_lookup import TransportLookup
from environment.surrogates.division import Division
from environment.condition.make_media import Media


DEFAULT_COLOR = [0.6, 0.4, 0.3]

class EnvironmentAgent(Outer):
    def build_state(self):
        lattice = {
            molecule: self.environment.lattice[index].tolist()
            for index, molecule in enumerate(self.environment.get_molecule_ids())}

        simulations = {
            agent_id: {
                'volume': simulation['state']['volume'],
                'width': self.environment.simulations[agent_id]['state'].get('width', self.environment.cell_radius*2), # TODO (Eran) initialize length, width in cellSimulation
                'length': self.environment.simulations[agent_id]['state'].get('length', self.environment.cell_radius*4),
                'color': simulation['state'].get('color', DEFAULT_COLOR),
                'location': self.environment.locations[agent_id][0:2].tolist(),
                'corner_location': self.environment.corner_locations[agent_id][0:2].tolist(),
                'orientation': self.environment.locations[agent_id][2],
                'parent_id': simulation.get('parent_id', '')}
            for agent_id, simulation in self.environment.simulations.iteritems()}

        return {
            'outer_id': self.agent_id,
            'agent_type': 'ecoli',
            'running': not self.paused,
            'time': self.environment.time(),
            'edge_length': self.environment.edge_length,
            'cell_radius': self.environment.cell_radius,   # TODO (Eran) -- remove this from environment config, should be an attribute of cellSimulations
            'lattice': lattice,
            'simulations': simulations}

    def update_state(self):
        self.send(
            self.topics['visualization_receive'],
            self.build_state(),
            print_send=False)

def boot_lattice(agent_id, agent_type, agent_config):
    media_id = agent_config.get('media_id', 'minimal')
    media = agent_config.get('media', {})
    print("Media condition: {}".format(media_id))
    if not media:
        make_media = Media()
        media = make_media.make_recipe(media_id)

    # # TODO (Eran) -- easier passing of medias to lattice. Put them in recipes?
    # media_id = 'toy'
    # media = {'A': 20.0,
    # 		 'E': 0.1,
    # 		 'D': 0.1,
    # 		 'F': 0.1,
    # 		 'H': 0.1,
    # 		 'O2': 0.1,
    # 		 }

    agent_config['concentrations'] = media
    environment = EnvironmentSpatialLattice(agent_config)

    return EnvironmentAgent(agent_id, agent_type, agent_config, environment)


# Metabolism surrogate initialize and boot
def initialize_metabolism(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return Metabolism(boot_config)

def boot_metabolism(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_metabolism)

    return inner


# minimal Chemotaxis surrogate initialize and boot
def initialize_chemotaxis_minimal(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return ChemotaxisMinimal(boot_config)

def boot_chemotaxis_minimal(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_chemotaxis_minimal)

    return inner


# Chemotaxis surrogate initialize and boot
def initialize_chemotaxis(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return Chemotaxis(boot_config)

def boot_chemotaxis(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_chemotaxis_minimal)

    return inner


# Endocrine surrogate initialize and boot
def initialize_endocrine(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return Endocrine(boot_config)

def boot_endocrine(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_endocrine)

    return inner


# Transport lookup minimal surrogate initialize and boot
def initialize_lookup_transport(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): essential options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization

    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return TransportLookup(boot_config)

def boot_lookup_transport(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'lookup': 'average',
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_lookup_transport)

    return inner


# Transport kinetic surrogate initialize and boot
def initialize_kinetic_transport(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): essential options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return TransportKinetics(boot_config)

def boot_kinetic_transport(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    state = {
        'lookup': 'average',
        'volume': 1.0,
        'environment_change': {}}
    agent_config['state'] = state
    options = {}

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        options,
        initialize_kinetic_transport)

    return inner


# division surrogate initialize and boot
def initialize_division(boot_config, synchronize_config):
    '''
    Args:
        boot_config (dict): options for initializing a simulation
        synchronize_config (dict): additional options that can be passed in for initialization
    Returns:
        simulation (CellSimulation): The actual simulation which will perform the calculations.
    '''
    boot_config.update(synchronize_config)
    return Division(boot_config)

def boot_division(agent_id, agent_type, agent_config):
    agent_id = agent_id
    outer_id = agent_config['outer_id']

    # initialize state and options
    volume = agent_config.get('volume', 1.0)

    # state of the cell that gets sent to the environment
    agent_config['state'] = {
        'volume': volume,
        'environment_change': {}}

    # boot_config previously called options, state sent to the cell upon boot
    boot_config = {
        'volume': volume
        }

    inner = Inner(
        agent_id,
        outer_id,
        agent_type,
        agent_config,
        boot_config,
        initialize_division)

    return inner


class BootEnvironment(BootAgent):
    def __init__(self):
        super(BootEnvironment, self).__init__()
        self.agent_types = {
            'lattice': boot_lattice,
            'chemotaxis-minimal': boot_chemotaxis_minimal,
            'chemotaxis': boot_chemotaxis,
            'metabolism': boot_metabolism,
            'endocrine': boot_endocrine,
            'kinetics': boot_kinetic_transport,
            'lookup': boot_lookup_transport,
            'division': boot_division,
            }

if __name__ == '__main__':
    boot = BootEnvironment()
    boot.execute()
