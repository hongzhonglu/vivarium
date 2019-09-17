from __future__ import absolute_import, division, print_function

import os
import shutil

from lens.actor.outer import Outer
from lens.actor.boot import BootAgent

from lens.environment.lattice import EnvironmentSpatialLattice
from lens.environment.make_media import Media


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
    working_dir = agent_config.get('working_dir', os.getcwd())
    media_id = agent_config.get('media_id', 'minimal')
    media = agent_config.get('media', {})
    print("Media condition: {}".format(media_id))
    if not media:
        make_media = Media()
        media = make_media.get_saved_media(media_id)

    output_dir = os.path.join(working_dir, 'out', agent_id)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)

    boot_config = {
        'output_dir': output_dir,
        'concentrations': media,
    }
    boot_config.update(agent_config)
    environment = EnvironmentSpatialLattice(boot_config)

    return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

def boot_glc_g6p(agent_id, agent_type, agent_config):
    working_dir = agent_config.get('working_dir', os.getcwd())

    media_id = 'GLC_G6P'
    make_media = Media()
    media = make_media.get_saved_media(media_id)

    print("Media condition: {}".format(media_id))
    output_dir = os.path.join(working_dir, 'out', agent_id)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)

    boot_config = {
        'output_dir': output_dir,
        'concentrations': media,
        'run_for': 10.0,
        'depth': 0.0001,  # 3000 um is default
        'edge_length': 10.0,
        'patches_per_edge': 1,
    }

    boot_config.update(agent_config)
    environment = EnvironmentSpatialLattice(boot_config)

    return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

def boot_glc_lct(agent_id, agent_type, agent_config):
    working_dir = agent_config.get('working_dir', os.getcwd())

    media_id = 'GLC_LCT'
    make_media = Media()
    media = make_media.get_saved_media(media_id)

    print("Media condition: {}".format(media_id))
    output_dir = os.path.join(working_dir, 'out', agent_id)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)

    boot_config = {
        'output_dir': output_dir,
        'concentrations': media,
    }
    boot_config.update(agent_config)
    environment = EnvironmentSpatialLattice(boot_config)

    return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

class BootEnvironment(BootAgent):
    def __init__(self):
        super(BootEnvironment, self).__init__()
        self.agent_types = {
            'lattice': boot_lattice,
            'sugar1': boot_glc_g6p,
            'sugar2': boot_glc_lct,
            }

def run():
    boot = BootEnvironment()
    boot.execute()

if __name__ == '__main__':
    run()
