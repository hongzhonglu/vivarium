'''
boot functions for initialize environment simulations:
Can be called from the command line
by specifying an environment "agent_type" and an outer-id:

>> python -m vivarium.environment.boot --type lattice --id lattice

boot functions for compartment simulations:
Can be called from the command line
after booting an environmental "outer" simulation and passing in the outer-id, the agent_type:

>> python -m vivarium.environment.boot --outer-id 'outer-id' --type 'agent_type'

Once an environment and agent are booted, they can be triggered with:

>> python -m vivarium.environment.control run --id 'outer-id'

'''

from __future__ import absolute_import, division, print_function

import os
import shutil

from vivarium.actor.inner import Inner
from vivarium.actor.outer import Outer
from vivarium.actor.boot import BootAgent
from vivarium.actor.emitter import get_emitter, configure_emitter

# environment
from vivarium.environment.lattice import EnvironmentSpatialLattice
from vivarium.environment.lattice_compartment import LatticeCompartment, generate_lattice_compartment
from vivarium.environment.make_media import Media
from vivarium.utils.units import units

# processes
from vivarium.processes.transport_lookup import TransportLookup
from vivarium.processes.BiGG_metabolism import BiGGMetabolism
from vivarium.processes.CovertPalsson2002_regulation import Regulation
from vivarium.processes.Kremling2007_transport import Transport
from vivarium.processes.growth import Growth
from vivarium.processes.protein_expression import ProteinExpression
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity
from vivarium.processes.membrane_potential import MembranePotential

# composites
from vivarium.composites.growth_division import compose_growth_division
from vivarium.composites.simple_chemotaxis import compose_simple_chemotaxis
from vivarium.composites.PMF_chemotaxis import compose_pmf_chemotaxis
from vivarium.composites.iFBA import compose_iFBA


DEFAULT_COLOR = [0.6, 0.4, 0.3]

def wrap_boot(initialize, initial_state):
    def boot(agent_id, agent_type, agent_config):
        initial_state.update(agent_config.get('declare', {}))
        agent_config['declare'] = initial_state  # 'declare' is for the environment

        return Inner(
            agent_id,
            agent_type,
            agent_config,
            initialize)

    return boot

def wrap_init_basic(make_process):
    def initialize(boot_config):
        boot_config.update({
            'emitter': {
                'type': 'database',
                'url': 'localhost:27017',
                'database': 'simulations',
                },
            })
        process = make_process(boot_config)  # 'boot_config', set in environment.control is the process' initial_parameters
        return generate_lattice_compartment(process, boot_config)

    return initialize

def wrap_init_composite(make_composite):
    def initialize(boot_config):
        # configure the emitter
        boot_config.update({
            'emitter': {
                'type': 'database',
                'url': 'localhost:27017',
                'database': 'simulations'},
            })

        # set up the the composite
        composite_config = make_composite(boot_config)
        processes = composite_config['processes']
        states = composite_config['states']
        options = composite_config['options']
        topology = options['topology']

        # update options
        emitter = configure_emitter(boot_config, processes, topology)
        options.update({'emitter': emitter})

        # create the compartment
        return LatticeCompartment(processes, states, options)

    return initialize

def wrap_boot_environment(intialize):
    def boot(agent_id, agent_type, agent_config):
        # get boot_config from initialize
        boot_config = intialize(agent_config)

        # paths
        working_dir = agent_config.get('working_dir', os.getcwd())
        output_dir = os.path.join(working_dir, 'out', agent_id)
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

        # emitter
        emitter_config = {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
            'experiment_id': agent_id}
        emitter = get_emitter(emitter_config)
        boot_config.update({
            'emitter': emitter,
            'output_dir': output_dir})

        # create the environment
        environment = EnvironmentSpatialLattice(boot_config)

        return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

    return boot

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
            for agent_id, simulation in self.environment.simulations.items()}

        return {
            'outer_id': self.agent_id,
            'agent_type': 'ecoli',
            'running': not self.paused,
            'time': self.environment.time(),
            'edge_length_x': self.environment.edge_length_x,
            'edge_length_y': self.environment.edge_length_y,
            'cell_radius': self.environment.cell_radius,   # TODO (Eran) -- remove this from environment config, should be an attribute of cellSimulations
            'lattice': lattice,
            'simulations': simulations}

    def update_state(self):
        self.send(
            self.topics['visualization_receive'],
            self.build_state(),
            print_send=False)

# Define environment initialization functions
def initialize_lattice(agent_config):
    # set up media
    media_id = agent_config.get('media_id', 'minimal')
    media = agent_config.get('media', {})
    if media:
        boot_config = {'concentrations': media}
    else:
        boot_config = {'media_id': media_id}
    boot_config.update(agent_config)

    return boot_config

def initialize_glc_g6p_small(agent_config):
    # set up media
    media_id = 'GLC_G6P'
    timeline_str = '0 {}, 1800 end'.format(media_id)  # (2hr*60*60 = 7200 s), (7hr*60*60 = 25200 s)
    boot_config = {
        'timeline_str': timeline_str,
        'run_for': 2.0,
        'depth': 1e-01, # 3000 um is default
        'edge_length_x': 1.0,
        'patches_per_edge_x': 1,
    }
    boot_config.update(agent_config)

    return boot_config

def initialize_custom_small(agent_config):
    # set up media
    media_id = 'custom'
    custom_media = {
        "ACET": 0.0,
        "CO+2": 100.0,
        "ETOH": 0.0,
        "FORMATE": 0.0,
        "GLYCEROL": 0.0,
        "LAC": 0.0,
        "LCTS": 3.4034,
        "OXYGEN-MOLECULE": 100.0,
        "PI": 100.0,
        "PYR": 0.0,
        "RIB": 0.0,
        "SUC": 0.0,
        "G6P": 1.3451,
        "GLC": 12.2087}

    new_media = {media_id: custom_media}
    timeline_str = '0 {}, 1200 end'.format(media_id)

    boot_config = {
        'timeline_str': timeline_str,
        'new_media': new_media,
        'run_for': 2.0,
        'depth': 1e-01, # 3000 um is default
        'edge_length_x': 1.0,
        'patches_per_edge_x': 1,
    }
    boot_config.update(agent_config)

    return boot_config

def initialize_glc_g6p(agent_config):
    timeline_str = '0 GLC_G6P, 3600 end'
    boot_config = {
        'timeline_str': timeline_str,
        'emit_fields': ['GLC', 'G6P']
    }
    boot_config.update(agent_config)
    return boot_config

def initialize_glc_lct(agent_config):
    timeline_str = '0 GLC_LCT, 3600 end'
    boot_config = {
        'timeline_str': timeline_str,
        'emit_fields': ['GLC', 'LCTS']
    }
    boot_config.update(agent_config)
    return boot_config

def initialize_glc_lct_shift(agent_config):
    timeline_str = '0 GLC_G6P, 1800 GLC_LCT, 3600 end'
    boot_config = {'timeline_str': timeline_str}
    boot_config.update(agent_config)
    return boot_config

def initialize_ecoli_core_glc(agent_config):
    timeline_str = '0 ecoli_core_GLC, 3600 end'
    boot_config = {
        'timeline_str': timeline_str,
        'emit_fields': ['glc__D_e']
    }
    boot_config.update(agent_config)
    return boot_config

def initialize_measp(agent_config):
    media_id = 'MeAsp_media'
    media = {'GLC': 20.0,  # assumes mmol/L
             'MeAsp': 1.0}
    new_media = {media_id: media}
    timeline_str = '0 {}, 3600 end'.format(media_id)  # (2hr*60*60 = 7200 s), (7hr*60*60 = 25200 s)
    boot_config = {
        'new_media': new_media,
        'timeline_str': timeline_str,
        'emit_fields': ['MeAsp'],
        'run_for': 1.0,
        'static_concentrations': True,
        'gradient': {
            'type': 'gaussian',
            'molecules': {
                'GLC': {
                    'center': [0.5, 0.5],
                    'deviation': 30.0},
                'MeAsp': {
                    'center': [0.5, 0.5],
                    'deviation': 30.0}
            }},
        'diffusion': 0.0,
        'rotation_jitter': 0.005,
        'edge_length_x': 50.0,
        'patches_per_edge_x': 40}
    boot_config.update(agent_config)

    return boot_config

def initialize_measp_long(agent_config):
    media_id = 'MeAsp_media'
    media = {'GLC': 20.0,  # assumes mmol/L
             'MeAsp': 1.0}
    new_media = {media_id: media}
    timeline_str = '0 {}, 3600 end'.format(media_id)  # (2hr*60*60 = 7200 s), (7hr*60*60 = 25200 s)
    boot_config = {
        'new_media': new_media,
        'timeline_str': timeline_str,
        'emit_fields': ['GLC','MeAsp'],
        'run_for': 1.0,
        'static_concentrations': True,
        'gradient': {
            'type': 'linear',
            'molecules': {
                'GLC': {
                    'center': [0.0, 0.0],
                    'slope': -1.0/150.0},
                'MeAsp': {
                    'center': [1.0, 1.0],
                    'slope': -1.0/150.0}
            }},
        'diffusion': 0.0,
        'rotation_jitter': 0.005,
        'edge_length_x': 200.0,
        'edge_length_y': 50.0,
        'patches_per_edge_x': 40}
    boot_config.update(agent_config)

    return boot_config


def initialize_measp_large(agent_config):

    ## Make media: GLC_G6P with MeAsp
    # get GLC_G6P media
    make_media = Media()
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

    emit_field = ['GLC', 'MeAsp']

    boot_config = {
        'timeline_str': timeline_str,
        'new_media': new_media,
        'run_for': 1.0,
        'cell_placement': [0.5, 0.5],  # place cells at center of lattice
        'emit_fields': emit_field,
        'static_concentrations': False,
        'gradient': {
            'type': 'linear',
            'molecules': {
                'GLC': {
                    'center': [0.0, 0.0],
                    'slope': -1.0 / 250.0},
                'MeAsp': {
                    'center': [1.0, 1.0],
                    'slope': -1.0 / 250.0}
            }},
        # 'diffusion': 0.001,
        'rotation_jitter': 0.005,
        'edge_length_x': 200.0,
        'patches_per_edge_x': 50}
    boot_config.update(agent_config)

    return boot_config

def initialize_measp_timeline(agent_config):
    # Endres and Wingreen (2006) use + 100 uM = 0.1 mmol for attractant. 0.2 b/c of dilution.
    timeline_str = '0 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '200 GLC 20.0 mmol 1 L + MeAsp 0.2 mmol 1 L, ' \
                   '600 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '1000 GLC 20.0 mmol 1 L + MeAsp 0.02 mmol 1 L, ' \
                   '1400 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '1600 end'

    boot_config = {
        'timeline_str': timeline_str,
        'run_for': 1.0,
        'static_concentrations': True,
        'diffusion': 0.0,
        'rotation_jitter': 0.005,
        'edge_length_x': 100.0,
        'patches_per_edge_x': 50,
    }
    boot_config.update(agent_config)

    return boot_config


class BootEnvironment(BootAgent):
    def __init__(self):
        super(BootEnvironment, self).__init__()
        self.agent_types = {
            # environments
            'lattice': wrap_boot_environment(initialize_lattice),
            'sugar1': wrap_boot_environment(initialize_glc_g6p),
            'sugar1_small': wrap_boot_environment(initialize_glc_g6p_small),
            'sugar2': wrap_boot_environment(initialize_glc_lct),
            'sugar_shift': wrap_boot_environment(initialize_glc_lct_shift),
            'ecoli_core_glc': wrap_boot_environment(initialize_ecoli_core_glc),
            'custom': wrap_boot_environment(initialize_custom_small),
            'measp': wrap_boot_environment(initialize_measp),
            'measp_long': wrap_boot_environment(initialize_measp_long),
            'measp_large': wrap_boot_environment(initialize_measp_large),
            'measp_timeline': wrap_boot_environment(initialize_measp_timeline),

            # basic compartments
            'lookup': wrap_boot(wrap_init_basic(TransportLookup), {'volume': 1.0}),
            'metabolism': wrap_boot(wrap_init_basic(BiGGMetabolism), {'volume': 1.0}),
            'regulation': wrap_boot(wrap_init_basic(Regulation), {'volume': 1.0}),
            'transport': wrap_boot(wrap_init_basic(Transport), {'volume': 1.0}),
            'growth': wrap_boot(wrap_init_basic(Growth), {'volume': 1.0}),
            'expression': wrap_boot(wrap_init_basic(ProteinExpression), {'volume': 1.0}),
            'receptor': wrap_boot(wrap_init_basic(ReceptorCluster), {'volume': 1.0}),
            'motor': wrap_boot(wrap_init_basic(MotorActivity), {'volume': 1.0}),
            'membrane_potential': wrap_boot(wrap_init_basic(MembranePotential), {'volume': 1.0}),

            # composite compartments
            'growth_division': wrap_boot(wrap_init_composite(compose_growth_division), {'volume': 1.0}),
            'chemotaxis': wrap_boot(wrap_init_composite(compose_simple_chemotaxis), {'volume': 1.0}),
            'pmf_chemotaxis': wrap_boot(wrap_init_composite(compose_pmf_chemotaxis), {'volume': 1.0}),
            'iFBA_ecoli_core': wrap_boot(wrap_init_composite(compose_iFBA), {'volume': 1.0}),
            }

def run():
    boot = BootEnvironment()
    boot.execute()

if __name__ == '__main__':
    run()