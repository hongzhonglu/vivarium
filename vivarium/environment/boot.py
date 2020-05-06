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
import copy

from vivarium.actor.inner import Inner
from vivarium.actor.outer import Outer
from vivarium.actor.boot import BootAgent
from vivarium.actor.control import DEFAULT_EMITTER_CONFIG
from vivarium.compartment.emitter import get_emitter, configure_emitter

# environment
from vivarium.environment.lattice import EnvironmentSpatialLattice
from vivarium.environment.lattice_compartment import LatticeCompartment, generate_lattice_compartment
from vivarium.environment.make_media import Media
from vivarium.utils.units import units

# processes
from vivarium.processes.transport_lookup import TransportLookup
from vivarium.processes.metabolism import Metabolism
from vivarium.processes.Kremling2007_transport import Transport
from vivarium.processes.growth import Growth
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Vladimirov2008_motor import MotorActivity
from vivarium.processes.membrane_potential import MembranePotential
from vivarium.processes.antibiotic_transport import AntibioticTransport
from vivarium.processes.transcription import Transcription

# composites
from vivarium.composites.master import compose_master
from vivarium.composites.glc_lct_shifter import compose_glc_lct_shifter
from vivarium.composites.growth_division import compose_growth_division
from vivarium.composites.chemotaxis_minimal import compose_simple_chemotaxis
from vivarium.composites.antibiotics import (
    compose_antibiotics,
)


DEFAULT_COLOR = [0.6, 0.4, 0.3]

def wrap_boot(initialize, initial_state):
    def boot(agent_id, agent_type, actor_config):
        initial_state.update(actor_config.get('declare', {}))
        actor_config['declare'] = initial_state  # 'declare' is for the environment

        return Inner(
            agent_id,
            agent_type,
            actor_config,
            initialize)

    return boot

def wrap_init_basic(make_process):
    def initialize(boot_config):
        process = make_process(boot_config)  # 'boot_config', set in environment.control is the process' initial_parameters
        return generate_lattice_compartment(process, boot_config)

    return initialize

def wrap_init_composite(make_composite):
    def initialize(boot_config):
        # set up the the composite
        composite_config = make_composite(boot_config)
        processes = composite_config['processes']
        derivers = composite_config['derivers']
        states = composite_config['states']
        options = composite_config['options']
        topology = options['topology']

        # update options
        emitter = configure_emitter(boot_config, processes, topology)
        options.update({'emitter': emitter})

        # create the compartment
        return LatticeCompartment(processes, derivers, states, options)

    return initialize

def wrap_boot_environment(intialize):
    def boot(agent_id, agent_type, actor_config):
        boot_config = copy.deepcopy(actor_config['boot_config'])

        # get boot_config from initialize
        boot_config = intialize(boot_config)

        # paths
        working_dir = actor_config.get('working_dir', os.getcwd())
        output_dir = os.path.join(working_dir, 'out', agent_id)
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

        emitter_config = boot_config.get('emitter', DEFAULT_EMITTER_CONFIG)
        emitter_config['experiment_id'] = agent_id

        emitter = get_emitter(emitter_config)
        boot_config.update({
            'emitter': emitter,
            'output_dir': output_dir})

        # create the environment
        environment = EnvironmentSpatialLattice(boot_config)

        return EnvironmentActor(agent_id, agent_type, actor_config, environment)

    return boot

class EnvironmentActor(Outer):
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


# Define environment initialization functions
def initialize_lattice(boot_config):

    lattice_config = {
        'name': 'lattice',
        'description': 'a standard lattice environment'}

    # set up media
    media_id = boot_config.get('media_id', 'minimal')
    media = boot_config.get('media', {})
    if media:
        lattice_config.update({'concentrations': media})
    else:
        lattice_config.update({'media_id': media_id})

    lattice_config.update(boot_config)
    return lattice_config

def initialize_glc_g6p_small(boot_config):
    # set up media
    media_id = 'GLC_G6P'
    timeline_str = '0 {}, 1800 end'.format(media_id)  # (2hr*60*60 = 7200 s), (7hr*60*60 = 25200 s)
    lattice_config = {
        'name': 'glc_g6p_small',
        'timeline_str': timeline_str,
        'run_for': 2.0,
        'depth': 1e-01, # 3000 um is default
        'edge_length_x': 1.0,
        'patches_per_edge_x': 1,
    }

    lattice_config.update(boot_config)
    return lattice_config

def initialize_custom_small(boot_config):
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

    lattice_config = {
        'timeline_str': timeline_str,
        'new_media': new_media,
        'run_for': 2.0,
        'depth': 1e-01, # 3000 um is default
        'edge_length_x': 1.0,
        'patches_per_edge_x': 1,
    }

    lattice_config.update(boot_config)
    return lattice_config

def initialize_glc_g6p(boot_config):
    timeline_str = '0 GLC_G6P, 3600 end'
    lattice_config = {
        'name': 'glc_g6p',
        'timeline_str': timeline_str,
        'emit_fields': ['GLC', 'G6P']
    }

    lattice_config.update(boot_config)
    return lattice_config

def initialize_glc_lct(boot_config):
    timeline_str = '0 GLC_LCT, 3600 end'
    lattice_config = {
        'name': 'glc_lct',
        'timeline_str': timeline_str,
        'emit_fields': ['GLC', 'LCTS']
    }

    lattice_config.update(boot_config)
    return lattice_config

def initialize_glc_lct_shift(boot_config):
    timeline_str = '0 GLC_G6P, 1800 GLC_LCT, 3600 end'
    lattice_config = {
        'name': 'glc_lct_shift',
        'timeline_str': timeline_str}

    lattice_config.update(boot_config)
    return lattice_config

def initialize_ecoli_core_glc(boot_config):
    timeline_str = '0 ecoli_core_GLC 1.0 L + lcts_e 1.0 mmol 0.1 L, 21600 end'

    lattice_config = {
        'name': 'ecoli_core',
        'timeline_str': timeline_str,
        'edge_length_x': 15.0,
        'patches_per_edge_x': 10,
        'run_for': 5.0,
        'diffusion': 1e-3,
        'depth': 2e-2,
        'translation_jitter': 1e-1,
        'emit_fields': [
            'co2_e',
            'o2_e',
            'glc__D_e',
            'lcts_e']}

    lattice_config.update(boot_config)
    return lattice_config

def initialize_measp(boot_config):
    media_id = 'MeAsp_media'
    media = {'GLC': 20.0,  # assumes mmol/L
             'MeAsp': 80.0}
    new_media = {media_id: media}
    timeline_str = '0 {}, 600 end'.format(media_id)
    lattice_config = {
        'name': 'measp',
        'new_media': new_media,
        'timeline_str': timeline_str,
        'emit_fields': ['MeAsp'],
        'run_for': 0.1,
        'static_concentrations': True,
        'emit_frequency': 20,
        'gradient': {
            'type': 'linear',
            'molecules': {
                'GLC': {
                    'center': [0.5, 0.5],
                    'slope': -1e-2},
                'MeAsp': {
                    'center': [0.5, 0.5],
                    'slope': -1e-2}
            }},
        'diffusion': 0.0,
        'edge_length_x': 50.0,
        'patches_per_edge_x': 40}

    lattice_config.update(boot_config)
    return lattice_config

def initialize_measp_long(boot_config):
    media_id = 'MeAsp_media'
    media = {'GLC': 0.01,  # assumes mmol/L
             'MeAsp': 0.01}
    new_media = {media_id: media}
    timeline_str = '0 {}, 100 end'.format(media_id)
    lattice_config = {
        'name': 'measp_long',
        'new_media': new_media,
        'timeline_str': timeline_str,
        'emit_fields': [],  # don't emit fields due to large lattice size
        'run_for': 0.1,  # high coupling between cell and env requires short exchange timestep
        'static_concentrations': True,
        'cell_placement': [0.1, 0.5],  # place cells at bottom of gradient
        'gradient': {
            'type': 'exponential',
            'molecules': {
                'GLC': {
                    'center': [0.0, 0.5],
                    'base': 1+6e-4},
                'MeAsp': {
                    'center': [0.0, 0.5],
                    'base': 1+6e-4}
            }},
        'jitter_force': 1e-2,
        'edge_length_x': 4000.0,
        'edge_length_y': 800.0,
        'patches_per_edge_x': 4000  # need high resolution for receptors to detect change
    }

    lattice_config.update(boot_config)
    return lattice_config

def initialize_antibiotic(boot_config):
    media_id = 'antibiotic_media'
    media = {
        'antibiotic': 1.0,  # mmol/L
    }
    new_media = {media_id: media}
    lattice_config = {
        'name': 'measp_long',
        'new_media': new_media,
        'timeline_str': '0 {}, 7200 end'.format(media_id),
        'emit_fields': ['antibiotic'],
        'run_for': 1.0,  # timestep, in sec
        'rotation_jitter': 0.005,
        'edge_length_x': 10.0,  # um
        'edge_length_y': 10.0,  # um
        'depth': 1.0,  # um
        # Patches discretize space for diffusion
        'patches_per_edge_x': 1,
    }

    boot_config.update(lattice_config)
    return boot_config


def initialize_measp_large(boot_config):
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
        'cell_placement': [0.5, 0.5],  # place cells at center of lattice
        'timeline_str': timeline_str,
        'new_media': new_media,
        'run_for': 0.1,
        'emit_fields': ['MeAsp', 'GLC'],
        'static_concentrations': False,
        'diffusion': 0.001,
        'edge_length_x': 200.0,
        'edge_length_y': 200.0,
        'patches_per_edge_x': 50}

    lattice_config.update(boot_config)
    return lattice_config

def initialize_measp_timeline(boot_config):
    # Endres and Wingreen (2006) use + 100 uM = 0.1 mmol for attractant. 0.2 b/c of dilution.
    timeline_str = '0 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '200 GLC 20.0 mmol 1 L + MeAsp 0.2 mmol 1 L, ' \
                   '600 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '1000 GLC 20.0 mmol 1 L + MeAsp 0.02 mmol 1 L, ' \
                   '1400 GLC 20.0 mmol 1 L + MeAsp 0.0 mmol 1 L, ' \
                   '1600 end'

    lattice_config = {
        'name': 'measp_timeline',
        'timeline_str': timeline_str,
        'run_for': 1.0,
        'static_concentrations': True,
        'diffusion': 0.0,
        'edge_length_x': 100.0,
        'patches_per_edge_x': 50,
    }

    lattice_config.update(boot_config)
    return lattice_config


class BootEnvironment(BootAgent):
    def __init__(self):
        super(BootEnvironment, self).__init__()
        self.agent_types = {
            # environments
            'lattice': wrap_boot_environment(initialize_lattice),
            'glc_g6p': wrap_boot_environment(initialize_glc_g6p),
            'glc_g6p_small': wrap_boot_environment(initialize_glc_g6p_small),
            'glc_lct': wrap_boot_environment(initialize_glc_lct),
            'glc_lct_shift': wrap_boot_environment(initialize_glc_lct_shift),
            'ecoli_core_glc': wrap_boot_environment(initialize_ecoli_core_glc),
            'custom': wrap_boot_environment(initialize_custom_small),
            'measp': wrap_boot_environment(initialize_measp),
            'measp_long': wrap_boot_environment(initialize_measp_long),
            'measp_large': wrap_boot_environment(initialize_measp_large),
            'measp_timeline': wrap_boot_environment(initialize_measp_timeline),
            'antibiotic_environment': wrap_boot_environment(initialize_antibiotic),

            # basic compartments
            'lookup': wrap_boot(wrap_init_basic(TransportLookup), {'volume': 1.0}),
            'metabolism': wrap_boot(wrap_init_basic(Metabolism), {'volume': 1.0}),
            'transport': wrap_boot(wrap_init_basic(Transport), {'volume': 1.0}),
            'growth': wrap_boot(wrap_init_basic(Growth), {'volume': 1.0}),
            'expression': wrap_boot(wrap_init_basic(MinimalExpression), {'volume': 1.0}),
            'receptor': wrap_boot(wrap_init_basic(ReceptorCluster), {'volume': 1.0}),
            'motor': wrap_boot(wrap_init_basic(MotorActivity), {'volume': 1.0}),
            'membrane_potential': wrap_boot(wrap_init_basic(MembranePotential), {'volume': 1.0}),
            'antibiotic_transport': wrap_boot(
                wrap_init_basic(AntibioticTransport), {'volume': 1.0}),
            'transcription': wrap_boot(wrap_init_basic(Transcription), {'volume': 1.0}),

            # composite compartments
            'master': wrap_boot(wrap_init_composite(compose_master), {'volume': 1.0}),
            'shifter': wrap_boot(wrap_init_composite(compose_glc_lct_shifter), {'volume': 1.0}),
            'growth_division': wrap_boot(wrap_init_composite(compose_growth_division), {'volume': 1.0}),
            'minimal_chemotaxis': wrap_boot(wrap_init_composite(compose_simple_chemotaxis), {'volume': 1.0}),
            'antibiotic_composite': wrap_boot(
                wrap_init_composite(compose_antibiotics),
                {'volume': 1.0},
            ),
        }

def run():
    boot = BootEnvironment()
    boot.execute()

if __name__ == '__main__':
    run()
