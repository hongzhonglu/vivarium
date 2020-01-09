from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.actor.process import initialize_state
import vivarium.utils.regulation_logic as rl

# processes
from vivarium.processes.derive_volume import DeriveVolume
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.BiGG_metabolism import BiGGMetabolism
from vivarium.processes.convenience_kinetics import ConvenienceKinetics


# BiGG model for metabolism
METABOLISM_FILE = os.path.join('models', 'e_coli_core.json')

# convenience kinetics configuration for transport
def get_transport_config():
    # stoichiometry needs to match metabolism -- perhaps this can be extracted?
    transport_reactions = {
        'GLCpts': {
            'stoichiometry': {
                'g6p_c_internal': 1.0,
                'glc__D_e_external': -1.0,
                'pep_c_internal': -1.0,
                'pyr_c_internal': 1.0},
            'is reversible': False,
            'catalyzed by': ['PTSG_internal']}}

    # very simplified PTS
    transport_kinetics = {
        'GLCpts': {
            'PTSG_internal': {
                'glc__D_e_external': 1e1,
                'pep_c_internal': 1e1,
                'kcat_f': 4e-1}}}

    transport_initial_state = {
        'internal': {
            'g6p_c': 1.0,
            'pep_c': 1.0,
            'pyr_c': 1.0,
            'PTSG': 1.0},
        'external': {
            'glc__D_e': 12.0},
        'fluxes': {
            'GLCpts': 1.0}}  # TODO -- initial fluxes should be set in kinetics process

    transport_roles = {
        'internal': ['g6p_c', 'pep_c', 'pyr_c', 'PTSG'],
        'external': ['glc__D_e'],
        }
    
    # return {
    #     'reactions': transport_reactions,
    #     'kinetic_parameters': transport_kinetics,
    #     'initial_state': transport_initial_state,
    #     'roles': transport_roles}
    return {
        'reactions': {},
        'kinetic_parameters': {},
        'initial_state': {
            'internal': {},
            'external': {}},
        'roles': {
            'internal': [],
            'external': []}}

def get_regulation():
    regulation = {
        'EX_lac__D_e': rl.build_rule('IF not (glc__D_e_external)'),
    }
    return regulation

# the composite function
def compose_kinetic_FBA(config):

    ## Declare the processes
    # ordering allows earlier processes to inform the configuration of later processes

    ## Transport
    transport_config = copy.deepcopy(config)
    transport_config.update(get_transport_config())
    transport = ConvenienceKinetics(transport_config)

    ## Metabolism
    # get target fluxes from transport, load in regulation function
    metabolism_config = copy.deepcopy(config)
    target_fluxes = transport.kinetic_rate_laws.reaction_ids
    regulation = get_regulation()

    metabolism_config.update({
        'model_path': METABOLISM_FILE,
        'constrained_reaction_ids': target_fluxes,
        'regulation': regulation})
    metabolism = BiGGMetabolism(metabolism_config)

    ## Division
    # get initial volume from metabolism
    division_config = copy.deepcopy(config)
    division_config.update({'initial_state': metabolism.initial_state})
    division = Division(division_config)

    # Other processes
    deriver = DeriveVolume(config)

    # Place processes in layers
    processes = [
        {'transport': transport},
        {'metabolism': metabolism},
        {'deriver': deriver,
        'division': division}
    ]

    ## Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds'},
        'metabolism': {
            'internal': 'cell',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds'},
        'division': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    ## Initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}


if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'kinetic_FBA_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # TODO -- load print emitter
    compartment = load_compartment(compose_kinetic_FBA)

    # settings for simulation and plot
    options = compose_kinetic_FBA({})['options']

    # define timeline
    timeline = [
        (0, {'environment': {
            'lac__D_e': 12.0}
        }),
        (500, {'environment': {
            'glc__D_e': 0.05}
        }),
        # (1000, {'environment': {
        #     'glc__D_e': 0.0000001}
        # }),
        # (3500, {})]
        (1000, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 25,
        'remove_zeros': True,
        'overlay': {'reactions': 'flux_bounds'},
        'show_state': [
            ('environment', 'glc__D_e'),
            ('environment', 'lac__D_e'),
            ('reactions', 'GLCpts'),
            ('reactions', 'EX_glc__D_e'),
            ('reactions', 'EX_lac__D_e'),
            ('cell', 'g6p_c'),
            ('cell', 'pep_c'),
            ('cell', 'pyr_c'),
            ('cell', 'PTSG'),
            ]}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
