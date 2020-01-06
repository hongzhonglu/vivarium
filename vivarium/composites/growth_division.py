from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.derive_volume import DeriveVolume
from vivarium.processes.growth import Growth
from vivarium.processes.division import Division, divide_condition, divide_state
from vivarium.processes.protein_expression import ProteinExpression
from vivarium.processes.convenience_kinetics import ConvenienceKinetics


# kinetics for transport
def get_transport_config():
    # stoichiometry needs to match metabolism -- perhaps this can be extracted?
    transport_reactions = {
        'GLC_transport': {
            'stoichiometry': {
                'GLC_external': -1.0,
                'GLC_internal': 1.0},
            'is reversible': False,
            'catalyzed by': ['transporter_internal']}}

    # very simplified PTS
    transport_kinetics = {
        'GLC_transport': {
            'transporter_internal': {
                'GLC_external': 6.0,  # km for GLC
                'kcat_f': 1e-2}}}

    transport_initial_state = {
        'internal': {
            'GLC': 1.0,
            'transporter': 1.0},
        'external': {
            'GLC': 12.0},
        'fluxes': {
            'GLC': 1.0}}

    transport_roles = {
        'internal': ['GLC', 'transporter'],
        'external': ['GLC'],
    }

    return {
        'reactions': transport_reactions,
        'kinetic_parameters': transport_kinetics,
        'initial_state': transport_initial_state,
        'roles': transport_roles}



def compose_growth_division(config):

    # declare the processes
    transport_config = get_transport_config()
    transport = ConvenienceKinetics(transport_config)
    growth = Growth(config)
    division = Division(config)
    expression = ProteinExpression(config)
    deriver = DeriveVolume(config)

    # place processes in layers
    processes = [
        {'transport': transport,
         'growth': growth,
         'expression': expression},
        {'deriver': deriver,
         'division': division}]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'exchange',
            'fluxes': 'null'},
        'growth': {
            'internal': 'cell'},
        'division': {
            'internal': 'cell'},
        'expression': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # initialize the states
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

    out_dir = os.path.join('out', 'tests', 'growth_division_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_growth_division)

    # settings for simulation and plot
    options = compose_growth_division({})['options']
    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 20}

    plot_settings = {
        'max_rows': 25}

    # saved_state = simulate_compartment(compartment, settings)
    saved_state = simulate_with_environment(compartment, settings)
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)