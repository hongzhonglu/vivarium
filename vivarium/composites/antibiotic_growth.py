from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.antibiotics import Antibiotics
from vivarium.processes.antibiotics import (
    DEFAULT_INITIAL_STATE as ANTIBIOTIC_DEFAULT_INITIAL_STATE,
)
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.deriver import Deriver
from vivarium.processes.division import (
    Division,
    divide_condition,
    divide_state,
)
from vivarium.processes.growth import Growth
from vivarium.processes.ode_expression import ODE_expression


def compose_antibiotic_growth(config):

    # Expression Config
    transcription_config = config.setdefault('transcription_rates', {})
    transcription_config.setdefault('antibiotic_exporter_RNA', 1e-3)
    translation_config = config.setdefault('translation_rates', {})
    translation_config.setdefault('antibiotic_exporter', 1.0)
    degradation_config = config.setdefault('degradation_rates', {})
    degradation_config.setdefault('antibiotic_exporter', 1.0)
    degradation_config.setdefault('antibiotic_exporter_RNA', 1e-3)
    protein_map = config.setdefault('protein_map', {})
    protein_map.setdefault(
        'antibiotic_exporter', 'antibiotic_exporter_RNA')

    initial_state_config = config.setdefault(
        'initial_state', ANTIBIOTIC_DEFAULT_INITIAL_STATE)
    internal_initial_config = initial_state_config.setdefault(
        'internal', {})
    internal_initial_config['antibiotic_exporter'] = 0.0

    # Deriver Config
    counted_molecules = config.setdefault('counted_molecules', [])
    if 'antibiotic_exporter' not in counted_molecules:
        counted_molecules.append('antibiotic_exporter')

    # Death Config
    checkers_config = config.setdefault('checkers', {})
    antibiotic_checker_config = checkers_config.setdefault(
        'antibiotic', {})
    antibiotic_checker_config.setdefault('antibiotic_threshold', 0.09)

    # declare the processes
    antibiotic_transport = Antibiotics(config)
    growth = Growth(config)
    expression = ODE_expression(config)
    deriver = Deriver(config)
    death = DeathFreezeState(config)
    division = Division(config)

    # place processes in layers
    processes = [
        {
            'antibiotic_transport': antibiotic_transport,
            'growth': growth,
            'expression': expression,
            'death': death,
        },
        {
            'deriver': deriver,
            'division': division,
        },
    ]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'antibiotic_transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'exchange',
            'fluxes': 'fluxes',
        },
        'growth': {
            'internal': 'cell'
        },
        'expression': {
            'counts': 'cell_counts',
            'internal': 'cell',
            'external': 'environment',
        },
        'deriver': {
            'state': 'cell',
            'counts': 'cell_counts',
            'prior_state': 'prior_state',
        },
        'division': {
            'internal': 'cell',
        },
        'death': {
            'internal': 'cell',
        },
    }

    # initialize the states
    states = initialize_state(
        processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'antibiotic_growth_composite',
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
        'divide_state': divide_state,
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}


def test_antibiotic_growth_composite():
    options = compose_antibiotic_growth({})['options']
    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-6,  # L
        'timestep': 1,
        'total_time': 2000,
    }
    config = {
        'transcription_rates': {
            'antibiotic_exporter_RNA': 1e-3,
        },
        'degradation_rates': {
            # Set for on the order of 100 RNAs at equilibrium
            'antibiotic_exporter_RNA': 1.0,
            # Set so exporter concentration reaches equilibrium
            'antibiotic_exporter': 1e-3,
        },
        'emitter': 'null',
        'checkers': {
            'antibiotic': {
                # Set so cell dies after first division
                'antibiotic_threshold': 0.09,
            },
        },
    }
    compartment = load_compartment(compose_antibiotic_growth, config)
    saved_state = simulate_with_environment(compartment, settings)
    return saved_state


if __name__ == '__main__':
    from vivarium.actor.process import (
        load_compartment,
        convert_to_timeseries,
        plot_simulation_output,
        simulate_with_environment,
    )
    out_dir = os.path.join('out', 'tests', 'antibiotic_growth_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {
        'max_rows': 25,
        'skip_roles': ['prior_state'],
    }

    saved_state = test_antibiotic_growth_composite()
    del saved_state[0]  # Delete first record, where everything is 0
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)
