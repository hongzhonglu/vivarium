from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.composition import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_with_environment,
    get_derivers,
)
from vivarium.compartment.process import (
    initialize_state,
    COMPARTMENT_STATE,
    load_compartment,
)
from vivarium.processes.antibiotics import Antibiotics
from vivarium.processes.antibiotics import (
    DEFAULT_INITIAL_STATE as ANTIBIOTIC_DEFAULT_INITIAL_STATE,
)
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.division import (
    Division,
    divide_condition,
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

    # Death Config
    checkers_config = config.setdefault('checkers', {})
    antibiotic_checker_config = checkers_config.setdefault(
        'antibiotic', {})
    antibiotic_checker_config.setdefault('antibiotic_threshold', 0.09)

    # declare the processes
    antibiotic_transport = Antibiotics(config)
    growth = Growth(config)
    expression = ODE_expression(config)
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
            'division': division,
        },
    ]

    # make the topology.
    # for each process, map process ports to compartment ports
    topology = {
        'antibiotic_transport': {
            'internal': 'cell',
            'external': 'environment',
            'exchange': 'exchange',
            'fluxes': 'fluxes',
            'global': 'global',
        },
        'growth': {
            'global': 'global',
        },
        'expression': {
            'counts': 'cell_counts',
            'internal': 'cell',
            'external': 'environment',
        },
        'division': {
            'global': 'global',
        },
        'death': {
            'internal': 'cell',
            'compartment': COMPARTMENT_STATE,
        },
    }

    # Add Derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])
    topology.update(derivers['deriver_topology'])

    # initialize the states
    states = initialize_state(
        processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'antibiotic_growth_composite',
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
    }

    return {
        'processes': processes,
        'states': states,
        'options': options}


def test_antibiotic_growth_composite():
    options = compose_antibiotic_growth({})['options']
    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
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
    out_dir = os.path.join('out', 'tests', 'antibiotic_growth_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {
        'max_rows': 25,
        'skip_ports': ['prior_state'],
    }

    saved_state = test_antibiotic_growth_composite()
    del saved_state[0]  # Delete first record, where everything is 0
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)
