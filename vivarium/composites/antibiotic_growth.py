from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state

# processes
from vivarium.processes.deriver import Deriver
from vivarium.processes.growth import Growth
from vivarium.processes.minimal_expression import MinimalExpression
from vivarium.processes.antibiotics import Antibiotics
from vivarium.processes.antibiotics import (
    DEFAULT_INITIAL_STATE as ANTIBIOTIC_DEFAULT_INITIAL_STATE,
)


def compose_antibiotic_growth(config):

    # Expression Config
    expression_rates_config = config.setdefault('expression_rates', {})
    expression_rates_config.setdefault('antibiotic_exporter', 0.5)

    initial_state_config = config.setdefault(
        'initial_state', ANTIBIOTIC_DEFAULT_INITIAL_STATE)
    initial_state_config['internal']['antibiotic_exporter'] = 0.0

    # Deriver Config
    counted_molecules = config.setdefault('counted_molecules', [])
    if 'antibiotic_exporter' not in counted_molecules:
        counted_molecules.append('antibiotic_exporter')

    # declare the processes
    antibiotic_transport = Antibiotics(config)
    growth = Growth(config)
    expression = MinimalExpression(config)
    deriver = Deriver(config)

    # place processes in layers
    processes = [
        {
            'antibiotic_transport': antibiotic_transport,
            'growth': growth,
            'expression': expression,
        },
        {
            'deriver': deriver,
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
            'internal': 'cell_counts',
            'external': 'environment',
        },
        'deriver': {
            'state': 'cell',
            'counts': 'cell_counts',
            'prior_state': 'prior_state',
        },
    }

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'antibiotic_growth_composite',
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
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
        'total_time': 20,
    }
    compartment = load_compartment(compose_antibiotic_growth, {})
    saved_state = simulate_with_environment(compartment, settings)
    return saved_state


if __name__ == '__main__':
    from vivarium.actor.process import load_compartment, convert_to_timeseries, plot_simulation_output, simulate_with_environment

    out_dir = os.path.join('out', 'tests', 'antibiotic_growth_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {
        'max_rows': 25,
        'skip_roles': ['prior_state', 'cell_counts'],
    }

    saved_state = test_antibiotic_growth_composite()
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)
