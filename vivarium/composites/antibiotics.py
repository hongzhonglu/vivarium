from __future__ import absolute_import, division, print_function

import math
import os

from vivarium.compartment.composition import (
    plot_simulation_output,
    simulate_with_environment,
    get_derivers,
    flatten_timeseries,
    save_timeseries,
    load_timeseries,
    REFERENCE_DATA_DIR,
    TEST_OUT_DIR,
    assert_timeseries_close,
    load_compartment)
from vivarium.compartment.process import (
    initialize_state,
    COMPARTMENT_STATE,
)
from vivarium.processes.antibiotic_transport import AntibioticTransport
from vivarium.processes.antibiotic_transport import (
    DEFAULT_INITIAL_STATE as ANTIBIOTIC_DEFAULT_INITIAL_STATE,
)
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.division import (
    Division,
    divide_condition,
)
from vivarium.processes.growth import Growth
from vivarium.processes.ode_expression import ODE_expression


NAME = 'antibiotics_composite'
NUM_DIVISIONS = 3
DIVISION_TIME = 2400  # seconds to divide


def compose_antibiotics(config):

    division_time = config.get('cell_cycle_division_time', 2400)

    # Expression Config
    transcription_config = config.setdefault('transcription_rates', {})
    transcription_config.setdefault('AcrAB-TolC_RNA', 1e-3)
    translation_config = config.setdefault('translation_rates', {})
    translation_config.setdefault('AcrAB-TolC', 1.0)
    degradation_config = config.setdefault('degradation_rates', {})
    degradation_config.setdefault('AcrAB-TolC', 1.0)
    degradation_config.setdefault('AcrAB-TolC_RNA', 1e-3)
    protein_map = config.setdefault('protein_map', {})
    protein_map.setdefault(
        'AcrAB-TolC', 'AcrAB-TolC_RNA')

    initial_state_config = config.setdefault(
        'initial_state', ANTIBIOTIC_DEFAULT_INITIAL_STATE)
    internal_initial_config = initial_state_config.setdefault(
        'internal', {})
    internal_initial_config['AcrAB-TolC'] = 0.0

    # Death Config
    checkers_config = config.setdefault('checkers', {})
    antibiotic_checker_config = checkers_config.setdefault(
        'antibiotic', {})
    antibiotic_checker_config.setdefault('antibiotic_threshold', 0.09)

    # Growth Config
    # Growth rate calculated so that 2 = exp(DIVISION_TIME * rate)
    # because division process divides once cell doubles in size
    config.setdefault('growth_rate', math.log(2) / division_time)

    # declare the processes
    antibiotic_transport = AntibioticTransport(config)
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
            'global': 'global',
        },
    }

    # add derivers
    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    all_processes = processes + derivers['deriver_processes']
    topology.update(derivers['deriver_topology'])

    # initialize the states
    states = initialize_state(
        all_processes,
        topology,
        config.get('initial_state', {}))

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
        'derivers': deriver_processes,
        'states': states,
        'options': options}


def run_antibiotics_composite():
    options = compose_antibiotics({})['options']
    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
        'environment_volume': 1e-5,  # L
        'timestep': 1,
        'total_time': DIVISION_TIME * NUM_DIVISIONS,
    }
    config = {
        'transcription_rates': {
            'AcrAB-TolC_RNA': 1e-3,
        },
        'degradation_rates': {
            # Set for on the order of 100 RNAs at equilibrium
            'AcrAB-TolC_RNA': 1.0,
            # Set so exporter concentration reaches equilibrium
            'AcrAB-TolC': 1e-3,
        },
        'checkers': {
            'antibiotic': {
                # Set so cell dies after first division
                'antibiotic_threshold': 10.0,
            },
        },
    }
    compartment = load_compartment(compose_antibiotics, config)
    saved_state = simulate_with_environment(compartment, settings)
    return saved_state


def test_antibiotics_composite_similar_to_reference():
    timeseries = run_antibiotics_composite()
    flattened = flatten_timeseries(timeseries)
    reference = load_timeseries(
        os.path.join(REFERENCE_DATA_DIR, NAME + '.csv'))
    assert_timeseries_close(flattened, reference)


def main():
    out_dir = os.path.join(TEST_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {
        'max_rows': 25,
        'skip_ports': ['prior_state'],
    }

    timeseries = run_antibiotics_composite()
    plot_simulation_output(timeseries, plot_settings, out_dir)
    save_timeseries(timeseries, out_dir)


if __name__ == '__main__':
    main()
