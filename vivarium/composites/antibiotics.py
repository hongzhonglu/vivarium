from __future__ import absolute_import, division, print_function

import math
import os

from vivarium.compartment.composition import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_with_environment,
    get_derivers,
    flatten_timeseries,
    save_timeseries,
    load_timeseries,
    REFERENCE_DATA_DIR,
    TEST_OUT_DIR,
    assert_timeseries_close,
)
from vivarium.compartment.process import (
    initialize_state,
    COMPARTMENT_STATE,
    load_compartment,
)
from vivarium.data.chromosomes.antibiotic_export_chromosome import (
    AntibioticExportChromosome,
)
from vivarium.data.nucleotides import nucleotides
from vivarium.data.amino_acids import amino_acids
from vivarium.processes.antibiotic_transport import AntibioticTransport
from vivarium.processes.antibiotic_transport import (
    DEFAULT_INITIAL_STATE as ANTIBIOTIC_DEFAULT_INITIAL_STATE,
)
from vivarium.processes.complexation import Complexation
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.degradation import RnaDegradation
from vivarium.processes.division import (
    Division,
    divide_condition,
)
from vivarium.processes.growth import Growth
from vivarium.processes.injector import Injector
from vivarium.processes.transcription import (
    UNBOUND_RNAP_KEY,
    Transcription,
)
from vivarium.processes.translation import (
    UNBOUND_RIBOSOME_KEY,
    Translation,
)
from vivarium.states.chromosome import Chromosome
from vivarium.utils.dict_utils import deep_merge_check


NAME = 'antibiotics_composite'
NUM_DIVISIONS = 3
DIVISION_TIME = 2400  # seconds to divide
INITIAL_MOLECULE_COUNT = 1e6


def compose_antibiotics(config):
    chromosome = AntibioticExportChromosome()
    plasmid = Chromosome(chromosome.config)
    sequences = plasmid.product_sequences()

    division_time = config.get('cell_cycle_division_time', 2400)

    transcription_config = {
        'sequence': chromosome.config['sequence'],
        'templates': chromosome.config['promoters'],
        'genes': chromosome.config['genes'],
        'transcription_factors': chromosome.transcription_factors,
        'promoter_affinities': chromosome.promoter_affinities,
        'polymerase_occlusion': 30,
        'elongation_rate': 50,
    }
    transcription = Transcription(transcription_config)

    translation_config = {
        'sequences': chromosome.operon_sequences,
        'templates': chromosome.transcript_templates,
        'concentration_keys': ['marA', 'acrR'],
        'transcript_affinities': chromosome.transcript_affinities,
        'elongation_rate': 22,
        'polymerase_occlusion': 50,
    }
    translation = Translation(translation_config)

    degradation_config = {
        'sequences': sequences,
        'catalysis_rates': {
            'endoRNAse': 0.01},
        'degradation_rates': {
            'transcripts': {
                'endoRNAse': {
                    transcript: 1e-23
                    for transcript in chromosome.config['genes']
                },
            },
        },
    }
    degradation = RnaDegradation(degradation_config)

    complexation_config = {
        'stoichiometry': chromosome.complexation_stoichiometry,
        'monomer_ids': chromosome.complexation_monomer_ids,
        'complex_ids': chromosome.complexation_complex_ids,
        'rates': chromosome.complexation_rates,
    }
    complexation = Complexation(complexation_config)

    init_molecule_counts = {}
    for nucleotide in nucleotides.values():
        init_molecule_counts[nucleotide] = INITIAL_MOLECULE_COUNT
    for amino_acid in amino_acids.values():
        init_molecule_counts[amino_acid] = INITIAL_MOLECULE_COUNT

    initial_state_config = config.setdefault(
        'initial_state', ANTIBIOTIC_DEFAULT_INITIAL_STATE)
    internal_initial_config = initial_state_config.setdefault(
        'internal', {})
    internal_initial_config['acrAB-tolC'] = 0.0

    # Injector Config
    injector_config = {
        'substrate_rate_map': {
            nucleotide: 100.0
            for nucleotide in nucleotides.values()
        }
    }
    injector_config['substrate_rate_map'].update({
        aa: 5.0 for aa in amino_acids.values()
    })
    injector = Injector(injector_config)

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
    death = DeathFreezeState(config)
    division = Division(config)

    # place processes in layers
    processes = [
        {
            'transcription': transcription,
            'translation': translation,
            'degradation': degradation,
            'complexation': complexation,
            'injector': injector,
        },
        {
            'antibiotic_transport': antibiotic_transport,
            'growth': growth,
            'death': death,
        },
        {
            'division': division,
        },
    ]

    # make the topology.
    # for each process, map process ports to compartment ports
    topology = {
        'transcription': {
            'chromosome': 'chromosome',
            'molecules': 'molecules',
            'proteins': 'cell',
            'transcripts': 'transcripts',
            'factors': 'concentrations',
        },
        'translation': {
            'ribosomes': 'ribosomes',
            'molecules': 'molecules',
            'transcripts': 'transcripts',
            'proteins': 'cell',
            'concentrations': 'concentrations',
        },
        'degradation': {
            'transcripts': 'transcripts',
            'proteins': 'cell',
            'molecules': 'molecules',
            'global': 'global',
        },
        'complexation': {
            'monomers': 'cell',
            'complexes': 'cell',
        },
        'injector': {
            'internal': 'molecules',
        },
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
        'division': {
            'global': 'global',
        },
        'death': {
            'internal': 'cell',
            'compartment': COMPARTMENT_STATE,
            'global': 'global',
        },
    }

    # Add Derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])
    topology.update(derivers['deriver_topology'])

    initial_state = {
        'molecules': init_molecule_counts,
        'transcripts': {
            gene: 0
            for gene in chromosome.config['genes']
        },
        'proteins': {
            'marA': 0,
            'acrR': 0,
            'endoRNAse': 1,
            UNBOUND_RIBOSOME_KEY: 10,
            UNBOUND_RNAP_KEY: 10,
        },
    }
    deep_merge_check(initial_state, config.get('initial_state', {}))

    # initialize the states
    states = initialize_state(
        processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'antibiotics_composite',
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


def run_antibiotics_composite(total_time):
    options = compose_antibiotics({})['options']
    settings = {
        'environment_port': options['environment_port'],
        'exchange_port': options['exchange_port'],
        'environment_volume': 1e-5,  # L
        'timestep': 1,
        'total_time': total_time,
    }
    config = {
        'transcription_rates': {
            'acrAB-tolC_RNA': 1e-3,
        },
        'degradation_rates': {
            # Set for on the order of 100 RNAs at equilibrium
            'acrAB-tolC_RNA': 1.0,
            # Set so exporter concentration reaches equilibrium
            'acrAB-tolC': 1e-3,
        },
        'emitter': 'null',
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
    saved_data = run_antibiotics_composite(100)
    timeseries = convert_to_timeseries(saved_data)
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

    saved_state = run_antibiotics_composite(
        DIVISION_TIME * NUM_DIVISIONS)
    del saved_state[0]  # Delete first record, where everything is 0
    timeseries = convert_to_timeseries(saved_state)
    plot_simulation_output(timeseries, plot_settings, out_dir)
    save_timeseries(timeseries, out_dir)


if __name__ == '__main__':
    main()
