from __future__ import absolute_import, division, print_function

import os
import argparse

from vivarium.core.tree import Compartment
from vivarium.core.composition import compartment_in_experiment




from vivarium.core.composition import (
    simulate_with_environment,
    plot_simulation_output,
    load_compartment)
from vivarium.parameters.parameters import (
    parameter_scan,
    get_parameters_logspace,
    plot_scan_results)

# processes
from vivarium.processes.meta_division import MetaDivision
from vivarium.processes.metabolism import (
    Metabolism,
    get_iAF1260b_config)
from vivarium.processes.convenience_kinetics import (
    ConvenienceKinetics,
    get_glc_lct_config)
from vivarium.processes.ode_expression import (
    ODE_expression,
    get_lacy_config)


def default_metabolism_config():
    config = get_iAF1260b_config()

    # set flux bond tolerance for reactions in ode_expression's lacy_config
    metabolism_config = {
        'moma': False,
        'tolerance': {
            'EX_glc__D_e': [1.05, 1.0],
            'EX_lcts_e': [1.05, 1.0]}}
    config.update(metabolism_config)
    return config


class TransportMetabolismExpression(Compartment):
    """
    TransportMetabolismExpression Compartment
    """

    defaults = {
        'transport': get_glc_lct_config(),
        'metabolism': default_metabolism_config(),
        'expression': get_lacy_config(),
        'division': {}}

    def __init__(self, config):
        self.transport_config = config.get('transport', self.defaults['transport'])
        self.metabolism_config = config.get('metabolism', self.defaults['metabolism'])
        self.expression_config = config.get('expression', self.defaults['expression'])
        self.division_config = config.get('division', self.defaults['division'])

    def generate_processes(self, config):
        # Transport
        # load the kinetic parameters
        transport = ConvenienceKinetics(config.get(
            'transport',
            self.transport_config))

        # Metabolism
        # get target fluxes from transport, and update constrained_reaction_ids
        metabolism_config = config.get(
            'metabolism',
            self.metabolism_config)
        target_fluxes = transport.kinetic_rate_laws.reaction_ids
        metabolism_config.update({'constrained_reaction_ids': target_fluxes})
        metabolism = Metabolism(metabolism_config)

        # Gene expression
        expression = ODE_expression(config.get(
            'expression',
            self.expression_config))

        # Division
        division_config = config.get(
            'division',
            self.division_config)
        # initial_mass = metabolism.initial_mass
        # division_config.update({'constrained_reaction_ids': target_fluxes})
        # TODO -- configure metadivision
        division = MetaDivision(division_config)

        return {
            'transport': transport,
            'metabolism': metabolism,
            'expression': expression,
            'division': division}

    def generate_topology(self, config):
        return {
        'transport': {
            'internal': 'cytoplasm',
            'external': 'environment',
            'exchange': 'null',  # metabolism's exchange is used
            'fluxes': 'flux_bounds',
            'global': 'global'},
        'metabolism': {
            'internal': 'cytoplasm',
            'external': 'environment',
            'reactions': 'reactions',
            'exchange': 'exchange',
            'flux_bounds': 'flux_bounds',
            'global': 'global'},
        'expression': {
            'counts': 'cytoplasm_counts',
            'internal': 'cytoplasm',
            'external': 'environment'},
        'division': {
            'global': 'global'}}


# simulate
def test_txp_mtb_ge(config={}, time=10):

    # configure the compartment
    compartment = TransportMetabolismExpression(config)

    # configure experiment
    experiment_settings = {
        'compartment': config}
    experiment = compartment_in_experiment(
        compartment,
        experiment_settings)

    import ipdb;
    ipdb.set_trace()

    # run experiment
    timestep = 1
    time = 0
    while time < end_time:
        experiment.update(timestep)
        time += timestep
    return experiment.emitter.get_data()

def simulate_txp_mtb_ge(config={}, out_dir='out'):

    # run simulation
    timeseries = test_txp_mtb_ge(config, 2520) # 2520 sec (42 min) is the expected doubling time in minimal media
    volume_ts = timeseries['global']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))

    # plot
    plot_settings = {
        'max_rows': 30,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux_bounds'},
        'skip_ports': [
            'prior_state', 'null'],
        'show_state': [
            ('reactions', 'EX_glc__D_e'),
            ('reactions', 'EX_lcts_e')]}
    plot_simulation_output(timeseries, plot_settings, out_dir)


# parameters
def scan_txp_mtb_ge():
    composite_function = compose_txp_mtb_ge

    # parameters to be scanned, and their values
    scan_params = {
        ('transport',
         'kinetic_parameters',
         'EX_glc__D_e',
         ('internal', 'EIIglc'),
         'kcat_f'):
            get_parameters_logspace(1e3, 1e6, 4),
        ('transport',
         'kinetic_parameters',
         'EX_lcts_e',
         ('internal', 'LacY'),
         'kcat_f'):
            get_parameters_logspace(1e3, 1e6, 4),
    }

    # metrics are the outputs of a scan
    metrics = [
        ('reactions', 'EX_glc__D_e'),
        ('reactions', 'EX_lcts_e'),
        ('global', 'mass')
    ]

    # define conditions
    conditions = [
        # {}, # default
        {
        'environment': {
            'glc__D_e': 12.0,
            'lcts_e': 10.0},
        'cytoplasm':{
            'LacY': 0.0}
        },
        {
        'environment': {
            'glc__D_e': 0.0,
            'lcts_e': 10.0},
        'cytoplasm':{
            'LacY': 1.0e-6}
        },
    ]

    ## TODO -- add targets
    # targets = {
    #     'global', 'growth_rate'
    # }

    # set up scan options
    timeline = [(10, {})]
    sim_settings = {
        'environment_port': 'environment',
        'exchange_port': 'exchange',
        'environment_volume': 1e-6,  # L
        'timeline': timeline}
    scan_options = {
        'simulate_with_environment': True,
        'simulation_settings': sim_settings}

    # run scan
    scan_config = {
        'composite': composite_function,
        'scan_parameters': scan_params,
        'conditions': conditions,
        'metrics': metrics,
        'options': scan_options}
    results = parameter_scan(scan_config)

    return results


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'txp_mtb_ge_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # run scan with python vivarium/compartments/txp_mtb_ge.py --scan
    parser = argparse.ArgumentParser(description='transport metabolism composite')
    parser.add_argument('--scan', '-s', action='store_true', default=False,)
    parser.add_argument('--run', '-r', action='store_true', default=False, )
    args = parser.parse_args()

    if args.scan:
        results = scan_txp_mtb_ge()
        plot_scan_results(results, out_dir)
    else:
        config = {}
        simulate_txp_mtb_ge(config, out_dir)
