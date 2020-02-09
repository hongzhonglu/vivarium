from __future__ import absolute_import, division, print_function

import os

from vivarium.actor.process import initialize_state, load_compartment, convert_to_timeseries, plot_simulation_output, \
    simulate_with_environment

# processes
from vivarium.processes.ode_expression import ODE_expression, get_flagella_expression
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Mears2014_flagella_activity import FlagellaActivity
from vivarium.processes.membrane_potential import MembranePotential
from vivarium.processes.convenience_kinetics import ConvenienceKinetics, get_glc_lct_config
from vivarium.processes.metabolism import Metabolism, get_e_coli_core_config
from vivarium.processes.deriver import Deriver
from vivarium.processes.division import Division, divide_condition, divide_state


def compose_pmf_chemotaxis(config):
    receptor_parameters = {'ligand': 'GLC'}
    receptor_parameters.update(config)

    # declare the processes
    transport = ConvenienceKinetics(config.get('transport', get_glc_lct_config()))
    metabolism = Metabolism(config.get('metabolism', get_e_coli_core_config()))
    expression = ODE_expression(config.get('expression', get_flagella_expression()))
    receptor = ReceptorCluster(config.get('receptor', receptor_parameters))
    flagella = FlagellaActivity(config.get('flagella', {}))
    PMF = MembranePotential(config.get('PMF', {}))
    deriver = Deriver(config.get('deriver', {}))
    division = Division(config.get('division', {}))

    # place processes in layers
    processes = [
        {'PMF': PMF},
        {'receptor': receptor,
         'transport': transport
         },
        {
         # 'metabolism': metabolism,
         'expression': expression},
        {'flagella': flagella},
        {'deriver': deriver,
         'division': division}]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cytoplasm'},
        'transport': {
            'exchange': 'exchange',
            'external': 'environment',
            'internal': 'cytoplasm',
            'fluxes': 'fluxes'},
        # 'metabolism': {
        #     'internal': 'cytoplasm',
        #     'external': 'environment',
        #     'reactions': 'reactions',
        #     'exchange': 'exchange',
        #     'flux_bounds': 'fluxes'},
        'expression' : {
            'counts': 'cell_counts',
            'internal': 'cytoplasm',  # todo -- hook this up with flagella
            'external': 'environment'},
        'flagella': {
            'internal': 'cytoplasm',
            'membrane': 'membrane',
            'flagella': 'flagella',
            'external': 'environment'},
        'PMF': {
            'external': 'environment',
            'membrane': 'membrane',
            'internal': 'cytoplasm'},
        'deriver': {
            'counts': 'cell_counts',
            'state': 'cytoplasm',
            'prior_state': 'prior_state'},
        'division': {
            'internal': 'cytoplasm'}}

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'name': 'PMF_chemotaxis_composite',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment_role': 'environment',
        'exchange_role': 'exchange',
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}



if __name__ == '__main__':

    out_dir = os.path.join('out', 'tests', 'PMF_chemotaxis')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    boot_config = {'emitter': 'null'}
    compartment = load_compartment(compose_pmf_chemotaxis, boot_config)

    # settings for simulation and plot
    options = compartment.configuration
    timeline = [(10, {})]

    settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux_bounds'},
        'skip_roles': [
            'prior_state', 'null']}

    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
