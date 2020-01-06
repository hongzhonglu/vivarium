from __future__ import absolute_import, division, print_function

import os

from vivarium.processes.metabolism import Metabolism
from vivarium.environment.make_media import Media
from vivarium.utils.units import units
import vivarium.utils.regulation_logic as rl

DATA_FILE = os.path.join('models', 'e_coli_core.json')
# DATA_FILE = os.path.join('models', 'iAF1260b.json')

def BiGGMetabolism(parameters):
    initial_state = get_initial_state()

    parameters['model_path'] = parameters.get('model_path', DATA_FILE)
    parameters['initial_state'] = parameters.get('initial_state', initial_state)

    return Metabolism(parameters)

def get_initial_state():
    # internal state
    mass = 1339 * units.fg
    density = 1100 * units.g/units.L
    volume = mass.to('g') / density
    internal = {
            'mass': mass.magnitude,  # fg
            'volume': volume.magnitude}

    # external state
    # TODO -- initial state is currently configured to e_coli_core, needs to be generalized
    make_media = Media()
    external = make_media.get_saved_media('ecoli_core_GLC')

    return {
        'internal': internal,
        'external': external}



# test functions
def test_metabolism():
    # configure process
    metabolism = BiGGMetabolism({})

    print('MODEL: {}'.format(metabolism.fba.model))
    print('REACTIONS: {}'.format(metabolism.fba.model.reactions))
    print('METABOLITES: {}'.format(metabolism.fba.model.metabolites))
    print('GENES: {}'.format(metabolism.fba.model.genes))
    print('COMPARTMENTS: {}'.format(metabolism.fba.model.compartments))
    print('SOLVER: {}'.format(metabolism.fba.model.solver))
    print('EXPRESSION: {}'.format(metabolism.fba.model.objective.expression))

    print(metabolism.fba.optimize())
    print(metabolism.fba.model.summary())
    print('internal: {}'.format(metabolism.fba.internal_reactions()))
    print('external: {}'.format(metabolism.fba.external_reactions()))
    print(metabolism.fba.reaction_ids())
    print(metabolism.fba.get_reactions())
    print(metabolism.fba.get_reaction_bounds())
    print(metabolism.fba.read_exchange_fluxes())

def kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux

    return rate



if __name__ == '__main__':
    from vivarium.processes.metabolism import simulate_metabolism, save_network
    from vivarium.actor.process import convert_to_timeseries, plot_simulation_output

    out_dir = os.path.join('out', 'tests', 'BiGG_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # toy functions
    def toy_transport():
        transport_kinetics = {
            # 'GLCptspp': kinetic_rate('glc__D_e', 1.5e0, 5),  # for model iAF1260b
            # 'EX_lac__D_e': kinetic_rate('lac__D_e', -5e-1, 8),
            # 'LACZpp': kinetic_rate('lac__D_e', 1.5e0, 5),
            # 'GLCpts': kinetic_rate('glc__D_e', 1.5e0, 5),  # for model e_coli_core
            # 'LACZ': kinetic_rate('lac__D_e', 1.5e0, 5),
        }
        return transport_kinetics

    def toy_regulation():
        regulation = {
            # 'EX_lac__D_e': rl.build_rule('IF not (glc__D_e_external)'),
            # 'EX_glc__D_e': rl.build_rule('IF (glc__D_e_external)'),
            # 'LACZpp': rl.build_rule('IF not (glc__D_e_external)'),
            # 'LACZ': rl.build_rule('IF not (glc__D_e_external)'),
        }
        return regulation

    # define process config
    config = {'model_path': DATA_FILE}

    # additional process-like transport and regulation functions
    transport = toy_transport()
    regulation = toy_regulation()
    config['constrained_reaction_ids'] = transport.keys()
    config['regulation'] = regulation

    # load metabolism model
    metabolism = BiGGMetabolism(config)

    # simulate model
    timeline = [
        # (0, {'external': {
        #     'glc__D_e': 12.0,
        #     'lac__D_e': 12.0}
        # }),
        # (100, {'external': {
        #     'glc__D_e': 0.0}
        # }),
        (2500, {})]

    simulation_config = {
        'process': metabolism,
        'timeline': timeline,
        'transport_kinetics': toy_transport(),
        'environment_volume': 1e-11}

    plot_settings = {
        'max_rows': 30,
        'remove_flat': True,
        'show_state': [
            ('external', 'glc__D_e'),
            ('internal', 'lcts_p'),
            ('external', 'lac__D_e'),
            ('flux_bounds', 'LACZ'),
            ('reactions', 'LACZ'),
            ('flux_bounds', 'LACZpp'),
            ('reactions', 'LACZpp'),
            ('flux_bounds', 'EX_lac__D_e'),
            ('reactions', 'EX_lac__D_e'),
            ('flux_bounds', 'EX_glc__D_e'),
            ('reactions', 'EX_glc__D_e'),
        ],
        'skip_roles': ['exchange'],
        'overlay': {'reactions': 'flux_bounds'}}

    saved_data = simulate_metabolism(simulation_config)
    del saved_data[0]  # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # make flux network from model
    network_config = {
        'process': metabolism,
        'total_time': 10,
        'transport_kinetics': toy_transport(),
        'environment_volume': 1e-11}
    save_network(network_config, out_dir)
