from __future__ import absolute_import, division, print_function

import os

from lens.processes.metabolism import Metabolism

DATA_FILE = os.path.join('models', 'e_coli_core.json')

def BiGGMetabolism(parameters):
    parameters['model_path'] = parameters.get('model_path', DATA_FILE)

    return Metabolism(parameters)

def kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def toy_transport_kinetics():
    transport_kinetics = {
        "GLCpts": kinetic_rate('glc__D_e', 1e1, 5),  # glucose
        "GLUt2r": kinetic_rate('glu__L_e', 1e1, 5),  # glucose
        # "PYRt2": kinetic_rate('pyr_e', 1e2, 5),
    }
    return transport_kinetics

def test_ecoli_core():
    # configure process
    metabolism = EcoliCoreMetabolism({})

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


if __name__ == '__main__':
    from lens.processes.metabolism import simulate_metabolism, plot_output, save_network

    out_dir = os.path.join('out', 'tests', 'BiGG_metabolism')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # get ecoli core metabolism model
    toy_config = {}
    toy_transport = toy_transport_kinetics()
    toy_config['constrained_reactions'] = toy_transport.keys()

    ecoli_core_metabolism = BiGGMetabolism(toy_config)

    # simulate model
    simulation_config = {
        'process': ecoli_core_metabolism,
        'total_time': 1000,
        'transport_kinetics': toy_transport_kinetics(),
        'environment_volume': 5e-15}
    saved_data = simulate_metabolism(simulation_config)
    plot_output(saved_data, out_dir)

    # make flux network from model
    save_network(ecoli_core_metabolism, 10, out_dir)
