from __future__ import absolute_import, division, print_function

from vivarium.actor.process import load_compartment, simulate_with_environment

# composites
from vivarium.composites.kinetic_FBA import compose_kinetic_FBA




def scan_kinetic_FBA():
    composite = compose_kinetic_FBA
    timestep = 1

    # a different kcat_f for each item.
    # TODO -- make it so you can scan a param without including all other params
    scan_params = [
        {'transport': {
            'kinetic_parameters': {
                'EX_glc__D_e': {
                    'PTSG_internal': {
                        'glc__D_e_external': 1e-1,
                        'pep_c_internal': None,
                        'kcat_f': -3e4
                    }}}}},
        {'transport': {
            'kinetic_parameters': {
                'EX_glc__D_e': {
                    'PTSG_internal': {
                        'glc__D_e_external': 1e-1,
                        'pep_c_internal': None,
                        'kcat_f': -3e3
                    }}}}},
        {'transport': {
            'kinetic_parameters': {
                'EX_glc__D_e': {
                    'PTSG_internal': {
                        'glc__D_e_external': 1e-1,
                        'pep_c_internal': None,
                        'kcat_f': -3e2
                    }}}}},
        {'transport': {
            'kinetic_parameters': {
                'EX_glc__D_e': {
                    'PTSG_internal': {
                        'glc__D_e_external': 1e-1,
                        'pep_c_internal': None,
                        'kcat_f': -3e1
                    }}}}},
        ]

    # set up environment
    options = compose_kinetic_FBA({})['options']
    timeline = [(10, {})]
    environment_settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    for new_params in scan_params:

        new_compartment = load_compartment(compose_kinetic_FBA, new_params)
        new_params = new_compartment.current_parameters()

        for process, params in new_params.items():
            print('{}: {}'.format(process, params))


        sim_out = simulate_with_environment(new_compartment, environment_settings)

    import ipdb; ipdb.set_trace()




if __name__ == '__main__':
    scan_kinetic_FBA()
