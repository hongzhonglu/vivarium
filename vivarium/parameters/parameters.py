from __future__ import absolute_import, division, print_function

import copy

import numpy as np

from vivarium.actor.process import load_compartment, simulate_with_environment

# composites
from vivarium.composites.kinetic_FBA import compose_kinetic_FBA


def nested_set(dict, keys, value, create_missing=True):
    d = dict
    for key in keys[:-1]:
        if key in d:
            d = d[key]
        elif create_missing:
            d = d.setdefault(key, {})
        else:
            return dict
    if keys[-1] in d or create_missing:
        d[keys[-1]] = value
    return dict


def parameter_scan(composite, scan_params):

    n_values = [len(v) for v in scan_params.values()]
    n_combinations = np.prod(np.array(n_values))

    # get initial parameters
    compartment = load_compartment(composite)
    initial_params = compartment.current_parameters()

    # make list of parameters to place in compartment
    # TODO -- need to make all combinations, not just 1 parameter at a time
    param_sets = []
    for param_id, values in scan_params.items():
        for value in values:
            params = copy.deepcopy(initial_params)
            new_params = nested_set(params, param_id, value)
            param_sets.append(new_params)

    # set up environment
    options = composite({})['options']
    timeline = [(2, {})]
    environment_settings = {
        'environment_role': options['environment_role'],
        'exchange_role': options['exchange_role'],
        'environment_volume': 1e-13,  # L
        'timeline': timeline}

    for params in param_sets:
        new_compartment = load_compartment(compose_kinetic_FBA, params)

        # new_params = new_compartment.current_parameters()
        # for process, params in new_params.items():
        #     print('{}: {}'.format(process, params))

        # TODO -- supress emitter. null emitter?
        sim_out = simulate_with_environment(new_compartment, environment_settings)

        # TODO -- what is the desired output?

        import ipdb;
        ipdb.set_trace()


def scan_kinetic_FBA():
    scan_params = {
        ('transport', 'kinetic_parameters', 'EX_glc__D_e', 'PTSG_internal', 'kcat_f'): [-3e4, -3e3, -3e2, -3e1]
    }

    parameter_scan(compose_kinetic_FBA, scan_params)



if __name__ == '__main__':
    scan_kinetic_FBA()
