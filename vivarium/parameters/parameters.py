from __future__ import absolute_import, division, print_function

import copy

import numpy as np

from vivarium.actor.process import load_compartment, simulate_compartment, simulate_with_environment

# composites
from vivarium.composites.kinetic_FBA import compose_kinetic_FBA


def get_nested(dict, keys):
    d = dict
    for key in keys[:-1]:
        if key in d:
            d = d[key]
    if keys[-1] in d:
        value = d[keys[-1]]
    return value

def set_nested(dict, keys, value, create_missing=True):
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

def parameter_scan(composite, scan_params, output_values):

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
            new_params = set_nested(params, param_id, value)
            param_sets.append(new_params)

    # simulation settings
    total_time = 10
    settings = {
        'timestep': 1,
        'total_time': total_time}

    results = []
    for params in param_sets:
        new_compartment = load_compartment(compose_kinetic_FBA, params)
        sim_out = simulate_compartment(new_compartment, settings)  # TODO -- supress emitter. null emitter?
        last_state = sim_out[total_time]

        output = []
        for output_value in output_values:
            output.append(get_nested(last_state, output_value))
        results.append(output)

    return results

def scan_kinetic_FBA():
    scan_params = {
        ('transport', 'kinetic_parameters', 'EX_glc__D_e', 'PTSG_internal', 'kcat_f'): [-3e4, -3e3, -3e2, -3e1]
    }

    output_values = [
        ('reactions', 'EX_glc__D_e')
    ]

    results = parameter_scan(compose_kinetic_FBA, scan_params, output_values)

    print('results: {}'.format(results))


    import ipdb;
    ipdb.set_trace()




if __name__ == '__main__':
    scan_kinetic_FBA()
