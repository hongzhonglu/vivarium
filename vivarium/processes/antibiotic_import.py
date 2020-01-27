'''Simulates antibiotic import'''


from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.actor.process import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_process_with_environment,
)
from vivarium.processes.convenience_kinetics import ConvenienceKinetics
from vivarium.utils.units import units


#: Default initial concentrations
DEFAULT_INITIAL_STATE = {
    'internal': {
        'antibiotic_importer': 1.0,  # Transporter
        'antibiotic': 0.0,
    },
    'external': {
        'antibiotic': 1.0,
    },
}

#: Default initial flux levels
DEFAULT_INITIAL_FLUXES = {
    'anitbiotic_import': 0.0,
}


class AntibioticImport(ConvenienceKinetics):
    def __init__(self, initial_parameters={}):
        if 'initial_state' not in initial_parameters:
            initial_state = DEFAULT_INITIAL_STATE
        else:
            initial_state = initial_parameters['initial_state']
        if 'fluxes' not in initial_state:
            initial_state['fluxes'] = DEFAULT_INITIAL_FLUXES
        parameters = {
            'reactions': {
                'antibiotic_import': {
                    'stoichiometry': {
                        ('internal', 'antibiotic'): 1,
                        ('external', 'antibiotic'): -1,
                    },
                    'is_reversible': False,
                    'catalyzed by': [('internal', 'antibiotic_importer')],
                },
            },
            'kinetic_parameters': {
                'antibiotic_import': {
                    ('internal', 'antibiotic_importer'): {
                        'kcat_f': 0.6,
                        ('external', 'antibiotic'): 0.6,
                    }
                }
            },
            'initial_state': initial_state,
            'roles': {
                'internal': ['antibiotic_importer'],
                'external': ['antibiotic'],
            },
        }

        super(AntibioticImport, self).__init__(parameters)

    def default_settings(self):
        default_settings = super(AntibioticImport, self).default_settings()
        default_settings.update(
            {
                'process_id': 'antibiotic_import',
                'emitter_keys': {
                    'internal': ['antibiotic'],
                    'external': ['antibiotic']
                },
            }
        )
        return default_settings


def test_antibiotic_import():
    process = AntibioticImport()
    settings = {
        'iterations': 20,
        'exchange_role': 'exchange',
        'environment_role': 'external',
        'environment_volume': 1e-15, # Units of L

    }
    return simulate_process_with_environment(process, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'antibiotic_import')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    saved_data = test_antibiotic_import()
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_settings = {}
    plot_simulation_output(timeseries, plot_settings, out_dir)
