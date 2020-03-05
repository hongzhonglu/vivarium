'''Simulates antibiotic import'''


from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.composition import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_process_with_environment,
)
from vivarium.processes.convenience_kinetics import ConvenienceKinetics


#: Default initial concentrations
DEFAULT_INITIAL_STATE = {
    'internal': {
        'porin': 1.0,  # Membrane pore through which antibiotics enter
        # EcoCyc ID: TRANS-CPLX-201
        'AcrAB-TolC': 1.0,  # Efflux pump complex
        'antibiotic': 0.0,
    },
    'external': {
        'antibiotic': 1.0,
    },
    'exchange': {
        'antibiotic': 0.0,
    },
}

#: Default initial flux levels
DEFAULT_INITIAL_FLUXES = {
    'antibiotic_import': 0.0,
    'antibiotic_export': 0.0,
}


class Antibiotics(ConvenienceKinetics):
    def __init__(self, initial_parameters={}):
        if 'initial_state' not in initial_parameters:
            initial_state = DEFAULT_INITIAL_STATE
        else:
            initial_state = initial_parameters['initial_state']
        if 'fluxes' not in initial_state:
            initial_state['fluxes'] = DEFAULT_INITIAL_FLUXES
        porin_kcat = initial_parameters.get('porin_kcat', 1e-4)
        porin_km = initial_parameters.get('porin_km', 0.6)
        acrABTolC_kcat = initial_parameters.get('acrABTolC_kcat', 2e-4)
        acrABTolC_km = initial_parameters.get('acrABTolC_km', 0.6)
        parameters = {
            'reactions': {
                'antibiotic_import': {
                    'stoichiometry': {
                        ('internal', 'antibiotic'): 1,
                        ('external', 'antibiotic'): -1,
                    },
                    'is_reversible': False,
                    'catalyzed by': [('internal', 'porin')],
                },
                'antibiotic_export': {
                    'stoichiometry': {
                        ('internal', 'antibiotic'): -1,
                        ('external', 'antibiotic'): 1,
                    },
                    'is_reversible': False,
                    'catalyzed by': [
                        ('internal', 'AcrAB-TolC')],
                },
            },
            'kinetic_parameters': {
                'antibiotic_import': {
                    ('internal', 'porin'): {
                        'kcat_f': porin_kcat,
                        ('external', 'antibiotic'): porin_km,
                    },
                },
                'antibiotic_export': {
                    ('internal', 'AcrAB-TolC'): {
                        'kcat_f': acrABTolC_kcat,
                        ('internal', 'antibiotic'): acrABTolC_km,
                    },
                },
            },
            'initial_state': initial_state,
            'ports': {
                'internal': [
                    'porin',
                    'AcrAB-TolC',
                    'antibiotic',
                ],
                'external': ['antibiotic'],
            },
        }

        super(Antibiotics, self).__init__(parameters)

    def default_settings(self):
        default_settings = super(Antibiotics, self).default_settings()
        default_settings.update(
            {
                'process_id': 'antibiotics',
                'emitter_keys': {
                    'internal': ['antibiotic'],
                    'external': ['antibiotic'],
                },
                'updaters': {
                    'exchange': {
                        'antibiotic': 'set',
                    },
                    'fluxes': {
                        'antibiotic_import': 'set',
                        'antibiotic_export': 'set',
                    },
                },
            },
        )
        return default_settings


def test_antibiotics():
    process = Antibiotics()
    settings = {
        'total_time': 4000,
        'exchange_port': 'exchange',
        'environment_port': 'external',
        'environment_volume': 1e-15,  # Units of L
    }
    return simulate_process_with_environment(process, settings)


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'antibiotics')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    saved_data = test_antibiotics()
    del saved_data[0]
    timeseries = convert_to_timeseries(saved_data)
    plot_settings = {}
    plot_simulation_output(timeseries, plot_settings, out_dir)
