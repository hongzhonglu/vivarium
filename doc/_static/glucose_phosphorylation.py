from vivarium.compartment.process import Process
from vivarium.compartment.composition import (
    simulate_process,
    plot_simulation_output,
)


class GlucosePhosphorylation(Process):

    default_parameters = {
        'k_cat': 2e-3,
        'K_ATP': 5e-2,
        'K_GLC': 4e-2,
    }

    def __init__(self, initial_parameters=None):
        ports = {
            'nucleoside_phosphates': ['ATP', 'ADP'],
            'cytoplasm': ['GLC', 'G6P', 'HK'],
        }
        parameters = GlucosePhosphorylation.default_parameters
        parameters.update(initial_parameters)
        super(GlucosePhosphorylation, self).__init__(
            ports, parameters)

    def next_update(self, timestep, states):
        # Get concentrations from state
        cytoplasm = states['cytoplasm']
        nucleoside_phosphates = states['nucleoside_phosphates']
        hk = cytoplasm['HK']
        glc = cytoplasm['GLC']
        atp = nucleoside_phosphates['ATP']

        # Get kinetic parameters
        k_cat = self.parameters['k_cat']
        k_atp = self.parameters['K_ATP']
        k_glc = self.parameters['K_GLC']

        # Compute reaction rate with michaelis-menten equation
        rate = k_cat * hk * glc * atp / (
            k_glc * k_atp + k_glc * atp + k_atp * glc + glc * atp)

        # Compute concentration changes from rate and timestep
        delta_glc = -rate * timestep
        delta_atp = -rate * timestep
        delta_g6p = rate * timestep
        delta_adp = rate * timestep

        # Compile changes into an update
        update = {
            'cytoplasm': {
                'GLC': delta_glc,
                'G6P': delta_g6p,
                # We exclude HK because it doesn't change
            },
            'nucleoside_phosphates': {
                'ATP': delta_atp,
                'ADP': delta_adp,
            },
        }

        return update

    def default_settings(self):
        default_state = {
            'cytoplasm': {
                'GLC': 1.0,
                'G6P': 0.0,
                'HK': 0.1,
            },
            'nucleoside_phosphates': {
                'ATP': 2.0,
                'ADP': 0.0,
            },
        }
        schema = {
            'cytoplasm': {
                'GLC': {
                    # accumulate means to add the updates
                    'updater': 'accumulate',
                },
                # accumulate is the default, so we don't need to specify
                # updaters for the rest of the variables
            },
        }
        emitter_keys = {
            # We want to track the substrates and products, but not HK
            'cytoplasm': ['GLC', 'G6P'],
            'nucleoside_phosphates': ['ATP', 'ADP'],
        }
        return {
            'state': default_state,
            'emitter_keys': emitter_keys,
            'schema': schema,
        }


if __name__ == '__main__':
    parameters = {
        'k_cat': 1.5,
    }
    my_process = GlucosePhosphorylation(parameters)

    settings = {
        'total_time': 10,
        'timestep': 0.1,
    }
    timeseries = simulate_process(my_process, settings)
    plot_simulation_output(timeseries, {}, './')
