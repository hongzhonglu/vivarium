from vivarium.actor import composition
from vivarium.actor.process import Process


class ToyLinearGrowthDeathProcess(Process):

    GROWTH_RATE = 1.0
    THRESHOLD = 5.0

    def __init__(self, initial_parameters={}):
        roles = {
            'compartment': ['processes'],
            'global': ['mass'],
        }
        super(ToyLinearGrowthDeathProcess, self).__init__(
            roles, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'global': {
                    'mass': 0.0
                }
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        mass = states['global']['mass']
        mass_grown = (
            ToyLinearGrowthDeathProcess.GROWTH_RATE * timestep)
        update = {
            'global': {'mass': mass_grown},
        }
        if mass > ToyLinearGrowthDeathProcess.THRESHOLD:
            update['compartment'] = {
                'processes': [],
            }
        return update


class TestSimulateProcess:

    def test_compartment_state_role(self):
        '''Check that compartment state roles are handled'''
        process = ToyLinearGrowthDeathProcess()
        settings = {
            'compartment_state_role': 'compartment',
        }
        saved_total_states = composition.simulate_process(
            process, settings)
        timeseries = composition.convert_to_timeseries(
            saved_total_states)
        expected_masses = [
            # Mass stops increasing the iteration after mass > 5 because
            # cell dies
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.0, 7.0]
        masses = timeseries['global']['mass']
        assert masses == expected_masses
