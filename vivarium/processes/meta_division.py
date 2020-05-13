from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Process


def divide_condition(compartment):
    division_port = compartment.division_port
    division = compartment.states[division_port].state_for(['division'])
    if division.get('division', 0) == 0:  # 0 means false
        divide = False
    else:
        divide = True
    return divide


class CountForever(object):
    def __init__(self, start=0, by=1):
        self.index = start
        self.by

    def generate(self):
        value = self.index
        self.index += self.by
        return value


class MetaDivision(Process):

    defaults = {
        'initial_state': {},
        'id_function': CountForever().generate}

    def __init__(self, initial_parameters={}):
        self.division = 0

        # must provide a compartment to generate new daughters
        self.compartment = initial_parameters['compartment']
        self.id_function = self.or_default(initial_parameters, 'id_function')
        self.cell_id = initial_parameters.get('cell_id', str(self.id_function()))

        ports = {
            'global': ['division'],
            'cells': ['*']}

        super(Division, self).__init__(ports, parameters)

    def default_settings(self):
        # default state
        default_state = {
            'global': {
                'division': 0}}

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'global': {
                'division': {
                    'updater': 'set',
                    'divide': 'zero'}}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def ports_schema(self):
        default = self.default_settings()
        initial = default['state']

        return {
            'global': {
                'division': {
                    '_default': False,
                    '_updater': 'set'}}
            'cells': {
                '*': {
                    'cell': {}}}}

    def next_update(self, timestep, states):
        division = states['global']['division']

        if division:
            harikari = [self.cell_id]

            # daughter_states = divide_state()

            daughter_ids = [
                self.id_function(), self.id_function()]

            daughter_updates = []
            
            for daughter_id in daughter_ids:
                compartment = self.compartment.generate({
                    'agent_id': daughter_id})
                daughter_updates['cells']['_generate'].append({
                    'path': (daughter_id, 'cell'),
                    'processes': compartment['processes'],
                    'topology': compartment['topology'],
                    # TODO: provide initial state})
                    'initial_state': {}})

            return {
                'cells': {
                    '_delete': harikari,
                    '_generate': daughter_updates}}
                        
