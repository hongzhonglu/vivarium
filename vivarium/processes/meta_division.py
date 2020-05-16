from __future__ import absolute_import, division, print_function

from vivarium.compartment.process import Deriver


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
        self.by = by

    def generate(self):
        value = self.index
        self.index += self.by
        return value


class MetaDivision(Deriver):

    defaults = {
        'initial_state': {},
        'id_function': CountForever(start=33).generate}

    def __init__(self, initial_parameters={}):
        self.division = 0

        # must provide a compartment to generate new daughters
        self.compartment = initial_parameters['compartment']
        self.id_function = self.or_default(initial_parameters, 'id_function')
        self.cell_id = initial_parameters.get('cell_id', str(self.id_function()))

        ports = {
            'global': ['divide'],
            'cells': ['*']}

        super(MetaDivision, self).__init__(ports, initial_parameters)

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
                'divide': {
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
                'divide': {
                    '_default': False,
                    '_updater': 'set'}},
            'cells': {
                '*': {}}}

    def next_update(self, timestep, states):
        divide = states['global']['divide']
        cells = states['cells']

        if divide:
            daughter_ids = [
                str(self.id_function()),
                str(self.id_function())]

            daughter_updates = []
            
            for daughter_id in daughter_ids:
                compartment = self.compartment.generate({
                    'agent_id': daughter_id})
                daughter_updates.append({
                    'daughter': daughter_id,
                    'path': (daughter_id, 'cell'),
                    'processes': compartment['processes'],
                    'topology': compartment['topology'],
                    'initial_state': {}})

            # initial state will be provided by division in the tree
            return {
                'cells': {
                    '_divide': {
                        'mother': self.cell_id,
                        'daughters': daughter_updates}}}
        else:
             return {}   
                        
