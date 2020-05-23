from __future__ import absolute_import, division, print_function

import numpy as np

from vivarium.compartment.process import Process


class GrowthProtein(Process):
 
    defaults = {
        'growth_rate': 0.0006,
        'global_deriver_key': 'global_deriver',
    }

    def __init__(self, initial_parameters={}):
        ports = {
            'internal': [
                'protein'],
            'global': [
                'volume',
                'divide']}

        self.growth_rate = self.or_default(initial_parameters, 'growth_rate')
        self.global_deriver_key = self.or_default(
            initial_parameters, 'global_deriver_key')

        parameters = {
            'growth_rate': self.growth_rate}
        parameters.update(initial_parameters)

        super(GrowthProtein, self).__init__(ports, parameters)

    def ports_schema(self):
        # default state
        # 1000 proteins per fg
        protein = 1339000  # (wet mass in fg)

        return {
            'internal': {
                'protein': {
                    '_default': protein,
                    '_divider': 'split',
                    '_emit': True,
                    '_properties': {
                        'mass': 1e-3}}},
            'global': {
                'volume': {
                    '_updater': 'set',
                    '_divider': 'split'},
                'divide': {
                    '_default': False,
                    '_updater': 'set'}}}

    def derivers(self):
        return {
            self.global_deriver_key: {
                'deriver': 'globals',
                'port_mapping': {
                    'global': 'global'},
                'config': {
                    'width': 1.11}}}

    def next_update(self, timestep, states):
        protein = states['internal']['protein']
        total_protein = protein * np.exp(self.parameters['growth_rate'] * timestep)
        new_protein = int(total_protein - protein)
        extra = total_protein - int(total_protein)

        # simulate remainder
        if np.random.random() < extra:
            new_protein += 1

        divide = np.random.random() < 0.02

        return {
            'internal': {
                'protein': new_protein},
            'global': {
                'divide': divide}}
