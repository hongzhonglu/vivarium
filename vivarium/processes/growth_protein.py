from __future__ import absolute_import, division, print_function

import numpy as np

from vivarium.compartment.process import Process


class GrowthProtein(Process):
    """The Growth :term:`process class` models exponential cell growth.

    The cell's mass :math:`m_{t + h}` at time :math:`t + h` for
    :term:`timestep` :math:`h` and with growth rate :math:`r` is modeled
    as:

    .. math::

        m_{t + h} = m_t e^{rh}

    Configuration Options:

    * ``growth_rate``: The cell's growth rate :math:`r`. This rate is
      0.0006 by default.

      .. todo:: Why is the rate 0.0006?

    Example Usage:

    >>> import math
    >>> TIMESTEP = 1.0  # in seconds
    >>> # growth rate chosen so mass doubles each timestep
    >>> configuration = {'growth_rate': math.log(2.0)}
    >>> growth_process = Growth(configuration)
    >>> state = growth_process.default_settings()['state']
    >>> state
    {'global': {'mass': 1339}}
    >>> update = growth_process.next_update(TIMESTEP, state)
    >>> update
    {'global': {'mass': 2678.0}}
    >>> update['global']['mass'] / state['global']['mass']
    2.0

    """

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

        super(Growth, self).__init__(ports, parameters)

    def ports_schema(self):
        # default state
        # 1000 proteins per fg
        protein = 1339000  # (wet mass in fg)

        return {
            'internal': {
                'protein': {
                    '_default': protein,
                    '_divider': 'split',
                    '_properties': {
                        'mass': 1e-3}}}
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
                'protein': new_protein}
            'global': {
                'divide': divide}}
