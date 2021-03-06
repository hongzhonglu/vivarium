from __future__ import absolute_import, division, print_function

import numpy as np

from vivarium.compartment.process import Process


class Growth(Process):
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
        'growth_rate': 0.0006}

    def __init__(self, initial_parameters={}):
        ports = {
            'global': ['mass', 'volume']}

        growth_rate = initial_parameters.get('growth_rate', self.defaults['growth_rate'])

        parameters = {
            'growth_rate': growth_rate}
        parameters.update(initial_parameters)

        super(Growth, self).__init__(ports, parameters)

    def default_settings(self):
        # default state
        mass = 1339  # (wet mass in fg)
        internal = {'mass': mass}
        default_state = {'global': internal}

        # default emitter keys
        default_emitter_keys = {'global': ['mass']}

        # schema
        schema = {
            'global': {
                'mass': {
                    'updater': 'set'}}}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        mass = states['global']['mass']
        new_mass = mass * np.exp(self.parameters['growth_rate'] * timestep)
        return {'global': {'mass': new_mass}}
