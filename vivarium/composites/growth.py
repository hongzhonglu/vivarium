from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.tree import process_derivers

# processes
from vivarium.processes.growth import Growth


def growth_compartment(config):


    # declare the processes
    growth = Growth(config.get('growth', {}))


    # place processes in layers
    processes = {
        'growth': growth
    }

    global_key = ['..', 'global']

    # make the topology.
    # for each process, map process ports to store ids
    topology = {
        'growth': {
            'global': global_key},
        }

    # add derivers
    derivers = process_derivers(processes, topology)
    processes.update(derivers['processes'])
    topology.update(derivers['topology'])  # add derivers to the topology

    return {
        'processes': processes,
        'topology': topology}