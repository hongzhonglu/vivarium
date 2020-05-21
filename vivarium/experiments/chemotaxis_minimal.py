from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.tree import (
    generate_state,
    Experiment)
from vivarium.compartment.composition import (
    make_agents,
    simulate_compartment,
    load_compartment,
    get_derivers
)

# compartments
from vivarium.composites.lattice_environment import Lattice
from vivarium.composites.chemotaxis_minimal import ChemotaxisMinimal



def chemotaxis_experiment(config={}):

    # get the environment
    environment = Lattice(config.get('environment', {}))
    processes = environment['processes']
    topology = environment['topology']

    # get the agents
    chemotaxis = ChemotaxisMinimal({
        'cells_key': ('..', 'agents')})
    agents = make_agents(range(count), chemotaxis, {})
    processes['agents'] = agents['processes']
    topology['agents'] = agents['topology']

    experiment = Experiment({
        'processes': processes,
        'topology': topology,
        'initial_state': config.get('initial_state', {})})

    import ipdb; ipdb.set_trace()

    return {}



if __name__ == '__main__':
    out_dir = os.path.join('out', 'experiments', 'minimal_chemotaxis')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    data = chemotaxis_experiment()
