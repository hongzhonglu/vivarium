from __future__ import absolute_import, division, print_function

from lens.actor.process import initialize_state

# processes
from lens.processes.derive_volume import DeriveVolume
from lens.processes.growth import Growth
from lens.processes.division import Division, divide_condition, divide_state
from lens.processes.protein_expression import ProteinExpression



def compose_growth_division(config):

    # declare the processes
    growth = Growth(config)
    division = Division(config)
    expression = ProteinExpression(config)
    deriver = DeriveVolume(config)

    # place processes in layers
    processes = [
        {'growth': growth,
         'division': division,
         'expression': expression},
        {'deriver': deriver}]

    # make the topology.
    # for each process, map process roles to compartment roles
    topology = {
        'growth': {
            'internal': 'cell'},
        'division': {
            'internal': 'cell'},
        'expression': {
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    # get environment ids, and make exchange_ids for external state
    # TODO -- remove this
    environment_ids = []

    options = {
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'environment': 'environment',
        'compartment': 'cell',
        'environment_ids': environment_ids,
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}
