'''
boot_compartments

boot functions for compartment simulations. Can be called from the command line
after booting an environmental "outer" simulation and passing in the outer-id, the agent_type:

>>python -m lens.boot_compartment --outer-id 'outer-id' --type 'agent_type'

For booting an environment, see see the top-level README or lens/environment/boot

Once an environment and agent are booted, they can be triggered with:

>>python -m lens.environment.control run --id 'outer-id'

'''

from __future__ import absolute_import, division, print_function

from lens.actor.inner import Inner
from lens.actor.boot import BootAgent
from lens.environment.lattice_compartment import generate_lattice_compartment

# processes
from lens.processes.transport_lookup import TransportLookup
from lens.processes.CovertPalsson2002_metabolism import Metabolism
from lens.processes.CovertPalsson2002_regulation import Regulation
from lens.processes.Kremling2007_transport import Transport
from lens.processes.growth import Growth
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity

# composites
from lens.composites.Covert2008 import initialize_covert2008
from lens.composites.growth_division import initialize_growth_division
from lens.composites.Vladimirov2008 import initialize_vladimirov2008

def wrap_boot(initialize, initial_state):
    def boot(agent_id, agent_type, agent_config):
        initial_state.update(agent_config.get('declare', {}))
        agent_config['declare'] = initial_state  # 'declare' is for the environment

        return Inner(
            agent_id,
            agent_type,
            agent_config,
            initialize)

    return boot

def wrap_initialize(make_process):
    def initialize(boot_config):
        boot_config.update({
            'emitter': {
                'type': 'database',
                'url': 'localhost:27017',
                'database': 'simulations',
                },
            'compartment_options':{
                'time_step': 1.0,
                },
            })
        process = make_process(boot_config)  # 'boot_config', set in environment.control is the process' initial_parameters
        return generate_lattice_compartment(process, boot_config)

    return initialize


class BootCompartment(BootAgent):
    def __init__(self):
        super(BootCompartment, self).__init__()
        self.agent_types = {
            # single process compartments
            'lookup': wrap_boot(wrap_initialize(TransportLookup), {'volume': 1.0}),
            'metabolism': wrap_boot(wrap_initialize(Metabolism), {'volume': 1.0}),
            'regulation': wrap_boot(wrap_initialize(Regulation), {'volume': 1.0}),
            'transport': wrap_boot(wrap_initialize(Transport), {'volume': 1.0}),
            'growth': wrap_boot(wrap_initialize(Growth), {'volume': 1.0}),
            'receptor': wrap_boot(wrap_initialize(ReceptorCluster), {'volume': 1.0}),
            'motor': wrap_boot(wrap_initialize(MotorActivity), {'volume': 1.0}),
            # composite compartments
            'growth_division': wrap_boot(initialize_growth_division, {'volume': 1.0}),
            'covert': wrap_boot(initialize_covert2008, {'volume': 1.0}),
            'vladimirov': wrap_boot(initialize_vladimirov2008, {'volume': 1.0})
        }

def run():
    boot = BootCompartment()
    boot.execute()

if __name__ == '__main__':
    run()
