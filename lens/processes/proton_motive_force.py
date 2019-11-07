from __future__ import absolute_import, division, print_function

import os

from lens.actor.process import Process
from lens.utils.units import units


class ProtonMotiveForce(Process):
    '''
    Need to add a boot method for this process to lens/environment/boot.py for it to run on its own
    '''
    def __init__(self, initial_parameters={}):

        parameters = {
            'F': 96485.33289, # (charge / mol)  Faraday constant
        }


        roles = {
            'internal': ['Q_in'],
            'membrane': ['V_m'],
            'external': [],
        }
        parameters.update(initial_parameters)

        super(ProtonMotiveForce, self).__init__(roles, parameters)

    def default_state(self):
        '''
        returns dictionary with:
        default_state = {
            'external': states (dict) -- external states ids with default initial values
            'internal': states (dict) -- internal states ids with default initial values
        '''

        return {
            'internal': {},
            'membrane': {},
            'external': {},
        }

    def default_emitter_keys(self):
        '''
        returns dictionary with:
        keys = {
            'internal': states (list), # a list of states to emit from internal
            'external': states (list), # a list of states to emit from external
        }
        '''
        keys = {
            'internal': [],
            'membrane': [],
            'external': [],
        }
        return keys

    def default_updaters(self):
        '''
        define the updater type for each state in roles.
        The default updater is to pass a delta,
        which is accumulated and passed to the environment at every exchange step'''
        keys = {
            'membrane': {
                'V_m': 'accumulate'}}
        return keys

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        Q_in = internal_state['Q_in']

        # V_m = electric potential difference
        # Q_in is the intracellular charge (in mole)
        # C is the membrane capacitance
        # F is the Faraday constant
        V_m = self.parameters['F'] * Q_in / C


        update = {
            # 'internal': {},
            'membrane': {'V_m': V_m},
            # 'external': {},
        }
        return update

def test_PMF():
    saved_state = {}
    timeline = [
        # (0, {'external': {
        #     'GLC': 1}
        # }),
        (500, {}),
    ]

    # configure process
    PMF = ProtonMotiveForce({})

    # get initial state and parameters
    state = PMF.default_state()
    saved_state = {'internal': {}, 'external': {}, 'time': []}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.iteritems():
                    state[key].update(change)

        update = PMF.next_update(timestep, state)
        saved_state['time'].append(time)

        # update external state
        for state_id, value in state['external'].iteritems():
            if state_id in saved_state['external'].keys():
                saved_state['external'][state_id].append(value)
            else:
                saved_state['external'][state_id] = [value]

        # update internal state from update
        for state_id, value in update['internal'].iteritems():
            if state_id in saved_state['internal'].keys():
                saved_state['internal'][state_id].append(value)
            else:
                saved_state['internal'][state_id] = [value]


    return saved_state

def plot_PMF(saved_state, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [float(t) / 3600 for t in saved_state['time']]  # convert to hours

    # make figure, with grid for subplots
    n_data = [len(saved_state[key].keys()) for key in data_keys]
    n_rows = sum(n_data)
    fig = plt.figure(figsize=(8, n_rows * 2.5))
    grid = plt.GridSpec(n_rows + 1, 1, wspace=0.4, hspace=1.5)

    # plot
    # TODO

    # save figure
    fig_path = os.path.join(out_dir, 'PMF')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_PMF()
    out_dir = os.path.join('out', 'proton_motive_force')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_PMF(saved_state, out_dir)
