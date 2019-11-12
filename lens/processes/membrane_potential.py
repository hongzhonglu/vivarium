from __future__ import absolute_import, division, print_function

import os
import numpy as np

from lens.actor.process import Process
from lens.utils.units import units


class MembranePotential(Process):
    '''
    Need to add a boot method for this process to lens/environment/boot.py for it to run on its own
    '''
    def __init__(self, initial_parameters={}):

        parameters = {
            'R': 8.314462618,  # (J * K^-1 * mol^-1) gas constant
            'F': 96485.33289, # (charge * mol^-1)  Faraday constant
            'p_K': 0.05,  # (unitless, relative permeability) membrane permeability of K
            'p_Na': 0.05,  # (unitless, relative permeability) membrane permeability of Na
            'p_Cl': 0.05,  # (unitless, relative permeability) membrane permeability of Cl
        }

        roles = {
            'internal': ['c_in'],
            'membrane': ['PMF', 'd_V', 'd_pH'],  # proton motive force (PMF), electrical difference (d_V), pH difference (d_pH)
            'external': ['c_out'],
        }
        parameters.update(initial_parameters)

        super(MembranePotential, self).__init__(roles, parameters)

    def default_state(self):
        # TODO -- get real concentrations.
        return {
            'internal': {'K': 30, 'Na': 10, 'Cl': 10},  # (mmol) http://book.bionumbers.org/what-are-the-concentrations-of-different-ions-in-cells/
            'external': {'K': 1, 'Na': 1, 'Cl': 1, 'T': 310.15}  # temperature in Kelvin
        }

    def default_emitter_keys(self):
        keys = {
            'membrane': ['d_V', 'd_pH', 'PMF'],
        }
        return keys

    def default_updaters(self):
        keys = {'membrane': {
            'd_V': 'set',
            'd_pH': 'set',
            'PMF': 'set',
        }}
        return keys

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        # parameters
        R = self.parameters['R']
        F = self.parameters['F']
        p_K = self.parameters['p_K']
        p_Na = self.parameters['p_Na']
        p_Cl = self.parameters['p_Cl']

        # state
        K_in = internal_state['K']
        Na_in = internal_state['Na']
        Cl_in = internal_state['Cl']
        K_out = external_state['K']
        Na_out = external_state['Na']
        Cl_out = external_state['Cl']

        T = external_state['T']  # temperature

        # Membrane potential. Goldman equation
        numerator = p_K * K_out + p_Na * Na_out + p_Cl * Cl_in
        denominator = p_K * K_in + p_Na * Na_in + p_Cl * Cl_out
        d_V = (R * T) / (F) * np.log(numerator / denominator)

        # pH difference
        d_pH = 0

        # proton motive force
        PMF = d_V + d_pH


        update = {
            'membrane': {
                'd_V': d_V,
                'd_pH': d_pH,
                'PMF': PMF,
            }}
        return update

def test_mem_potential():
    timeline = [
        (0, {'external': {
            'Na': 1}
        }),
        (100, {'external': {
            'Na': 2}
        }),
        (500, {}),
    ]

    # configure process
    mp = MembranePotential({})

    # get initial state and parameters
    state = mp.default_state()
    saved_state = {'internal': {}, 'external': {}, 'membrane': {}, 'time': []}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.iteritems():
                    state[key].update(change)

        update = mp.next_update(timestep, state)
        saved_state['time'].append(time)

        # update external state
        for role in ['internal', 'external']:
            for state_id, value in state[role].iteritems():
                if state_id in saved_state[role].keys():
                    saved_state[role][state_id].append(value)
                else:
                    saved_state[role][state_id] = [value]

        # update membrane state from update
        for state_id, value in update['membrane'].iteritems():
            if state_id in saved_state['membrane'].keys():
                saved_state['membrane'][state_id].append(value)
            else:
                saved_state['membrane'][state_id] = [value]

    return saved_state

def plot_mem_potential(saved_state, out_dir='out'):
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

    # plot data
    plot_idx = 0
    for key in data_keys:
        for mol_id, series in sorted(saved_state[key].iteritems()):
            ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)

            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            ax.set_xlabel('time (hrs)')

            if key is 'internal':
                ax.set_yticks([0.0, 1.0])
                ax.set_yticklabels(["False", "True"])

            plot_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'membrane_potential')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_mem_potential()
    out_dir = os.path.join('out', 'membrane_potential')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_mem_potential(saved_state, out_dir)
