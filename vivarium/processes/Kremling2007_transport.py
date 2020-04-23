from __future__ import absolute_import, division, print_function

import os
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from vivarium.compartment.process import Process
from vivarium.utils.flux_conversion import millimolar_to_counts, counts_to_millimolar
from vivarium.environment.make_media import Media

DEFAULT_PARAMETERS = {
    # enzyme synthesis
    'k1': 0.00001,  # enzyme synthesis k1 [\mu mo / gDW h]
    'k2': 0.0001,  # enzyme synthesis k2 [\mu mo / gDW h]
    'k3': 0.00016,  # enzyme synthesis k3 [\mu mo / gDW h]
    'K1': 3000,  # K1 [\mu mo / gDW]
    'K2': 2800,  # K2 [\mu mo / gDW]
    'K3': 15000,  # K3 [\mu mo / gDW]

    # enzyme degration
    'kd': 0.4,  # enzyme degradation k_d[1 / h]

    'm': 1,  # m
    'n': 2,  # n
    'x0': 0.1,  # EIIA_0 [\mu mo / gDW]

    'kg6p': 2.8e6,  # Glc 6P transporter k_g6p[1 / h]
    'Kg6p': 0.1,  # Glc 6P transporter K_g6p[g / l]

    'kptsup': 2.7e8,  # Glc transporter k_pts_up[1 / h]
    'Kglc': 0.12,  # Glc transporter K_glc[g / l]
    'Keiiap': 12,  # Glc transporter K_EIIAP[\mu mol / gDW]

    'klac': 5.4e5,  # Lac transporter k_lac [1 / h]
    'Km_lac': 0.13,  # Km for lac  [g / l]
    'Kieiia': 5.0,  # K_IEIIA: inhibition of Lac transporter by EIIA[-]

    'kgly': 2.80e4,  # [1 / h]
    'kpyk': 9.39e5,  # [1 / (\mu mol / gDW) ^ 3 h]
    'kpdh': 5.50e3,  # [1 / h]
    'kpts': 1.86e5,  # [1(\mu mol / gDW) g]

    # 'K_pts': 0.7,  # K_pts
    'km_pts': 0.7 * 1.86e5,  # K_pts * kpts

    'Y': 1.0e-4,  # Yglc [g TM / mu mol Substrate]

    'mw1': 2.602e-4,  # molecular weight g6p [g Substrate / mu mol Substrate]
    'mw2': 1.802e-4,  # molecular weight glc [g Substrate / mu mol Substrate]
    'mw3': 3.423e-4,  # molecular weight lac [g Substrate / mu mol Substrate]

    'Y1_sim': 6.2448e-05,  # Yg6p [gDW /\mu mol]
    'Y2_sim': 1.0e-4,  # Yglc [gDW /\mu mol]
    'Y3_sim': 9.2421e-05,  # Ylac [gDW /\mu mol]
    'Y4_sim': 1.0e-04,  # Mtl[gDW /\mu mol]

    'K': 0.4,  # K[\mu mol / gDW]
    'kb': 600,  # k_bias[-]
    'ksyn': 3.2623e3,  # k_syn[-] rate of protein synthesis

    'KI': 1 / 8000,  # K_I: Glc 6P transporter inhibits Glc transp.syn.[\mu mol / gDW]
}

MOLECULAR_WEIGHTS = {
    'ACET': 60.05,
    'CO+2': 44.0095,
    'ETOH': 46.06844,
    'FORMATE': 46.0254,
    'GLYCEROL': 92.09382,
    'LAC': 90.08,
    'LCTS': 342.3,
    'OXYGEN-MOLECULE': 31.9988,
    'PI': 94.973,
    'PYR': 88.06,
    'RIB': 150.13,
    'SUC': 118.09,
    'G6P': 260.136,
    'GLC': 180.16,
}



class Transport(Process):

    defaults = {
        'target_fluxes': [],  # ['glc__D_e', 'GLCpts', 'PPS', 'PYK']
        'parameters': DEFAULT_PARAMETERS
    }


    def __init__(self, initial_parameters={}):
        self.dt = 0.01  # timestep for ode integration (seconds)
        self.target_fluxes = initial_parameters.get('target_fluxes', self.defaults['target_fluxes'])

        default_settings = self.default_settings()
        default_state = default_settings['state']
        internal_state = default_state['internal']
        external_state = default_state['external']

        ports = {
            'external': list(external_state.keys()),
            'exchange': list(external_state.keys()),
            'internal': list(internal_state.keys()),
            'fluxes': self.target_fluxes,
            'global': ['volume']}

        parameters = self.defaults['parameters']
        parameters.update(initial_parameters)

        super(Transport, self).__init__(ports, parameters)

    def default_settings(self):

        # default state
        # TODO -- select state based on media
        glc_g6p = True
        glc_lct = False

        make_media = Media()

        # TODO -- don't use this if/else, select state based on media
        if glc_g6p:
            external = make_media.get_saved_media('GLC_G6P')
            internal = {
                'mass': 0.032,  # [g / l]
                'LACZ': 0.0,  # absent in GLC_G6P condition
                'UHPT': 0.0003,  # [\mu mol gDCW] enz g6p
                'PTSG': 0.007,  # [\mu mol gDCW] enz glc
                'G6P': 0.2057,  # [\mu mol gDCW]
                'PEP': 2.0949,  # [\mu mol gDCW]
                'PYR': 2.0949,  # [\mu mol gDCW]
                'XP': 0.0038,  # [fraction phosphorylation]
            }
        elif glc_lct:
            external = make_media.get_saved_media('GLC_LCT')
            internal = {
                'mass': 0.032,  # [g / l]
                'LACZ': 1e-5,  # [\mu mol gDCW] enz g6p
                'UHPT': 0.0,  # absent in GLC_G6P condition
                'PTSG': 0.001,  # [\mu mol gDCW] enz glc
                'G6P': 0.1,  # [\mu mol gDCW]
                'PEP': 0.05,  # [\mu mol gDCW]
                'PYR': 0.1,  # [\mu mol gDCW]
                'XP': 0.01,  # [fraction phosphorylation]
            }
        self.environment_ids = list(external.keys())

        default_state = {
            'internal': internal,
            'external': external,
            'exchange': {state_id: 0.0 for state_id in self.environment_ids},
            'fluxes': {},
            'global': {'volume': 1}}

        # default emitter keys
        default_emitter_keys = {
            'internal': ['mass', 'UHPT', 'PTSG', 'G6P', 'PEP', 'PYR', 'XP'],
            'external': ['G6P', 'GLC', 'LAC'],
            'fluxes': self.target_fluxes}

        # schema
        set_internal = ['mass', 'UHPT', 'LACZ', 'PTSG', 'G6P', 'PEP', 'PYR', 'XP']
        internal_schema = {
            state_id: {
                'updater': 'set',
                'divide': 'set'}
            for state_id in set_internal}
        fluxes_schema = {
            state_id: {
                'updater': 'set',
                'divide': 'set'}
            for state_id in self.target_fluxes}
        exchange_schema = {
            mol_id: {
                'updater': 'accumulate'}
            for mol_id in self.environment_ids}

        schema = {
            'internal': internal_schema,
            'fluxes': fluxes_schema,
            'exchange': exchange_schema}

        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):

        def model(state, t):
            '''
            Model of sugar transport based on Kremling, Bettenbrock, & Gilles. (2007).
            Analysis of global control of Escherichia coli carbohydrate uptake
            '''

            biomass = state[state_keys.index('mass')]  # mass

            # transporters
            UHPT = state[state_keys.index('UHPT')]
            LACZ = state[state_keys.index('LACZ')]
            PTSG = state[state_keys.index('PTSG')]  # transporter Glc (transporter2)

            # metabolites
            G6P = state[state_keys.index('G6P')]  # Glc 6P
            PEP = state[state_keys.index('PEP')]  # Pep
            PYR = state[state_keys.index('PYR')]  # Pyruvate
            XP = state[state_keys.index('XP')]  # EIIAP (phosphorylated PTS Protein)

            # external sugars
            GLC_e = state[state_keys.index('GLC[e]')]
            G6P_e = state[state_keys.index('G6P[e]')]
            LCTS_e = state[state_keys.index('LCTS[e]')]

            # TODO --- G6P_e?
            G6P_e_present = G6P > 0.01
            if G6P_e_present:
                sugar1 = G6P_e
                transporter1 = UHPT
            else:
                sugar1 = LCTS_e
                transporter1 = LACZ

            ## uptake rates
            # sugar1 uptake
            if G6P_e_present:
                # G6P uptake. eqn 40
                uptake1 = p['kg6p'] * (transporter1 * sugar1) / (p['Kg6p'] + sugar1)
            else:
                # Lactose uptake. eqn 45
                uptake1 = p['klac'] * (transporter1 * sugar1) / (
                            p['Km_lac'] + sugar1 * (1 + ((p['x0'] - XP) / p['x0']) / p['Kieiia']))

            # PTS uptake. eqn 38
            uptake2 = p['kptsup'] * XP * (PTSG * GLC_e) / (
                    p['Kglc'] * p['Keiiap'] * p['x0'] + GLC_e * p['Keiiap'] * p['x0'] + XP * p['Kglc'] + XP * GLC_e)

            # enzyme synthesis. signmoid eqn 37
            # Hill coefficient is high (n=6), indicating narrow range of input
            if G6P_e_present:
                # UHPT synthesis. eqn 42
                synthesis1 = p['k1'] * (p['kb'] + p['ksyn'] * XP ** 6 / (XP ** 6 + p['K'] ** 6)) * uptake1 / (
                            p['K1'] + uptake1)

                # PTSG synthesis. eqn 43
                synthesis2 = p['k2'] * (p['KI'] / (transporter1 + p['KI'])) * (
                            p['kb'] + p['ksyn'] * XP ** 6 / (XP ** 6 + p['K'] ** 6)) * uptake2 / (p['K2'] + uptake2)
            else:
                synthesis1 = p['k3'] * (p['kb'] + p['ksyn'] * XP ** 6 / (XP ** 6 + p['K'] ** 6)) * uptake1 / (
                            p['K3'] + uptake1)
                synthesis2 = p['k2'] * (p['kb'] + p['ksyn'] * XP ** 6 / (XP ** 6 + p['K'] ** 6)) * uptake2 / (
                            p['K2'] + uptake2)

            # rates
            rgly = p['kgly'] * G6P  # Glycolyse. eqn 10
            rpdh = p['kpdh'] * PYR  # pdh. eqn 11
            rpts = p['kpts'] * PEP * (p['x0'] - XP) - p['km_pts'] * PYR * XP  # PTS rate. eqn 12
            f = (G6P ** p['n']) * PEP ** p[
                'm']  # represent different model variants.  TODO -- are these variants modeled?
            rpyk = p['kpyk'] * PEP * f  # pyk. eqn 13

            # growth rate. eqn 44
            if G6P_e_present:
                mu = p['Y1_sim'] * uptake1 + p['Y2_sim'] * uptake2
            else:
                mu = p['Y3_sim'] * uptake1 + p['Y2_sim'] * uptake2

            # iFBA code for modifying dPEP
            # TODO -- implement or remove code below
            # if (nargin >= 7 & & ~isempty(varargin{1}) & & ~isempty(varargin{2})):
            #     ppc_rate = (FBA_primal(FBA_rxnInd.PPC) - FBA_primal(FBA_rxnInd.PCKA)) * 1e3;
            #     if (nargin == 7 | | (nargin == 8 & & varargin{4}~=0))
            #         mu = options.iFBA_growRateScale * FBA_primal(FBA_rxnInd.VGRO);
            # else:
            #     ppc_rate = 0;

            # get derivatives
            dbiomass = mu * biomass  # mass

            # sugar uptake and transporter synthesis
            if G6P_e_present:
                dLCTS_e = 0.0
                dLACZ = 0.0
                dG6P_e = -p['mw1'] * uptake1 * biomass
                dUHPT = synthesis1 - (p['kd'] + mu) * transporter1  # transporter synthesis. eqn 46
            else:
                dLCTS_e = -p['mw3'] * uptake1 * biomass
                dLACZ = synthesis1 - (p['kd'] + mu) * transporter1  # transporter synthesis. eqn 46
                dG6P_e = 0.0
                dUHPT = 0.0

            dGLC_e = -p['mw2'] * uptake2 * biomass  # GLC[e]
            dPTSG = synthesis2 - (p['kd'] + mu) * PTSG  # transporter synthesis. eqn 46

            # metabolism
            dG6P = uptake1 + uptake2 - rgly  # eqn 2
            dPEP = 2 * rgly - rpyk - rpts  # - ppc_rate  # eqn 3
            dPYR = rpyk + rpts - rpdh  # eqn 4
            dXP = rpts - uptake2  # eqn 5

            # save to numpy array
            update = np.zeros_like(state)
            update[state_keys.index('mass')] = dbiomass
            update[state_keys.index('UHPT')] = dUHPT  # transporter1 changes with condition
            update[state_keys.index('LACZ')] = dLACZ  # transporter1 changes with condition
            update[state_keys.index('PTSG')] = dPTSG
            update[state_keys.index('G6P')] = dG6P
            update[state_keys.index('PEP')] = dPEP
            update[state_keys.index('PYR')] = dPYR
            update[state_keys.index('XP')] = dXP
            update[state_keys.index('GLC[e]')] = dGLC_e
            update[state_keys.index('G6P[e]')] = dG6P_e    # sugar1 changes with condition
            update[state_keys.index('LCTS[e]')] = dLCTS_e  # sugar1 changes with condition

            # flux states, for metabolism constraint
            # TODO: rup2 is PTS uptake, rup1 is G6P uptake or Lactose uptake,
            update[state_keys.index('GLCpts')] = uptake2
            update[state_keys.index('PPS')] = uptake2  # v_pts used for both PTS uptake and PEP->PYR reaction
            update[state_keys.index('PYK')] = rpyk
            update[state_keys.index('glc__D_e')] = dG6P_e

            return update

        # set up state and parameters for odeint
        timestep_hours = timestep / 3600
        dt_hours = self.dt / 3600
        p = self.parameters
        t = np.arange(0, timestep_hours, dt_hours)

        # get states
        volume = states['global']['volume'] * 1e-15  # convert volume fL to L
        combined_state = {
            'mass': states['internal']['mass'],  # mass
            'UHPT': states['internal']['UHPT'],
            'LACZ': states['internal']['LACZ'],
            'PTSG': states['internal']['PTSG'],  # transporter Glc
            'G6P': states['internal']['G6P'],  # Glc 6P
            'PEP': states['internal']['PEP'],  # Pep
            'PYR': states['internal']['PYR'],  # Pyruvate
            'XP': states['internal']['XP'],  # EIIAP(Pts Protein)
            'GLC[e]': states['external']['GLC'],  # Glc external
            'G6P[e]': states['external']['G6P'],
            'LCTS[e]': states['external']['LCTS'],

            # flux states, for metabolism constraint
            'GLCpts': 0.0,
            'PPS': 0.0,
            'PYK': 0.0,
            'glc__D_e': 0.0,
        }
        state_keys = list(combined_state.keys())
        state_init = np.asarray(list(combined_state.values()))

        # run ode model for t time, get back full solution
        solution = odeint(model, state_init, t)

        # get updates
        internal_update = {}
        external_update = {}
        fluxes = {}
        for state_idx, state_id in enumerate(state_keys):
            if state_id in ['GLC[e]', 'G6P[e]', 'LCTS[e]']:
                # delta counts for external state
                # Note: converting concentrations to counts loses precision
                initial_conc = solution[0, state_idx]
                final_conc = solution[-1, state_idx]
                delta_conc = final_conc - initial_conc

                delta_count = millimolar_to_counts(delta_conc, volume)
                external_update[state_id.replace('[e]','')] = delta_count

            elif state_id in self.target_fluxes:
                # set targets
                # TODO -- units?
                mean_flux = np.mean(solution[:, state_idx])
                fluxes[state_id] = mean_flux

            elif state_id in list(states['internal'].keys()):
                # set internal directly
                internal_update[state_id] = solution[-1, state_idx]

        return {
            'exchange': external_update,
            'internal': internal_update,
            'fluxes': fluxes}

# test and analysis of process
def test_transport(sim_time = 10):
    # Kremling 2007 runs sim for 7.5 hours

    # media for glucose/lactose diauxic growth
    GLC_LCT_shift = {
        'internal': {
            'mass': 0.032,
            'UHPT': 1e-5,
            'PTSG': 0.001,
            'G6P': 0.1,
            'PEP': 0.05,
            'PYR': 0.1,
            'XP': 0.01,
        },
        'external': {
            'GLC': 0.22,
            'G6P': 0.0,
            'LCTS': 1.165,
        }}

    # make the timeline
    timeline = [
        # (0, GLC_LCT_shift),
        (sim_time, {}),
    ]

    # get process, initial state, and saved state
    transport = Transport({})
    settings = transport.default_settings()
    state = settings['state']
    saved_state = {'internal': {}, 'external': {}, 'time': []}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.items():
                    state[key].update(change)

        # get update and apply to state
        update = transport.next_update(timestep, state)
        saved_state['time'].append(time)
        state['internal'].update(update['internal'])

        # use exchange to update external state, reset exchange
        volume = state['global']['volume'] * 1e-15  # convert volume fL to L
        for mol_id, delta_count in update['exchange'].items():
            delta_conc = counts_to_millimolar(delta_count, volume)
            state['external'][mol_id] += delta_conc
            state['exchange'][mol_id] = 0

        # save state
        for port in ['internal', 'external']:
             for state_id, value in state[port].items():
                 if state_id in saved_state[port].keys():
                     saved_state[port][state_id].append(value)
                 else:
                     saved_state[port][state_id] = [value]

    return saved_state

def kremling_figures(saved_state, out_dir='out'):

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [float(t) / 3600 for t in saved_state['time']]  # convert to hours

    # figure 11
    G6P = saved_state['internal']['G6P']
    PEP = saved_state['internal']['PEP']

    # figure 12A -- external state
    biomass = saved_state['internal']['mass']
    GLC_e = saved_state['external']['GLC']
    G6P_e = saved_state['external']['G6P']
    LCTS_e = saved_state['external']['LCTS']

    # figure 12B -- transporters
    PTSG = saved_state['internal']['PTSG']
    LACZ = saved_state['internal']['LACZ']
    UHPT = saved_state['internal']['UHPT']

    # figure 12C -- degree of phosphorylation
    XP = saved_state['internal']['XP']

    # plot results
    n_cols = 1
    n_rows = 4
    plt.figure(figsize=(n_cols * 6, n_rows * 2))

    ax1 = plt.subplot(n_rows, n_cols, 1)
    ax2 = plt.subplot(n_rows, n_cols, 2)
    ax3 = plt.subplot(n_rows, n_cols, 3)
    ax4 = plt.subplot(n_rows, n_cols, 4)

    # figure 11
    ax1.plot(time_vec, G6P, '--', label='G6P')
    ax1.plot(time_vec, PEP, label='PEP')
    ax1.set_xlabel('time (hrs)')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # figure 12A -- external state
    ax2.plot(time_vec, biomass, label='biomass')
    ax2.plot(time_vec, GLC_e, label='GLC_e')
    ax2.plot(time_vec, G6P_e, label='G6P_e')
    ax2.plot(time_vec, LCTS_e, label='LCTS_e')
    ax2.set_xlabel('time (hrs)')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # figure 12B -- transporters
    ax3.plot(time_vec, PTSG, label='PTSG')
    ax3.plot(time_vec, LACZ, label='LACZ')
    ax3.plot(time_vec, UHPT, label='UHPT')
    ax3.set_xlabel('time (hrs)')
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # figure 12C -- degree of phosphorylation
    ax4.plot(time_vec, XP, label='XP')
    ax4.set_xlabel('time (hrs)')
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # save figure
    fig_path = os.path.join(out_dir, 'kremling_fig11')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


def plot_all_state(saved_state, out_dir='out'):

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [float(t) / 3600 for t in saved_state['time']]  # convert to hours

    # make figure, with grid for subplots
    n_rows = 10
    n_data = [len(saved_state[key].keys()) for key in data_keys]
    n_cols = int(math.ceil(sum(n_data) / float(n_rows)))

    fig = plt.figure(figsize=(n_cols * 8, n_rows * 2.5))
    grid = plt.GridSpec(n_rows + 1, n_cols, wspace=0.4, hspace=1.5)

    # plot data
    row_idx = 0
    col_idx = 0
    for key in data_keys:
        for mol_id, series in sorted(saved_state[key].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            ax.set_xlabel('time (hrs)')

            row_idx += 1
            if row_idx > n_rows:
                row_idx = 0
                col_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'transport')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'Kremling2007_transport')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # run simulation
    saved_state = test_transport(7.5*60*60)  # (7.5*60*60)

    kremling_figures(saved_state, out_dir)
    # plot_all_state(saved_state, out_dir)
