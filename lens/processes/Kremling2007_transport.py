from __future__ import absolute_import, division, print_function

import os
import numpy as np
from scipy.integrate import odeint

from lens.actor.process import Process
from lens.utils.flux_conversion import millimolar_to_counts, counts_to_millimolar
from lens.environment.make_media import Media
from lens.utils.dict_utils import merge_dicts


DEFAULT_PARAMETERS = {
    # enzyme synthesis
    'k1': 0.00001,  # enzyme synthesis k1[\mu mo / gDW h]
    'k2': 0.0001,  # enzyme synthesis k2[\mu mo / gDW h]
    'k3': 0.00016,  # enzyme synthesis k3[\mu mo / gDW h]
    'K1': 3000,  # K1[\mu mo / gDW]
    'K2': 2800,  # K2[\mu mo / gDW]
    'K3': 15000,  # K3[\mu mo / gDW]
    # enzyme degration
    'kd': 0.4,  # enzyme degradation k_d[1 / h]

    'm': 1,  # m[-]
    'n': 2,  # n[-]
    'x0': 0.1,  # EIIA_0[\mu mo / gDW]

    'kg6p': 2.8e6,  # Glc 6P transporter k_g6p[1 / h]
    'Kg6p': 0.1,  # Glc 6P transporter K_g6p[g / l]

    'kptsup': 2.7e8,  # Glc transporter k_pts_up[1 / h]
    'Kglc': 0.12,  # Glc transporter K_glc[g / l]
    'Keiiap': 12,  # Glc transporter K_EIIAP[\mu mol / gDW]

    'klac': 5.4e5,  # Lac transporter k_lac[1 / h]
    'Klac': 0.13,  # Lac transporter K_lac[g / l]
    'Kieiia': 5.0,  # K_IEIIA: inhibition of Lac transporter by EIIA[-]

    'kgly': 2.80e4,  # k_gly[1 / h]
    'kpyk': 9.39e5,  # k_pyk[1 / (\mu mol / gDW) ^ 3 h]
    'kpdh': 5.50e3,  # k_pdh[1 / h]
    'kpts': 1.86e5,  # k_pts[1(\mu mol / gDW) g]
    'Kpts': 0.7,  # K_pts[-]
    'km_pts': 0.7 * 1.86e5,  # params['Kpts'] * params['kpts']  #

    'Y': 1.0e-4,  # Yglc[g TM / mu mol Substrate]

    'mw1': 2.602e-4,  # molecular weight g6p[g Substrate / mu mol Substrate]
    'mw2': 1.802e-4,  # molecular weight glc[g Substrate / mu mol Substrate]
    'mw3': 3.423e-4,  # molecular weight lac[g Substrate / mu mol Substrate]
    'Y1_sim': 6.2448e-05,  # Yg6p[gDW /\mu mol]
    'Y2_sim': 1.0e-4,  # Yglc[gDW /\mu mol]
    'Y3_sim': 9.2421e-05,  # Ylac[gDW /\mu mol]

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
    def __init__(self, initial_parameters={}):
        self.dt = 0.001  # timestep for ode integration (seconds)

        default_settings = self.default_settings()
        default_state = default_settings['state']
        internal_state = default_state['internal']
        external_state = default_state['external']

        roles = {
            'external': external_state.keys(),
            'exchange': external_state.keys(),
            'internal': internal_state.keys()}
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(Transport, self).__init__(roles, parameters)

    def default_settings(self):

        # default state
        # TODO -- select state based on media
        glc_g6p = True
        glc_lct = False

        make_media = Media()
        if glc_g6p:
            external = make_media.get_saved_media('GLC_G6P')
            internal = {
                'mass': 0.032,  # [g / l]  #1.339e-12,  # g
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
        self.environment_ids = external.keys()

        default_state = {
            'internal': merge_dicts([internal, {'volume': 1}]),  # TODO -- get volume with deriver?
            'external': external,
            'exchange': {state_id: 0.0 for state_id in self.environment_ids}
        }

        # default emitter keys
        default_emitter_keys = {
            'internal': ['mass', 'UHPT', 'PTSG', 'G6P', 'PEP', 'PYR', 'XP'],
            'external': ['G6P', 'GLC', 'LAC']}

        # default updaters
        default_updaters = {
            'internal': {state_id: 'set'
                         for state_id in ['mass', 'UHPT', 'PTSG', 'G6P', 'PEP', 'PYR', 'XP']},  # reactions set values directly
            'external': {},  # reactions set values directly
            'exchange': {mol_id: 'accumulate' for mol_id in self.environment_ids}}  # all external values use default 'delta' udpater


        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters}

        return default_settings

    def next_update(self, timestep, states):

        def model(state, t):
            '''
            Model of sugar transport based on Kremling, Bettenbrock, & Gilles. (2007).
            Analysis of global control of Escherichia coli carbohydrate uptake
            '''

            biomass = state[state_keys.index('mass')]  # biomass
            sugar2 = state[state_keys.index('GLC[e]')]  # GLC external
            transporter2 = state[state_keys.index('PTSG')]  # transporter GLC
            metabolite1 = state[state_keys.index('G6P')]  # Glc6P
            metabolite2 = state[state_keys.index('PEP')]  # PEP
            metabolite3 = state[state_keys.index('PYR')]  # Pyruvate
            protein_p = state[state_keys.index('XP')]  # EIIAP (Pts Protein)

            g6pext_present = state[state_keys.index('G6P')] > 0.01
            if g6pext_present:
                sugar1 = state[state_keys.index('G6P[e]')]  # G6P external
                transporter1 = state[state_keys.index('UHPT')]  # transporter for G6P
            else:
                sugar1 = state[state_keys.index('LCTS[e]')]  # LCTS external
                transporter1 = state[state_keys.index('LACZ')]  # transporter for LCTS

            ## Original rates.
            if g6pext_present:
                # G6P uptake
                rup1 = p['kg6p'] * (transporter1 * sugar1) / (p['Kg6p'] + sugar1)
            else:
                # Lactose uptake
                rup1 = p['klac'] * (transporter1 * sugar1) / (p['Klac'] + sugar1 * (1 + ((p['x0'] - protein_p) / p['x0']) / p['Kieiia']))

            # PTS uptake. eqn 38
            rup2 = p['kptsup'] * protein_p * (transporter2 * sugar2) / (p['Kglc'] * p['Keiiap'] * p['x0'] + sugar2 * p['Keiiap'] * p['x0'] + protein_p * p['Kglc'] + protein_p * sugar2)

            # enzyme syn 1 / 2
            if g6pext_present:
                rsyn1 = p['k1'] * (
                        p['kb'] + p['ksyn'] * protein_p**6. / (protein_p**6 + p['K'] ** 6)) * rup1 / (p['K1'] + rup1)
                rsyn2 = p['k2'] * (p['KI'] / (transporter1 + p['KI'])) * (p['kb'] + p['ksyn'] * protein_p**6. / (protein_p**6 + p['K']**6)) * rup2 / (p['K2'] + rup2)
            else:
                rsyn1 = p['k3'] * (p['kb'] + p['ksyn'] * (protein_p)**6. / (protein_p**6 + p['K']**6)) * rup1 / (p['K3'] + rup1)
                rsyn2 = p['k2'] * (p['kb'] + p['ksyn'] * (protein_p)**6. / (protein_p**6 + p['K']**6)) * rup2 / (p['K2'] + rup2)

            ## Rates
            rgly = p['kgly'] * metabolite1  # Glycolyse
            f = (metabolite1**p['n']) * metabolite2**p['m']  #
            rpyk = p['kpyk'] * metabolite2 * f  # Pyk
            rpts = p['kpts'] * metabolite2 * (p['x0'] - protein_p) - p['km_pts'] * metabolite3 * protein_p  # PTS rate
            rpdh = p['kpdh'] * metabolite3  # Pdh

            # growth rate. eqn 44
            if g6pext_present:
                mu = p['Y1_sim'] * rup1 + p['Y2_sim'] * rup2
            else:
                mu = p['Y3_sim'] * rup1 + p['Y2_sim'] * rup2

            ppc_rate = 0  # TODO -- figure out code below
            # if (nargin >= 7 & & ~isempty(varargin{1}) & & ~isempty(varargin{2})):
            #     FBA_primal = varargin{1}
            #     FBA_rxnInd = varargin{2}
            #     options = varargin{3}
            #     ppc_rate = (FBA_primal(FBA_rxnInd.PPC) - FBA_primal(FBA_rxnInd.PCKA)) * 1e3;
            #     if (nargin == 7 | | (nargin == 8 & & varargin{4}~=0))
            #         mu = options.iFBA_growRateScale * FBA_primal(FBA_rxnInd.VGRO);
            # else:
            #     ppc_rate = 0;

            # initialize sugar1, transporter1
            dG6P_e = 0.0
            dLCTS_e = 0.0
            dUHPT = 0.0
            dLACZ = 0.0

            if g6pext_present:
                dG6P_e = -p['mw1'] * rup1 * biomass
                dUHPT = rsyn1 - (p['kd'] + mu) * transporter1  # transporter. eqn 46
            else:
                dLCTS_e = -p['mw3'] * rup1 * biomass
                dLACZ = rsyn1 - (p['kd'] + mu) * transporter1  # transporter. eqn 46

            dbiomass = mu * biomass  # mass
            dGLC_e = -p['mw2'] * rup2 * biomass  # GLC[e]
            dPTSG = rsyn2 - (p['kd'] + mu) * transporter2  # transporter. eqn 46
            dG6Pdt = rup1 + rup2 - rgly # eqn 2
            dPEPdt = 2 * rgly - rpyk - rpts - ppc_rate  # eqn 3
            dPYRdt = rpyk + rpts - rpdh   # eqn 4
            dprotein_p = rpts - rup2  # eqn 5 EIIAP (Pts Protein)

            # save to numpy array
            dx = np.zeros_like(state)
            dx[state_keys.index('mass')] = dbiomass
            dx[state_keys.index('UHPT')] = dUHPT  # transporter1 changes with condition
            dx[state_keys.index('LACZ')] = dLACZ  # transporter1 changes with condition
            dx[state_keys.index('PTSG')] = dPTSG
            dx[state_keys.index('G6P')] = dG6Pdt
            dx[state_keys.index('PEP')] = dPEPdt
            dx[state_keys.index('PYR')] = dPYRdt
            dx[state_keys.index('XP')] = dprotein_p
            dx[state_keys.index('GLC[e]')] = dGLC_e
            dx[state_keys.index('G6P[e]')] = dG6P_e    # sugar1 changes with condition
            dx[state_keys.index('LCTS[e]')] = dLCTS_e  # sugar1 changes with condition

            return dx

        # set up state and parameters for odeint
        timestep_hours = timestep / 3600
        dt_hours = self.dt / 3600
        p = self.parameters
        t = np.arange(0, timestep_hours, dt_hours)

        # get states
        volume = states['internal']['volume'] * 1e-15  # convert volume fL to L
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
        }
        state_keys = combined_state.keys()
        state_init = np.asarray(combined_state.values())

        # run ode model for t time, get back full solution
        solution = odeint(model, state_init, t)

        internal_update = {}
        external_update = {}
        for state_idx, state_id in enumerate(state_keys):
            initial_conc = solution[0, state_idx]
            final_conc = solution[-1, state_idx]
            delta_conc = final_conc - initial_conc

            if state_id in ['GLC[e]', 'G6P[e]', 'LCTS[e]']:
                # Get delta counts for external state
                # Note: converting concentrations to counts can lose precision
                delta_count = millimolar_to_counts(delta_conc, volume)
                external_update[state_id.replace('[e]','')] = delta_count

            elif state_id is 'volume':
                # skip volume
                continue
            else:
                # set internal directly
                internal_update[state_id] = final_conc

        # TODO -- make an optional 'external' update to set external concentrations directly.  set updater to null in composition.
        return {
            'exchange': external_update,
            'internal': internal_update}

# test and analysis of process
def test_transport():
    timeline = [
        # (0, {'external': {
        #     'GLC': 1}
        # }),
        (7.5*60*60, {}),  # Kremling runs sim for 7.5 hours
    ]

    # configure process
    transport = Transport({})

    # get initial state and parameters
    settings = transport.default_settings()
    state = settings['state']
    saved_state = {'internal': {}, 'exchange': {}, 'time': []}

    # run simulation
    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.iteritems():
                    state[key].update(change)

        # get update
        update = transport.next_update(timestep, state)

        # apply update to state
        saved_state['time'].append(time)
        state['internal'].update(update['internal'])

        # apply external update
        volume = state['internal']['volume'] * 1e-15  # convert volume fL to L
        for mol_id, delta_count in update['exchange'].iteritems():
            delta_conc = counts_to_millimolar(delta_count, volume)
            state['exchange'][mol_id] += delta_conc

        # save state
        for role in ['internal', 'exchange']:
             for state_id, value in state[role].iteritems():
                 if state_id in saved_state[role].keys():
                     saved_state[role][state_id].append(value)
                 else:
                     saved_state[role][state_id] = [value]

    return saved_state

def plot_transport(saved_state, out_dir='out'):
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
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            ax.set_xlabel('time (hrs)')

            plot_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'transport')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_transport()
    out_dir = os.path.join('out', 'tests', 'Kremling2007_transport')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_transport(saved_state, out_dir)
