from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.integrate import odeint

from lens.actor.process import Process
from lens.utils.flux_conversion import millimolar_to_counts

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

# GLC-G6P media
GLC_G6P_EXTERNAL = {
    # from covert2002
    'ACET': 0.0,
    'CO+2': 100.0,  # "units.mmol / units.L"
    'ETOH': 0.0,
    'FORMATE': 0.0,
    'GLYCEROL': 0.0,
    'LAC': 0.0,
    'LCTS': 0.0,
    'OXYGEN-MOLECULE': 100.0,  # "units.mmol / units.L"
    'PI': 100.0,  # "units.mmol / units.L"
    'PYR': 0.0,
    'RIB': 0.0,
    'SUC': 0.0,
    # from kremling2007
    'G6P': 1.3451,  # [m mol/L]
    'GLC': 12.2087,  # [m mol/L]
}
GLC_G6P_INTERNAL = {
    'Biomass': 0.032,  # [g / l]
    'UHPT': 0.0003,  # [\mu mol gDCW] enz g6p
    'PTSG': 0.007,  # [\mu mol gDCW] enz glc
    'G6P': 0.2057,  # [\mu mol gDCW]
    'PEP': 2.0949,  # [\mu mol gDCW]
    'PYR': 2.0949,  # [\mu mol gDCW]
    'XP': 0.0038,  # [fraction phosphorylation]
}

# GLC-LCT media
GLC_LCT_EXTERNAL = {
    'G6P': 0,  # [m mol / L]
    'GLC': 1.2209,  # [m mol / L]
    'RIB': 0,  # [m mol / L]
    'GLYCEROL': 0,  # [m mol / L]
    'SUC': 0,  # [m mol / L]
    'PYR': 0,  # [m mol / L]
    'LAC': 0,  # [m mol / L]
    'LCTS': 3.4034,  # [m mol / L]
    'FORMATE': 0,  # [m mol / L]
    'ETOH': 0,  # [m mol / L]
    'ACET': 0,  # [m mol / L]
    'PI': 100,  # [m mol / L]
    'CO+2': 100,  # [m mol / L]
    'OXYGEN-MOLECULE': 100,  # [m mol / L]
}
GLC_LCT_INTERNAL = {
    'Biomass': 0.032,  # [g / l]
    'LACZ': 1e-5,  # [\mu mol gDCW] enz g6p
    'PTSG': 0.001,  # [\mu mol gDCW] enz glc
    'G6P': 0.1,  # [\mu mol gDCW]
    'PEP': 0.05,  # [\mu mol gDCW]
    'PYR': 0.1,  # [\mu mol gDCW]
    'XP': 0.01,  # [fraction phosphorylation]
}


def merge_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


class Transport(Process):
    def __init__(self, initial_parameters={}):
        self.exchange_key = initial_parameters['exchange_key']
        self.dt = 0.01  # timestep for ode integration

        default_state = self.default_state()
        internal_state = default_state['internal']
        external_state = default_state['external']
        self.environment_ids = default_state['environment_ids']

        roles = {
            'external': external_state.keys(),
            'internal': internal_state.keys()}
        parameters = DEFAULT_PARAMETERS
        parameters.update(initial_parameters)

        super(Transport, self).__init__(roles, parameters)

    def default_state(self):
        glc_g6p = True
        glc_lct = False

        if glc_g6p:
            external = GLC_G6P_EXTERNAL
            internal = GLC_G6P_INTERNAL
        elif glc_lct:
            external = GLC_LCT_EXTERNAL
            internal = GLC_LCT_INTERNAL

        environment_ids = external.keys()
        external_changes = [key + self.exchange_key for key in environment_ids]

        # declare the states
        external_molecules = merge_dicts(external, {key: 0 for key in external_changes})
        internal_molecules = merge_dicts(internal, {'volume': 1}) # fL TODO -- get volume with deriver?

        return {
            'environment_deltas': external_changes,
            'environment_ids': environment_ids,
            'external': external_molecules,
            'internal': internal_molecules}

    def default_emitter_keys(self):
        keys = {
            'internal': ['Biomass', 'UHPT', 'PTSG', 'G6P', 'PEP', 'PYR', 'XP'],
            'external': ['G6P','GLC']
        }
        return keys


    def next_update(self, timestep, states):

        def model(state, t):
            '''
            Model of sugar transport based on Kremling, Bettenbrock, & Gilles. (2007).
            Analysis of global control of Escherichia coli carbohydrate uptake
            '''

            G6Pext_present = state[state_keys.index('G6Pext_present')]
            XX = state[state_keys.index('Biomass')]  # biomass
            S1 = state[state_keys.index('S1')]  # G6P or LCTS, depending on environment
            S2 = state[state_keys.index('GLC[e]')]  # GLC external
            E1 = state[state_keys.index('E1')]  # UHPT or LACZ, depending on environment
            E2 = state[state_keys.index('PTSG')]  # transporter GLC
            M1 = state[state_keys.index('G6P')]  # Glc6P
            M2 = state[state_keys.index('PEP')]  # PEP
            M3 = state[state_keys.index('PYR')]  # Pyruvate
            XP = state[state_keys.index('XP')]  # EIIAP (Pts Protein)

            ## Original rates. kremling_rates.m
            if G6Pext_present:
                # G6P uptake
                rup1 = p['kg6p'] * (E1 * S1) / (p['Kg6p'] + S1)
            else:
                # Lactose uptake
                rup1 = p['klac'] * (E1 * S1) / (p['Klac'] + S1 * (1 + ((p['x0'] - XP) / p['x0']) / p['Kieiia']))

            rup2 = p['kptsup'] * XP * (E2 * S2) / (p['Kglc'] * p['Keiiap'] * p['x0'] + S2 * p['Keiiap'] * p['x0'] + XP * p['Kglc'] + XP * S2)

            # enzyme syn 1 / 2
            if G6Pext_present:
                rsyn1 = p['k1'] * (
                        p['kb'] + p['ksyn'] * (XP)**6. / (XP**6 + p['K'] ** 6)) * rup1 / (p['K1'] + rup1)
                rsyn2 = p['k2'] * (p['KI'] / (E1 + p['KI'])) * (p['kb'] + p['ksyn'] * (XP)**6. / (XP**6 + p['K']**6)) * rup2 / (p['K2'] + rup2)
            else:
                rsyn1 = p['k3'] * (p['kb'] + p['ksyn'] * (XP)**6. / (XP**6 + p['K']**6)) * rup1 / (p['K3'] + rup1)
                rsyn2 = p['k2'] * (p['kb'] + p['ksyn'] * (XP)**6. / (XP**6 + p['K']**6)) * rup2 / (p['K2'] + rup2)

            ## Rates. kremling_rates.m
            rgly = p['kgly'] * M1  # Glycolyse
            f = (M1**p['n']) * M2**p['m']  #
            rpyk = p['kpyk'] * M2 * f  # Pyk
            rpts = p['kpts'] * M2 * (p['x0'] - XP) - p['km_pts'] * M3 * XP  # PTS rate
            rpdh = p['kpdh'] * M3  # Pdh

            # additional rates for iFBA
            if G6Pext_present:
                mu = p['Y1_sim'] * rup1 + p['Y2_sim'] * rup2  # growth rate
            else:
                mu = p['Y3_sim'] * rup1 + p['Y2_sim'] * rup2  # growth rate

            ppcRate = 0  # TODO -- figure out code below
            # if (nargin >= 7 & & ~isempty(varargin{1}) & & ~isempty(varargin{2})):
            #     FBA_primal = varargin{1}
            #     FBA_rxnInd = varargin{2}
            #     options = varargin{3}
            #     ppcRate = (FBA_primal(FBA_rxnInd.PPC) - FBA_primal(FBA_rxnInd.PCKA)) * 1e3;
            #     if (nargin == 7 | | (nargin == 8 & & varargin{4}~=0))
            #         mu = options.iFBA_growRateScale * FBA_primal(FBA_rxnInd.VGRO);
            # else:
            #     ppcRate = 0;

            ## ODE model. kremling_ode_model.m
            if G6Pext_present:
                carbo1 = -p['mw1'] * rup1 * XX
            else:
                carbo1 = -p['mw3'] * rup1 * XX
            carbo2 = -p['mw2'] * rup2 * XX

            dBiomass = mu * XX  # biomass
            dS1 = carbo1 # G6P[e] / LCTS[e] TODO -- in LCTS, this is dLCTS[e]
            dS2 = carbo2 # GLC[e]
            dE1 = rsyn1 - (p['kd'] + mu) * E1  # transporter. eqn 46
            dE2 = rsyn2 - (p['kd'] + mu) * E2  # transporter. eqn 46
            dG6Pdt = rup1 + rup2 - rgly # eqn 2
            dPEPdt = 2 * rgly - rpyk - rpts - ppcRate  # eqn 3
            dPYRdt = rpyk + rpts - rpdh   # eqn 4
            dXP = rpts - rup2  # eqn 5 EIIAP (Pts Protein)

            # save to numpy array
            dx = np.zeros_like(state)
            dx[state_keys.index('Biomass')] = dBiomass
            dx[state_keys.index('S1')] = dS1  # S1_key changes with condition
            dx[state_keys.index('GLC[e]')] = dS2
            dx[state_keys.index('E1')] = dE1  # E1_key changes with condition
            dx[state_keys.index('PTSG')] = dE2
            dx[state_keys.index('G6P')] = dG6Pdt
            dx[state_keys.index('PEP')] = dPEPdt
            dx[state_keys.index('PYR')] = dPYRdt
            dx[state_keys.index('XP')] = dXP

            return dx


        p = self.parameters
        t = np.arange(0, timestep, self.dt)

        # get states
        # s1, e1 keys depend on G6P[e] presence
        G6Pext_present = states['external'].get('G6P',0)>0.01
        if G6Pext_present:
            s1_key = 'G6P'  # Glc 6P external
            e1_key = 'UHPT'  # transporter Glc 6P
        else:
            s1_key = 'LCTS'  # LCTS external
            e1_key = 'LACZ'  # transporter LCTS

        combined_state = {
            'G6Pext_present': G6Pext_present,
            'Biomass': states['internal']['Biomass'],  # biomass
            'E1': states['internal'][e1_key],
            'PTSG': states['internal']['PTSG'],  # transporter Glc
            'G6P': states['internal']['G6P'],  # Glc 6P
            'PEP': states['internal']['PEP'],  # Pep
            'PYR': states['internal']['PYR'],  # Pyruvate
            'XP': states['internal']['XP'],  # EIIAP(Pts Protein)
            'S1': states['external'][s1_key],
            'GLC[e]': states['external']['GLC'],  # Glc external
        }
        state_keys = combined_state.keys()
        state_init = np.asarray(combined_state.values())

        # run ode model for t time, get back full solution
        solution = odeint(model, state_init, t)


        environment_delta_counts = {}
        internal_update = {}

        # TODO -- get difference between initial and final state

        for state_idx, state_id in enumerate(state_keys):
            series = solution[:, state_idx]

            import ipdb; ipdb.set_trace()

        # for env_id in self.environment_ids:
        #     flux = dx['external'].get(env_id, 0.) * timestep_hrs # flux is in mmol/L/hr
        #
        #     if np.isnan(flux):
        #         environment_delta_counts[env_id] = 0.0
        #     else:
        #         # TODO -- should this be an accumulate_deltas?
        #         environment_delta_counts[env_id] = millimolar_to_counts(flux, volume)
        #
        # internal_update = {mol_id: dx['internal'].get(mol_id,0.) * timestep_hrs for mol_id in states['internal'].iterkeys()}

        # TODO -- 1 sec timestep too high -- values drop to (-)

        return {
            'external': {mol_id + self.exchange_key: delta for mol_id, delta in environment_delta_counts.iteritems()},
            'internal': internal_update}
