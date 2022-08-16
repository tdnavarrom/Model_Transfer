import pandas as pd
import numpy as np
from scipy import integrate

def heat_index(T):
    i  = (T/5)**1.514
    return i

def thornthwaite(T, dates):
    T = np.where(T<0, 0, T)
#     print(T.shape, dates.shape)
    assert T.shape == dates.shape
    months = dates.month
    monthly_means = []
    for m_ix in np.unique(months):
        ixs = np.in1d(months, m_ix)
        monthly_means.append(np.nanmean(T[ixs]))
    # months, counts = np.unique(months, return_counts=True)
    # monthly_normals_daily = np.repeat(a=monthly_means, repeats=counts)
    i = np.nanmean(heat_index(np.array(monthly_means)))
    a = (6.75e-7 * i**3) - (7.71e-5 * i**2) + (1.79e-2 * i) + 0.49
    e = 1.6*(10*T/i)**a # e in centimeters
    return np.where(np.isnan(e), 0, e) #*10 # return in millimeters

def uh_4_full(vol, d_base):
    d_base = np.ceil(d_base)
    tt = np.arange(d_base)
    ff = 0.5 / (0.5*(0.5*d_base)**2)
    d50 = 0.5*d_base
    tri = lambda t: np.max([0, ff*(t-d50)*np.sign(d50-t) + ff*d50])
    UH = np.zeros(tt.shape[0])
    for t in range(1, tt.shape[0]+1):
        UH[t-1] = integrate.quad(tri, t-1, t )[0]
    tmp_diff = 1 - np.nansum(UH)
    tmp_weight = UH / np.nansum(UH)
    UH = UH + (tmp_weight*tmp_diff)
    return vol * UH

class HydrologicModel(object):
    print("running hbv \n")
    def __init__(self, start_date, end_date, timestep, precip, temp, parameters):
        self.start = start_date
        self.end = end_date
        self.timestep = timestep
        self.precip = precip
        self.temp = temp
        self.parameters = parameters
        self.dates = pd.date_range(start=self.start,
                                   periods=self.precip.shape[0],
                                   freq=self.timestep)
        # print("parameters in the HBV model", parameters)
        # assign parameter values
        self.TT = 0 # parameters[0]
        self.TTI = 6 # parameters[1]
        self.CFR = 0.4 # parameters[2]
        self.CFMAX = 0.4 # parameters[3]
        self.TTM = 0 # parameters[4]
        self.WHC = 0.34 # parameters[5]
        self.CFLUX = parameters[0]
        self.FC = parameters[1]
        self.LP = parameters[2]
        self.BETA = parameters[3]
        self.K0 = parameters[4]
        self.ALPHA = parameters[5]
        self.PERC = parameters[6]
        self.K1 = parameters[7]
        self.MAXBAS = parameters[8]


        # create storage vectors
        self.store_S1 = np.zeros(self.dates.shape) # WC
        self.store_S2 = np.zeros(self.dates.shape) # SP
        self.store_S3 = np.zeros(self.dates.shape) # SM
        self.store_S4 = np.zeros(self.dates.shape) # UZ
        self.store_S5 = np.zeros(self.dates.shape) # LZ
        self.states = np.stack([self.store_S1,
                                self.store_S2,
                                self.store_S3,
                                self.store_S4,
                                self.store_S5], axis=-1)
        # create flux vectors
        self.flux_water_in = np.zeros(self.dates.shape) # water input
        self.flux_sf   = np.zeros(self.dates.shape) # snowfall
        self.flux_rf   = np.zeros(self.dates.shape) # rainfall
        self.flux_refr = np.zeros(self.dates.shape) # refreeze
        self.flux_melt = np.zeros(self.dates.shape) # snowmelt
        self.flux_in = np.zeros(self.dates.shape) # infiltration
        self.flux_se = np.zeros(self.dates.shape) # excess water
        self.flux_cf = np.zeros(self.dates.shape) # capillary
        self.flux_r = np.zeros(self.dates.shape) # recharge
        self.flux_ea = np.zeros(self.dates.shape) # actual evaporation
        self.flux_q0 = np.zeros(self.dates.shape) # direct runoff
        self.flux_perc = np.zeros(self.dates.shape) # percolation
        self.flux_q1 = np.zeros(self.dates.shape) # baseflow
        self.flux_q = np.zeros(self.dates.shape) # flow
        self.energy_pet = np.zeros(self.dates.shape) # PET

    def snow(self, flux_rf, flux_sf, T_current, WC_state, SP_state, snow_parameters):
        TTM = snow_parameters[0]
        CFMAX = snow_parameters[1]
        CFR = snow_parameters[2]
        WHC = snow_parameters[3]
        if T_current > TTM:
            if CFMAX*(T_current - TTM) < SP_state + flux_sf:
                MELT = CFMAX * (T_current - TTM)
            else:
                MELT = SP_state + flux_sf
            SP_new = SP_state + flux_sf - MELT
            WC_int = WC_state + MELT + flux_rf
            REFR = 0.0
        else:
            if CFR*CFMAX*(TTM-T_current) < WC_state + flux_rf:
                REFR = CFR*CFMAX*(TTM-T_current)
            else:
                REFR = WC_state + flux_rf
            SP_new = SP_state + flux_sf + REFR
            WC_int = WC_state - REFR + flux_rf
            MELT = 0.0
        if WC_int > WHC*SP_new:
            IN = WC_int - WHC*SP_new
            WC_new = WHC*SP_new
        else:
            IN = 0
            WC_new = WC_int
        return IN, WC_new, SP_new, MELT, REFR

    def soil(self, flux_in, flux_ep, SM_state, UZ_state, parameters):
        # flux_in = 2000
        LP = parameters[0]
        BETA = parameters[1]
        FC = parameters[2]
        CFLUX = parameters[3]
        # print("sm_state", SM_state, type(SM_state))
        # print("flux_in", flux_in, type(flux_in))
        # print("fc", FC, type(FC))
        # val = SM_state + flux_in - FC
        # print("val", val, type(val))
        ponding = np.max([ 0, SM_state + flux_in - FC ])
        _in = flux_in - ponding
        # print("in", _in)
        RECH = ((SM_state/FC)**BETA) * _in
        # print("RECH", RECH)
        # print("Ep", flux_ep)
        # print("state sm", SM_state)
        # print("LP", LP)
        # print("FC", FC)
        # print()
        EA = self.evaporation(Ep=flux_ep, SM_state=SM_state, FC=FC, LP=LP)
        # np.min([ flux_ep, flux_ep*(SM_state/(LP*FC)) ])
        # print("ea", EA)
        # use ponding to attain energy limit
        # free surface et
        extra_et = np.min([ flux_ep, ponding ])
        # print(ponding, extra_et, EA)
        QDR = ponding #- extra_et # ponding balance to become direct runoff
        # print("QDR", QDR)
        C_FLUX = np.min( [UZ_state, CFLUX * ((FC - SM_state)/FC)] )
        # print("CFLUX", C_FLUX)
        SM_new = np.max([ 0, SM_state + _in - RECH + C_FLUX - EA])
        UZ_int = np.max([0, UZ_state + RECH - C_FLUX + QDR ])
        # print("UZ_int1", UZ_int)
        return SM_new, UZ_int, QDR, RECH, EA, C_FLUX, extra_et

    def response(self, UZ_int, LZ_state, parameters):
        PERC = parameters[0]
        ALPHA = parameters[1]
        K0 = parameters[2]
        K1 = parameters[3]
        flux_P = np.min([UZ_int, PERC])
        LZ_int = LZ_state + flux_P
        UZ_int = np.max([ 0, UZ_int - flux_P ])
        # print("UZ_int2", UZ_int)
        Q0 = np.min([UZ_int, K0*(UZ_int**(1.+ALPHA))])
        # print("Q0", Q0)
        Q1 = np.min([LZ_int, K1*LZ_int])
        UZ_new = np.max([0, UZ_int - Q0])
        # print(UZ_new, UZ_int)
        LZ_new = np.max([0, LZ_int - Q1])
        return Q0, Q1, flux_P, UZ_new, LZ_new


    def snowfall(self, P, T_current, TT, TTI):
        return np.min([P, np.max([0, P*((TT + 0.5*TTI - T_current)/(TTI))])])

    def rainfall(self, P, T_current, TT, TTI):
        return np.min([P, np.max([0, P*((T_current - (TT-0.5*TTI))/(TTI))])])

    def infiltration(self, INPUTS, SP_state, WC_state, WHC):
        return np.sum(INPUTS) * np.max([0, np.sign(WC_state - WHC*SP_state)])

    def excess_water(self, WC_state, SP_state, WHC):
        return WC_state - WHC*SP_state * np.max([0, np.sign(WC_state - WHC*SP_state)])

    def capillary_rise(self, SM_state, FC, CFLUX, UZ_state):
        return np.min([ UZ_state, CFLUX * ( 1 - (SM_state/FC) ) ])

    def evaporation(self, Ep, SM_state, FC, LP):
        return np.min([Ep*(SM_state/(FC*LP)), Ep, SM_state])

    def recharge(self, INPUTS, SM_state, FC, BETA):
        # INPUT is infiltration + excess_WATER
        return np.sum(INPUT) * (np.max([0, SM_state])/FC)**BETA

    def interflow(self, UZ_state, K0, ALPHA): # Q0
        return np.min([ K0*UZ_state**(1+ALPHA), np.max([0, UZ_state])])

    def percolation(self, UZ_state, C):
        return np.min([UZ_state, C])

    def baseflow(self, LZ_state, K1): # Q1
        return LZ_state * K1

    def mass_balance(self, INPUTS, OUTPUTS, INIT_STATE, verbose=False):
        if verbose:
            print("IN", INPUTS, np.nansum(INPUTS))
            print("OUT", OUTPUTS, np.nansum(OUTPUTS))
            print("INIT",INIT_STATE)
            print(np.max([0, INIT_STATE + np.nansum(INPUTS) - np.nansum(OUTPUTS)]))
            print()
            NEW_STATE = INIT_STATE + np.nansum(INPUTS) - np.nansum(OUTPUTS)
            assert ( NEW_STATE - INIT_STATE ) - (np.nansum(INPUTS) - np.nansum(OUTPUTS)) < 1e-6
        return np.max([0, INIT_STATE + np.sum(INPUTS) - np.sum(OUTPUTS)])

    def step_run(self, old_states, input_p, input_t, timestep):
        """
        THIS FUNCTION WILL
        RUN ONE FORWARD STEP
        OF THE MODEL
        """
        # set store to old_states
        S1_old = old_states[0] # WC_state #.*self.x1
        S2_old = old_states[1] # SP_state
        S3_old = old_states[2] # SM_state
        S4_old = old_states[3] # UZ_state
        S5_old = old_states[4] # LZ_state
        # print("soil_states", S1_old, S2_old, S3_old, S4_old, S5_old)

        P = input_p
        T = input_t
        CURRENT_DATE = pd.to_datetime([self.dates[timestep]])
        PET = thornthwaite(T=np.array([T]), dates=CURRENT_DATE)
        self.energy_pet[timestep] = PET[0]
        # print(PET)
        # partition precip into rain and snow
        RF = self.rainfall(P=P, T_current=T, TT=self.TT, TTI=self.TTI)
        SF = self.snowfall(P=P, T_current=T, TT=self.TT, TTI=self.TTI)
        # run snow module
        IN, S1_new, S2_new, MELT, REFR = self.snow(flux_rf=RF,
                                              flux_sf=SF,
                                              T_current=T,
                                              WC_state=S1_old,
                                              SP_state=S2_old,
                                              snow_parameters=[self.TTM, self.CFMAX, self.CFR, self.WHC])
        # run top soil module
        S3_new, UZ_int, QDR, RECH, EA, C_FLUX, extra_et = self.soil(flux_in=IN,
                                                            flux_ep=PET[0],
                                                            SM_state=S3_old,
                                                            UZ_state=S4_old,
                                                            parameters=[self.LP, self.BETA, self.FC, self.CFLUX])
        Q0, Q1, PERC_FLUX, S4_new, S5_new = self.response( UZ_int=UZ_int,
                                               LZ_state=S5_old,
                                               parameters=[self.PERC, self.ALPHA, self.K0, self.K1])
        # UPDATE STATE VECTORS
        self.states[timestep, :] = [S1_new, S2_new, S3_new, S4_new, S5_new]
        # UPDATE FLUX VECTORS
        self.flux_water_in[timestep] = input_p
        self.flux_sf[timestep] = SF
        self.flux_rf[timestep] = RF
        self.flux_refr[timestep] = REFR
        self.flux_melt[timestep] = MELT
        self.flux_in[timestep] = IN
        self.flux_ea[timestep] = EA
        self.flux_r[timestep] = RECH
        self.flux_cf[timestep] = C_FLUX
        self.flux_q0[timestep] = Q0; # print("q0:", Q0)
        self.flux_q1[timestep] = Q1; # print("q1:", Q1)
        self.flux_perc[timestep] = PERC_FLUX
        self.flux_se[timestep] = extra_et
        self.flux_q[timestep] = Q1+Q0 # final runoff
        return

    def simulate(self):
        for t in range(self.dates.shape[0]):
            if t == 0:
                old_states = np.zeros(shape=self.states.shape[1])
            else:
                old_states = self.states[t-1]
#             print(old_states)
            T_current = self.temp[t]
            P = self.precip[t]
            self.step_run(old_states=old_states, input_p=P, input_t=T_current, timestep=t)
        UH = uh_4_full(vol=1, d_base=self.MAXBAS)
        self.flux_q = np.convolve(UH, self.flux_q, 'same')
        self.flux_q = np.where(np.isnan(self.flux_q), 0, self.flux_q)
        # print(self.flux_q)
        self.flux_ea = np.where(np.isnan(self.flux_ea), 0, self.flux_ea)
        return
