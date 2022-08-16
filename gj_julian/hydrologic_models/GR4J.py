import pandas as pd
import numpy as np
from scipy import integrate
import math
from math import tanh

def uh_1_half(vol, d_base):
    d_base = int(np.ceil(d_base))
    tt = d_base
    SH = (np.arange(0, tt+1)/d_base)**(5/2)
    UH = np.diff(SH)
    return vol * UH

def uh_2_full(vol, d_base):
    d_base = int(np.ceil(d_base)); # print(d_base)
    tt = d_base * 2
    SH = np.arange(1, tt+1)/np.ceil(d_base)
    SH[:d_base] = 0.5*(SH[:d_base] ** (5./2))
    SH[d_base:] = 0.5*(2 - SH[d_base:])**(5./2)
    UH = SH/np.nansum(SH)
    return vol * UH

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

def heat_index(T):
    i  = (T/5)**1.514
    return i

def thornthwaite(T, dates):
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
    return e #*10 # return in millimeters

class GR4J(object):
    def __init__(self, start_date, end_date, timestep, precip, temp, parameters):
        self.start = start_date
        self.end = end_date
        self.timestep = timestep
        self.precip = precip
        self.temp = temp
        assert np.all(np.isfinite(temp))
        self.parameters = parameters
        self.dates = pd.date_range(start=self.start,
                                   end=self.end,
                                   freq=self.timestep); # print(self.dates.shape[0])
        # assign parameter values
        self.x1 = parameters[0]
        self.x2 = parameters[1]
        self.x3 = parameters[2]
        self.x4 = parameters[3]

        # create storage vectors
        self.store_S = np.zeros(self.dates.shape)
        self.store_R = np.zeros(self.dates.shape)
        self.store_S1 = np.zeros(self.precip.shape) # SNOW PACK
        self.store_S2 = np.zeros(self.precip.shape) # LIQUID CONTENT OF SNOW PACK
        # create flux vectors
        self.flux_ps   = np.zeros(self.dates.shape)#
        self.flux_pq   = np.zeros(self.dates.shape)
        self.flux_es   = np.zeros(self.dates.shape)#
        self.flux_perc = np.zeros(self.dates.shape)#
        self.flux_q9   = np.zeros(self.dates.shape)#
        self.flux_q1   = np.zeros(self.dates.shape)#
        self.flux_qr   = np.zeros(self.dates.shape)#
        self.flux_q   = np.zeros(self.dates.shape)#
        self.flux_f    = np.zeros(self.dates.shape)#
        self.error_S = np.zeros(self.dates.shape)
        self.error_R = np.zeros(self.dates.shape)
        # prepare unit hydrograph
        # self.uh_half = uh_1_half(1, self.x4)
        # self.uh_full = uh_2_full(1, self.x4*2)
        # initialize routing vectors

    def snow(self, flux_rf, flux_sf, T_current, WC_state, SP_state):
        TTM = 0
        CFMAX = 0.4
        CFR = 0.4
        WHC = 0.34
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


    def snowfall(self, P, T_current):
        TT = 0
        TTI = 6
        return np.min([P, np.max([0, P*((TT + 0.5*TTI - T_current)/(TTI))])])

    def rainfall(self, P, T_current):
        TT = 0
        TTI = 6
        return np.min([P, np.max([0, P*((T_current - (TT-0.5*TTI))/(TTI))])])

    def step_run(self, old_states, input_p, input_t, timestep):
        S_old, R_old, S1_old, S2_old = old_states

        CURRENT_DATE = pd.to_datetime([self.dates[timestep]])
        PET = thornthwaite(T=np.array([input_t]),dates=CURRENT_DATE)
        # partition precip into rain and snow
        RF = self.rainfall(P=input_p, T_current=input_t)
        SF = self.snowfall(P=input_p, T_current=input_t)
        # run snow module
        IN, S1_new, S2_new, MELT, REFR = self.snow(flux_rf=RF,
                                                   flux_sf=SF,
                                                   T_current=input_t,
                                                   WC_state=S1_old,
                                                   SP_state=S2_old)
        PN = np.max([IN - PET, 0])
        EN = np.max([PET - IN, 0])

        if S_old > self.x1:
            self.error_S[t] = S_old - self.x1
            S_old = self.x1*0.9999
        if R_old > self.x3:
            self.error_R[t] = R_old - self.x3
            R_old = self.x3*0.9999
        PS = (self.x1*(1-(S_old/self.x1)**2)*tanh(PN/self.x1))/(1+((S_old/self.x1)*tanh(PN/self.x1)))
        ES = (S_old*(2-S_old/self.x1)*tanh(EN/self.x1)) / (1 + ((1-(S_old/self.x1))*tanh(EN/self.x1)))
        S_int = S_old - ES + PS
        perc = S_int*(1-(1+((4*S_int)/(9*self.x1))**4)**(-1/4))
        S_new = S_int - perc
        PR = perc + (PN-PS)
        # compute fluxes based on current state
        F = self.x2 * (R_old/self.x3)**3.5
        Q9 = uh_1_half(vol= 0.9*(PR), d_base=self.x4); # print(Q9.shape)
        Q1 = uh_2_full(vol= 0.1*(PR), d_base=self.x4); # print(Q1)
        # add routed fluxes to timeseries
        if timestep + Q9.shape[0] > self.flux_q9.shape[0]:
            stop = self.flux_q9[timestep:].shape[0]; # print(stop)
            self.flux_q9[timestep:] += Q9[:stop]
        else:
            self.flux_q9[timestep:timestep+Q9.shape[0]] += Q9 # , self.flux_q9[t:t+Q9.shape[0]] )
        if timestep+Q1.shape[0] > self.flux_q1.shape[0]:
            stop = self.flux_q1[timestep:].shape[0]; # print(stop)
            self.flux_q1[timestep:] += Q1[:stop]
        else:
            self.flux_q1[timestep:timestep+Q1.shape[0]] += Q1

        R_int = np.max([0, R_old + self.flux_q9[timestep] + F])
        Qr = R_int*(1-(1+(R_int/self.x3)**4)**(-1/4))
        R_new = R_int - Qr
        # print((R_new - R_old)-(Q9 - Qr + F))
        assert abs((R_new - R_old)-(self.flux_q9[timestep] - Qr + F)) <= 0.00000001, \
        print("t: {}\nqr: {:.04f}\nr_old: {:.04f}\nr_new: {:.04f}\nf: {:.04f}\nq9: {:.04f}\nx2: {:.04f}\nx3: {:.04f}".\
              format(timestep, Qr, R_old, R_new, F, self.flux_q9[timestep], self.x2, self.x3))
        # compute the final outflow
        Qt = Qr + np.max([0, self.flux_q1[timestep] + F])
        # add the final flow to the flux vector
        # add new state to the state vector
        # UPDATE STATES
        self.flux_q[timestep] = Qt # FINAL RUNOFF
        self.flux_es[timestep] = ES
        self.flux_ps[timestep] = PS
        self.flux_perc[timestep] = perc
        self.flux_f[timestep] = F
        self.store_S[timestep] = S_new
        self.store_R[timestep] = R_new
        self.store_S1[timestep] = S1_new
        self.store_S2[timestep] = S2_new
        self.flux_qr[timestep] = Qr
        return

    def simulate(self):
        for t in range(self.dates.shape[0]):
            if t == 0:
                old_states = [0.000, 0, 0, 0]
            else:
                old_states = [self.store_S[t-1],
                              self.store_R[t-1],
                              self.store_S1[t-1],
                              self.store_S2[t-1]]
            T_current = self.temp[t]
            P = self.precip[t]
            self.step_run(old_states=old_states, input_p=P, input_t=T_current, timestep=t)
        self.flux_q = np.where(np.isnan(self.flux_q), 0, self.flux_q)
        return
