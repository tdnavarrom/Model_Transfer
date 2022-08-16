import pandas as pd
import numpy as np
from scipy import integrate

# hymod routine copied from Bart Nijssen
# https://github.com/bartnijssen/pythonlib/blob/master/hymod.py

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
    T = np.where(T<0, 0, T)
    assert T.shape == dates.shape
    months = dates.month
    monthly_means = []
    for m_ix in np.unique(months):
        ixs = np.in1d(months, m_ix)
        monthly_means.append(np.nanmean(T[ixs]))
    i_sum = np.nanmean(heat_index(np.array(monthly_means)))
    a = (6.75e-7 * i_sum**3) - (7.71e-5 * i_sum**2) + (1.79e-2 * i_sum) + 0.49
    e = 1.6*(10*T/i_sum)**a # e in centimeters
    return np.where(np.isnan(e), 0, e) #*10 # return in millimeters

class HYMOD:
    """Simple 5 parameter model"""

    def __init__(self, start_date, end_date, timestep, precip, temp, parameters):
        """Initiate a hymod instance"""
        self.start = start_date
        self.end = end_date
        self.timestep = timestep
        self.precip = precip
        self.temp = temp
        self.parameters = parameters
        self.dates = pd.date_range(start=self.start,
                                   periods=self.precip.shape[0],
                                   freq=self.timestep)
        self.parameters = parameters
        # assign parameter values
        self.CMAX = parameters[0]
        self.B = parameters[1]
        self.A = parameters[2]
        self.KF = parameters[3]
        self.KS = parameters[4]

        self.SMAX = self.CMAX / (1. + self.B)
        # INITIALIZE STORES
        self.store_SM = np.zeros(self.dates.shape) # soil moisture
        self.store_S = np.zeros(self.dates.shape)
        self.store_S1 = np.zeros(self.precip.shape) # SNOW PACK
        self.store_S2 = np.zeros(self.precip.shape) # LIQUID CONTENT OF SNOW PACK
        # slowflow reservoir
        self.flux_fast = [0, 0, 0] # fastflow reservoirs
        self.error = 0
        self.flux_ea   = np.zeros(self.dates.shape) # evaporation
        self.flux_pet  = np.zeros(self.dates.shape)
        self.flux_q    = np.zeros(self.dates.shape)
        self.flux_ps   = np.zeros(self.dates.shape)
        self.flux_er   = np.zeros(self.dates.shape)
        self.flux_qs   = np.zeros(self.dates.shape)
        self.flux_into_S   = np.zeros(self.dates.shape)

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
        SM_old, S_old, S1_old, S2_old = old_states

        CURRENT_DATE = pd.to_datetime([self.dates[timestep]])
        PET = thornthwaite(T=np.array([input_t]),dates=CURRENT_DATE)
        # partition precip into rain and snow
        RF = self.rainfall(P=input_p, T_current=input_t)
        SF = self.snowfall(P=input_p, T_current=input_t)
        # run snow module
        IN_, S1_new, S2_new, MELT, REFR = self.snow(flux_rf=RF,
                                                   flux_sf=SF,
                                                   T_current=input_t,
                                                   WC_state=S1_old,
                                                   SP_state=S2_old)
#         print(PET)
        IN = IN_
        if SM_old > self.SMAX:
            self.error += SM_old - 0.999 * self.SMAX
            SM_old = 0.999 * self.SMAX
        cprev = self.CMAX * (1 - np.power((1-((self.B+1)*SM_old/self.CMAX)), (1/(self.B+1))))
        ER1 = np.maximum(IN + cprev - self.CMAX, 0.0) # effective rainfal part 1
        IN -= ER1
        dummy = np.minimum(((cprev + IN)/self.CMAX), 1)
        SM_int = (self.CMAX/(self.B+1)) * (1 - np.power((1-dummy), (self.B+1))) # new state
        ER2 = np.maximum(IN-(SM_int-SM_old), 0) # effective rainfall part 2
        evap = np.minimum(SM_int, SM_int/self.SMAX * PET) # actual ET is linearly related to the soil moisture state
        SM_new = SM_int-evap # update state
        UQ = ER1 + self.A * ER2 # quickflow contribution
        US_ = (1 - self.A) * ER2 # slowflow contribution
        for i in range(3):
            self.flux_fast[i] = (1-self.KF) * self.flux_fast[i] + (1-self.KF) * UQ # forecast step
            UQ = (self.KF/(1-self.KF)) * self.flux_fast[i]
        S_new = (1-self.KS) * S_old + (1-self.KS) * US_
        US = (self.KS/(1-self.KS)) * S_new
        Q = UQ + US

        self.flux_q[timestep] = Q # FINAL RUNOFF
        self.store_SM[timestep] = SM_new
        self.store_S[timestep] = S_new
        self.store_S1[timestep] = S1_new
        self.store_S2[timestep] = S2_new
        self.flux_ea[timestep] = evap
        self.flux_pet[timestep] = PET
        self.flux_ps[timestep] = IN_
        self.flux_er[timestep] = ER1+ER2
        self.flux_qs[timestep] = US
        self.flux_into_S[timestep] = US_
        return

    def simulate(self):
        for t in range(self.dates.shape[0]):
            if t == 0:
                old_states = [0.000, 0, 0, 0]
            else:
                old_states = [self.store_SM[t-1],
                              self.store_S[t-1],
                              self.store_S1[t-1],
                              self.store_S2[t-1]]
            # print(old_states)
            T_current = self.temp[t]
            P = self.precip[t]
            self.step_run(old_states=old_states, input_p=P, input_t=T_current, timestep=t)
        self.flux_q = np.where(np.isnan(self.flux_q), 0, self.flux_q)
        return
