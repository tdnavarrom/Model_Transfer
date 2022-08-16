import pandas as pd
import numpy as np
import math

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
    i = np.nanmean(heat_index(np.array(monthly_means)))
    a = (6.75e-7 * i**3) - (7.71e-5 * i**2) + (1.79e-2 * i) + 0.49
    e = 1.6*(10*T/i)**a # e in centimeters
    return np.where(np.isnan(e), 0, e) #*10 # return in millimeters

class HydrologicModel(object):
    print("running sacsma \n")
    def __init__(self, start_date, end_date, timestep, precip, temp, parameters):
        self.start = start_date
        self.end = end_date
        self.timestep = timestep; print(timestep)
        self.precip = precip
        self.temp = temp
        self.parameters = parameters
        self.dates = pd.date_range(start=self.start,
                                   end=self.end,
                                   freq=self.timestep)
        # Capacity Thresholds
        self.uztwm  =  parameters[0] # Upper zone tension water capacity [mm]
        self.uzfwm  =  parameters[1] # Upper zone free water capacity [mm]
        self.lztwm  =  parameters[2] # Lower zone tension water capacity [mm]
        self.lzfpm  =  parameters[3] # Lower zone primary free water capacity [mm]
        self.lzfsm  =  parameters[4] # Lower zone supplementary free water capacity [mm]
        # Recession Parameters
        self.uzk    =  parameters[5] # Upper zone free water lateral depletion rate [1/day]
        self.lzpk   =  parameters[6] # Lower zone primary free water depletion rate [1/day]
        self.lzsk   =  parameters[7] # Lower zone supplementary free water depletion rate [1/day]
        # Percolation
        self.zperc  =  parameters[8] # Percolation demand scale parameter [-]
        self.rexp   =  parameters[9] # Percolation demand shape parameter [-]: exponent of the percolation equation
        self.pfree  =  parameters[10] # Percolating water split parameter (decimal fraction): Fraction of water percolating from upper zone directly to lower zone free water storage
        # Impervious area
        self.pctim  =  parameters[11] # Impervious fraction of the watershed area (decimal fraction)
        self.adimp  =  parameters[12] # Additional impervious areas (decimal fraction)
        # Others
        self.riva   =  0.0 # parameters[13] # Riparian vegetation area (decimal fraction)
        self.side   =  0.0 # parameters[14] # The ratio of deep recharge to channel base flow [-]
        self.rserv  =  0.3 # parameters[15] # Fraction of lower zone free water not transferrable to lower zone tension water (decimal fraction)
        # create storage vectors
        # UPPER ZONE STATES
        self.uztwc_tot = np.zeros(self.dates.shape) #
        self.uzfwc_tot = np.zeros(self.dates.shape) #
        # LOWER ZONES
        self.lztwc_tot = np.zeros(self.dates.shape) #
        self.lzfsc_tot = np.zeros(self.dates.shape) #
        self.lzfpc_tot = np.zeros(self.dates.shape) #
        # ADDITIONAL IMPERVIOUS ZONE
        self.adimc_tot = np.zeros(self.dates.shape)
        # SETUP THE ARRAY OF SOIL STATES
        self.states = np.stack([self.uztwc_tot,
                                self.uzfwc_tot,
                                self.lztwc_tot,
                                self.lzfsc_tot,
                                self.lzfpc_tot,
                                self.adimc_tot], axis=-1)
        self.soil_moisture_fraction = np.zeros(self.dates.shape)
        # create flux vectors
        self.flux_ea   = np.zeros(self.dates.shape) # evaporation
        self.flux_pr   = np.zeros(self.dates.shape)
        self.flux_pet   = np.zeros(self.dates.shape)
        self.surf_tot   = np.zeros(self.dates.shape) # surface
        self.base_tot = np.zeros(self.dates.shape) # baseflow
        self.flux_q = np.zeros(self.dates.shape)

    def snow(self, flux_rf, flux_sf, T_current, WC_state, SP_state, snow_parameters):
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
        """
        THIS FUNCTION WILL
        RUN ONE FORWARD STEP
        OF THE MODEL
        """
        # Initial Storage States
        # SAC-SMA
        uztwc = old_states[0] # Upper zone tension water storage
        uzfwc = old_states[1] # Upper zone free water storage
        lztwc = old_states[2] # Lower zone tension water storage
        lzfsc = old_states[3] # Lower zone supplementary free water storage
        lzfpc = old_states[4] # Upper zone primary free water storage
        adimc = old_states[5] # Additional impervious area storage
        # print("soil states", uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)
        # PERFORMS SIMULATION USING THE SAC-SMA COUPLED WITH SNOW17
        thres_zero = 0.00001 # Threshold to be considered as zero
        parea = 1 - self.adimp - self.pctim
        P = input_p
        T = input_t
        CURRENT_DATE = pd.to_datetime([self.dates[timestep]])
        PET = thornthwaite(T=np.array([T]), dates=CURRENT_DATE)
        edmnd = PET
        # print(PET)
        # partition precip into rain and snow
        pr = self.rainfall(P=P, T_current=T)
        # SF = self.snowfall(P=P, T_current=T)
        # run module
        # ET(1)
        # ET(1), ET from Upper zone tension water storage
        et1 = edmnd * uztwc/self.uztwm
        red = edmnd - et1 # residual ET demand
        uztwc = uztwc - et1
        # ET(2), ET from upper zone free water storage
        et2 = 0
        if uztwc <= 0: # in case et1 > uztws, no water in the upper tension water storage
            et1 = et1 + uztwc # et1 = uztwc
            uztwc = 0
            red = edmnd - et1
            if uzfwc < red: # when upper zone free water content is less than residual ET
                et2 = uzfwc # all content at upper zone free water zone will be gone as ET
                uzfwc = 0
                red = red - et2
                if uztwc < thres_zero: uztwc = 0
                if uzfwc < thres_zero: uzfwc = 0
            else: # when upper zone free water content is larger than residual ET
                et2 = red # all residual ET will be gone as ET
                uzfwc = uzfwc - et2
                red = 0
        else: # all maximum et (et1) are consumed at uztwc, so no et from uzfwc (et2=0)
            # There's possibility that upper zone free water ratio exceeds upper zone tension water ratio
            # If so, free water is transferred to tension water storage
            if (uztwc / self.uztwm) < (uzfwc / self.uzfwm):
                uzrat = (uztwc + uzfwc) / (self.uztwm + self.uzfwm)
                uztwc = self.uztwm * uzrat
                uzfwc = self.uzfwm * uzrat
            if uztwc < thres_zero: uztwc = 0
            if uzfwc < thres_zero: uzfwc = 0
        # ET(3), ET from Lower zone tension water storage when residual ET > 0
        et3 = red * lztwc / (self.uztwm + self.lztwm) # according to this calculation, residual ET is always bigger than ET(3)
        lztwc = lztwc - et3
        if lztwc < 0: # et3 cannot exceed lztws
            et3 = et3 + lztwc # et3 = lztwc
            lztwc = 0
        # ????? why is the red not updated? not updated red is used later for ET(5) calculation
        # Water resupply from Lower free water storages to Lower tension water storage
        saved = self.rserv * (self.lzfpm + self.lzfsm)
        ratlzt = lztwc / self.lztwm
        ratlz = (lztwc + lzfpc + lzfsc - saved) / (self.lztwm + self.lzfpm + self.lzfsm - saved)
        if ratlzt < ratlz: # water is first taken from supplementary water storage for resupply
            del_ = (ratlz - ratlzt) * self.lztwm
            lztwc = lztwc + del_ # Transfer water from lzfss to lztws
            lzfsc = lzfsc - del_
            if lzfsc < 0:  # if tranfer exceeds lzfsc then remainder comes from lzfps
                lzfpc = lzfpc + lzfsc
                lzfsc = 0
        if lztwc < thres_zero: lztwc = 0
        # ET(5), ET from additional impervious (ADIMP) area
        et5 = et1 + (red + et2) * (adimc - et1 - uztwc) / (self.uztwm + self.lztwm) # ????? no idea where this come from, I think there's a possibility that et5 can be negative values
        adimc = adimc - et5
        if adimc < 0: # et5 cannot exceed adims
            et5 = et5 + adimc # et5 = adimc
            adimc = 0
        et5 = et5 * self.adimp
        # Time interval available moisture in excess of uztw requirements
        twx = pr + uztwc - self.uztwm
        if twx < 0: # all moisture held in uztw- no excess
            uztwc = uztwc + pr
            twx = 0
        else: # moisture available in excess of uztw storage
            uztwc = self.uztwm
        # for now twx is excess rainfall after filling the uztwc
        adimc = adimc + pr - twx
        # Compute Impervious Area Runoff
        roimp = pr * self.pctim
        # Initialize time interval sums
        sbf = 0 # Sum of total baseflow(from primary and supplemental storages)
        ssur = 0 # Sum of surface runoff
        sif = 0 # Sum of interflow
        sperc = 0 # Time interval summation of percolation
        sdro = 0 # Sum of direct runoff from the additional impervious area
        # Determine computational time increments for the basic time interval
        ninc = math.floor(1.0 + 0.2*(uzfwc+twx)) # Number of time increments that interval is divided into for further soil-moisture accountng
        dinc = 1 / ninc # Length of each increment in days
        pinc = twx / ninc # Amount of available moisture for each increment
        # Compute free water depletion fractions for the time increment (basic depletions are for one day)
        duz = 1 - (1 - self.uzk)**dinc
        dlzp = 1 - (1 - self.lzpk)**dinc
        dlzs = 1 - (1 - self.lzsk)**dinc
        # Start incremental for-loop for the time interval
        for n in range(1, int(ninc+1)):
            adsur = 0 # Amount of surface runoff. This will be updated.
            # Compute direct runoff from self.adimp area
            ratio = (adimc - uztwc) / self.lztwm
            if ratio < 0: ratio = 0
            addro = pinc*(ratio**2) # Amount of direct runoff from the additional impervious area
            # Compute baseflow and keep track of time interval sum.
            bf_p = lzfpc * dlzp # Baseflow from free water primary storage
            lzfpc = lzfpc - bf_p
            if lzfpc <= 0.0001:
                bf_p = bf_p + lzfpc
                lzfpc = 0
            sbf = sbf + bf_p
            bf_s = lzfsc * dlzs # Baseflow from free water supplemental storage
            lzfsc = lzfsc - bf_s
            if lzfsc <= 0.0001:
                bf_s = bf_s + lzfsc
                lzfsc = 0
            sbf = sbf + bf_s # Total Baseflow from primary and supplemental storages
            # Compute PERCOLATION- if no water available then skip.
            if (pinc + uzfwc) <= 0.01:
                uzfwc = uzfwc + pinc
            else:
                percm = self.lzfpm * dlzp + self.lzfsm * dlzs # Limiting drainage rate from the combined saturated lower zone storages
                perc = percm * uzfwc / self.uzfwm
                defr = 1.0 - (lztwc + lzfpc + lzfsc)/(self.lztwm + self.lzfpm + self.lzfsm) # DEFR is the lower zone moisture deficiency ratio
                if defr < 0: defr = 0
                perc = perc * (1.0 + self.zperc * (defr**self.rexp))
                # Note. . . percolation occurs from uzfws before pav is added
                if perc >= uzfwc: # Percolation rate exceeds uzfws
                    perc = uzfwc
                uzfwc = uzfwc - perc # Percolation rate is less than uzfws.
                check = lztwc + lzfpc + lzfsc + perc - self.lztwm - self.lzfpm - self.lzfsm
                if check > 0: # Check to see if percolation exceeds lower zone deficiency.
                    perc = perc - check
                    uzfwc = uzfwc + check
                sperc = sperc + perc # SPERC is the time interval summation of PERC
                # Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
                del_ = uzfwc * duz # The amount of interflow
                sif = sif + del_
                uzfwc = uzfwc - del_
                # Distribute percolated water into the lower zones. Tension water
                # must be filled first except for the PFREE area. PERCT is
                # percolation to tension water and PERCF is percolation going to
                # free water.
                perct = perc * (1.0 - self.pfree) # Percolation going to the tension water storage
                if (perct + lztwc) <= self.lztwm:
                    lztwc = lztwc + perct
                    percf = 0.0 # Pecolation going to th lower zone free water storages
                else:
                    percf = lztwc + perct - self.lztwm
                    lztwc = self.lztwm
                # Distribute percolation in excess of tension requirements among the free water storages.
                percf = percf + (perc * self.pfree)
                if percf != 0:
                    hpl = self.lzfpm / (self.lzfpm + self.lzfsm) # Relative size of the primary storage as compared with total lower zone free water storages.
                    # Relative fullness of each storage.
                    ratlp = lzfpc / self.lzfpm
                    ratls = lzfsc / self.lzfsm
                    fracp = hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls) # The fraction going to primary
                    if fracp > 1.0: fracp = 1.0
                    percp = percf * fracp # Amount of the excess percolation going to primary
                    percs = percf - percp # Amount of the excess percolation going to supplemental
                    lzfsc = lzfsc + percs
                    if lzfsc > self.lzfsm:
                        percs = percs - lzfsc + self.lzfsm
                        lzfsc = self.lzfsm
                    lzfpc = lzfpc + percf - percs
                    if lzfpc >= self.lzfpm: # Check to make sure lzfps does not exceed self.lzfpm
                        excess = lzfpc - self.lzfpm
                        lztwc = lztwc + excess
                        lzfpc = self.lzfpm
                # Distribute PINC between uzfws and surface runoff
                if pinc != 0:
                    if (pinc + uzfwc) <= self.uzfwm: # check if pinc exceeds self.uzfwm
                        uzfwc = uzfwc + pinc # no surface runoff
                    else:
                        sur = pinc + uzfwc - self.uzfwm # Surface runoff
                        uzfwc = self.uzfwm
                        ssur = ssur + (sur * parea)
                        # ADSUR is the amount of surface runoff which comes from
                        # that portion of self.adimp which is not currently generating
                        # direct runoff. ADDRO/PINC is the fraction of self.adimp
                        # currently generating direct runoff.
                        adsur = sur * (1.0 - addro / pinc)
                        ssur = ssur + adsur * self.adimp
            adimc = adimc + pinc - addro - adsur
            if adimc > (self.uztwm + self.lztwm):
                addro = addro + adimc - (self.uztwm + self.lztwm)
                adimc = self.uztwm + self.lztwm
            sdro  = sdro + (addro * self.adimp) # Direct runoff from the additional impervious area
            if adimc < thres_zero:  adimc = 0
        # Compute sums and adjust runoff amounts by the area over which they are generated.
        # EUSED is the ET from PAREA which is 1.0 - self.adimp - self.pctim
        eused = et1 + et2 + et3
        sif = sif * parea
        # Separate channel component of baseflow from the non-channel component
        tbf = sbf * parea # TBF is the total baseflow
        bfcc = tbf * (1.0 / (1.0 + self.side)) # BFCC is baseflow, channel component
        # Ground flow and Surface flow
        base = bfcc # Baseflow and Interflow are considered as Ground inflow to the channel
        surf = roimp + sdro + ssur + sif # Surface flow consists of Direct runoff and Surface inflow to the channel
        # ET(4)- ET from riparian vegetation.
        et4 = (edmnd - eused) * self.riva # no effect if self.riva is set to zero
        # Check that adims >= uztws
        if adimc < uztwc:
            adimc = uztwc
        # Total inflow to channel for a timestep
        tot_outflow = surf + base - et4
        if tot_outflow < 0:
            et4 = surf + base
            tot_outflow = 0
            surf = 0
            base = 0
        else: # surf and base need to be updated
            surf_remainder = surf - et4
            surf = max(0, surf_remainder)
            if surf_remainder < 0:
                base = base + surf_remainder
        # Compute total evapotransporation - TET
        eused = eused * parea
        tet = eused + et4 + et5
        # UPDATE STATE VECTORS
        self.uztwc_tot[timestep] = uztwc #
        self.uzfwc_tot[timestep] = uzfwc #
        # LOWER ZONES
        self.lztwc_tot[timestep] = lztwc #
        self.lzfsc_tot[timestep] = lzfsc #
        self.lzfpc_tot[timestep] = lzfpc #
        # ADDITIONAL IMPERVIOUS ZONE
        self.adimc_tot[timestep] = adimc
        # UPDATE STATE VECTORS
        self.states[timestep, :] = [uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc]
        self.soil_moisture_fraction[timestep] =  (uztwc+uzfwc) / (self.uztwm + self.uzfwm)
        # UPDATE FLUX VECTORS
        self.flux_pr[timestep]  = pr
        self.flux_pet[timestep] = PET
        self.flux_ea[timestep]  = tet # evaporation
        self.surf_tot[timestep] = surf # surface
        self.base_tot[timestep] = base # baseflow
        self.flux_q[timestep]   = tot_outflow
        return

    def simulate(self):
        for t in range(self.dates.shape[0]):
            if t == 0:
                old_states = np.array([0, 0, 0, 0, 0, 0]) # np.zeros(shape=self.states.shape[1])
            else:
                old_states = self.states[t-1]
#             print(old_states)
            T_current = self.temp[t]
            P = self.precip[t]
            self.step_run(old_states=old_states, input_p=P, input_t=T_current, timestep=t)
        self.flux_q = np.where(np.isnan(self.flux_q), 0, self.flux_q)
        self.flux_q = np.where(self.flux_q < 0, 0, self.flux_q)
        self.flux_ea = np.where(self.flux_ea < 0, 0, self.flux_ea)
        return
