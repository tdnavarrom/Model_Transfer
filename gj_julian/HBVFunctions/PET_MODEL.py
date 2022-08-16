#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 05:05:11 2018

@author: chinedum
"""

import numpy as np
import pandas as pd
import datetime

def days_of_year(start_date, end_date):

    # PLEASE INPUT DATES AS YEAR, MONTH, DATE
    # ASSERT START DATE IS GREATER THAN END DATE
    #####################################
    ## PSEUDO CODE
    ## get list of dates between start and end dates
    ## calculate day of year for each date
    #####################################
    dates=pd.date_range(start=start_date, end=end_date, freq='D')
    #print(start_date, end_date, len(dates))
    return np.array([dd.timetuple().tm_yday for dd in dates])

def CBM_hours_of_daylight_model(start_date, end_date, latitude):
    p = 0.833 # constant based on Astronomical calendar
    J = days_of_year(start_date, end_date)
    #print(J)
    theta = 0.2163108 + (2 * np.arctan(0.9671396 * np.arctan(0.00860 * (J - 186))))
    phi = np.arcsin(0.39795 * np.cos(theta))

    numerator = np.sin(p * np.pi / 180) + (np.sin(latitude * np.pi / 180)*np.sin(phi))
    denominator = np.cos(latitude * np.pi / 180) * np.cos(phi)

    hours_of_daylight = 24 - ( (24 / np.pi) * (np.arccos(numerator / denominator)) )
    #print(hours_of_daylight)
    return hours_of_daylight


def PET_Hammon(start_date, end_date, tavg, latitude, hammon_coeff=1.2):
    assert type(tavg) == np.ndarray, 'please input temprature as a numpy ndarray'
    # print(hammon_coeff)
    assert np.any(np.minimum(hammon_coeff, 3.5)), 'check Hammon Coef'
    # assert tavg.shape == latitude.shape, "temp and latitude not same shape"
    # print("latitude shape {}".format(latitude.shape))
    # print("temp shape {}".format(tavg.shape))

    esat = 6.108 * np.exp((17.26939 * tavg)/(tavg + 237.3))
    # print("esat", esat)

    rhosat = 216.7 * esat / (tavg + 237.3)
    # print("rhosat", rhosat)
    # print(start_date, end_date, latitude)
    Ld = CBM_hours_of_daylight_model(start_date, end_date, latitude)[0]
    # print("Ld", Ld)

    PET = 0.1651 * (Ld/12.) * rhosat * hammon_coeff
    # return PET in mm/day
    return PET


def ref_ET(start_date, end_date, tavg, latitude):
    Ld = CBM_hours_of_daylight_model(start_date, end_date, latitude)
    p = Ld/24.0
    ET0 = p*((0.46*tavg) + 8)
    return ET0
