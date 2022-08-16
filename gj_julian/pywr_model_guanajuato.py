# IMPORT OTHER FUNCTIONS
import pandas as pd
import numpy as np
import datetime
import os

# IMPORT PYWR MODULE FUNCTIONS
from pywr.core import Model

# THESE ARE THE FUNCTIONS THAT WILL COMPUTE THE HBV SOIL RESERVOIRS
from pywr_hydrologic_model_class import HydrologicModelClass, GetHydrologicModelValue
HydrologicModelClass.register()
GetHydrologicModelValue.register()
# THESE ARE THE FUNCTIONS THAT WILL COMPUTE THE IRRIGAITON REQUIREMENTS,
# AND CROP WATER USE
from pywr_irrigation_requirements_model_class import CropCoefficients, PETDemand
CropCoefficients.register()
PETDemand.register()
# FUNCTIONS FOR RESERVOIRS PET
from pywr_reservoir_pet_model_class import ReservoirPETDemand
ReservoirPETDemand.register()
# FUNCTIONS FOR RESERVOIR VOLUME TO AREA
from pywr_reservoir_volume_area_class import VolumeAreaCurve
VolumeAreaCurve.register()
# FUNCTIONS FOR RETURN FLOWS

# model = Model.load("./ToyModelJSON/gj_trial_observed.json")
model = Model.load("./gj_json_files/gj_outlet_section.json")

# run, forest, run!
start_date = model.timestepper.start
end_date = model.timestepper.end
freq = model.timestepper.freq
dates = pd.date_range(start=start_date, end=end_date, freq=freq)
print("dates shape", dates.shape)

ST = datetime.datetime.now()
for i in range(dates.shape[0]):
    model.step()
TT = datetime.datetime.now() - ST
print(TT)
results = model.to_dataframe()
results.columns = results.columns.get_level_values(0)
results.to_csv("../Output/gj_output_test_1014.csv")
