# IMPORT OTHER FUNCTIONS
import pandas as pd
import numpy as np
import datetime
import os
import json
from pathlib import Path

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


def recursive_search(folder_path:Path, project:str):
    for file in folder_path.iterdir():
        #print(file.name)
        json_creator(file.name, project)


def json_creator(file_name:str, project:str):

    preffix_name = file_name.split(".")[0]

    with open("solis_json_files/template/SOLIS.json") as file:
        json_temp = json.load(file)
    # print(json_temp.keys())
    # print(json_temp["parameters"].keys())
    # print(json_temp["parameters"]["flow_solis_sim"])
    # print(json_temp["parameters"]["flow_solis_sim"]["url"].split("/"))
    temp_file = json_temp["parameters"]["flow_solis_sim"]["url"].split("/")
    temp_file[-1] = file_name
    new_file_name = "/".join(temp_file)
    print("File generated: " + new_file_name)
    json_temp["parameters"]["flow_solis_sim"]["url"] = new_file_name
    
    file_path = "solis_json_files/"+preffix_name+".json"

    with open(file_path, "w") as outfile:
        outfile.write(json.dumps(json_temp,indent=4))


def run_model(model_json_path:Path):

    for file in model_json_path.iterdir():
        
        if not file.is_dir():
            
            preffix = file.name.split(".")[0]
            print("Executing model with: " + file.name)

            file_path_str = str(model_json_path/file.name)

            model = Model.load(file_path_str)

            #model = Model.load("./solis_json_files/solis_inflow2.json")
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
            results.to_csv("../Output/"+preffix+"_test.csv")



input_folder_path = "Input"
project = "Solis"
json_files_input = "solis_json_files"

main_path = Path().resolve()
print(main_path)

flow_solis_sim = "hydrology_results"
flow_solis_sim_path = main_path.parent / input_folder_path / project / flow_solis_sim
recursive_search(flow_solis_sim_path, project)

model_json_path = main_path / json_files_input
run_model(model_json_path)

