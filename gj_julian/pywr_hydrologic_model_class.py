# IMPORT OTHER FUNCTIONS
import pandas as pd
import numpy as np
# IMPORT PYWR MODULE FUNCTIONS
from pywr.parameters import Parameter

# SETUP THE CLASS FOR THE PARAMETER
class HydrologicModelClass(Parameter):
    def __init__(self, model, folder, title, catchment_name, factor, irr, hydrologic_model_type):
        super().__init__(model)
        # import the right hydrologic model
        self._model_type = hydrologic_model_type
        self._value = factor
        self._model = model
        self._name = catchment_name
        self._title = title
        self._folder = folder
        self._irr = irr
        self._unit = model.timestepper.freq
        self._input_data_folder = "{}/{}/catchments/{}".format(folder, title, catchment_name)

        # READ THE PRECIP DATA
        self._precip_data = pd.read_csv("{}/precip.csv".format(self._input_data_folder),
                                        header=[0], parse_dates=True, index_col=[0]).resample("{}".format(self._unit[1])).mean()
        # READ THE TEMP DATA
        self._temp_data = pd.read_csv("{}/temp.csv".format(self._input_data_folder),
                                        header=[0], parse_dates=True, index_col=[0]).resample("{}".format(self._unit[1])).mean()
        # READ THE PARAMETERS
        if self._unit[1] == 'M':
            self._parameters = pd.read_csv("{}/{}_parameters_monthly.csv".format(self._input_data_folder, self._model_type),
                                                header=None)
        else:
            self._parameters = pd.read_csv("{}/{}_parameters.csv".format(self._input_data_folder, self._model_type),
                                                header=None)
        # READ THE HRU DATA
        self._hru_data = pd.read_csv("{}/hru_info.csv".format(self._input_data_folder),
                                                header=[0])
        # GET THE BASIN AREA SQUARE METERS
        self._hru_area = self._hru_data.area.values[0]

        print("parameter_set", self._parameters.values[0])
        # INSTANTIATE THE HYDROLOGIC MODEL
        self._timestep = model.timestepper.freq

        start_date = model.timestepper.start.date()
        end_date = model.timestepper.end.date()
        HydrologicModel = self.import_hydrologic_model(self._model_type)
        self._hydrologic_model_instance = HydrologicModel.HydrologicModel(start_date=start_date,
                                                            end_date=end_date,
                                                            timestep=self._timestep,
                                                            precip=self._precip_data,
                                                            temp=self._temp_data,
                                                            parameters=self._parameters.values[0])

    def import_hydrologic_model(self, model_type):
        if model_type == "SACSMA":
            import hydrologic_models.SACSMA as HydrologicModel
        elif model_type == "HBV9P":
            import hydrologic_models.HBV9P_96 as HydrologicModel
        else:
            print("ERROR: PLEASE SELECT HBV, OR SACSMA")
            return
        return HydrologicModel

    def get_next_date(self, timestep, model):
        """ This function will read the model timestepper
            and will
        """
        # print("input model", model)
        day = timestep.day
        month = timestep.month
        year = timestep.year
        unit = model.timestepper.freq
        current_date = np.datetime64("{:04}-{:02}-{:02}".format(year, month, day), str(unit[1]))
        # print("current date", current_date)
        # print("timedelta value: {}\ntimedelta unit: {}".format(unit[0], unit[1]))
        end_date = current_date + np.timedelta64(int(unit[0]), str(unit[1]))
        return end_date

    def value(self, timestep, scenario_index):
        # run the model here
        self.run_model(timestep)
        return 0.0

    def run_model(self, timestep):
        # print("ix", timestep.index)
        current_date = "{:04}-{:02}-{:02}".format(timestep.year, timestep.month, timestep.day)
        next_date = self.get_next_date(timestep, self._model)
        # print(scenario_index.global_id)
        # print(self._value)
        current_precip = self._precip_data.loc[current_date].values[0]
        current_temp = self._temp_data.loc[current_date].values[0]
        # print(current_precip, current_temp, "\n")
        parameter_set = self._parameters.values[0, :]
        # get_state = #
        ix = timestep.index
        if ix == 0:
            soil_states = np.zeros(self._hydrologic_model_instance.states.shape[1])
        else:
            soil_states = self._hydrologic_model_instance.states[ix-1]
        if self._irr == "True":
            # search the demand node for the water that is delivered
            irr_delivery_parameter = "{}_demand".format(self._name)
            deliveries = self._model.nodes[irr_delivery_parameter].prev_flow[0]
            deliveries_mm = deliveries * (1000) / self._hru_area
            # print("irrigated deliveries", deliveries_mm)
        else:
            deliveries_mm = 0
        water_input = current_precip+deliveries_mm
        # RUN THE HYDROLOGIC MODEL FOR THE CURRENT TIMESTEP
        # print("stepping through")
        # print("current_precip", current_precip, "\ncurrent_temp", current_temp, "\ncurrent_date", current_date)
        self._hydrologic_model_instance.step_run(old_states=soil_states, input_p=water_input, input_t=current_temp, timestep=ix)
        return

    def get_url(self):
        return self._folder

    def get_pet_demand(self, timestep):
        pet_demand = self._hydrologic_model_instance.flux_pet[timestep.index-1]
        # sm_fraction = self._hydrologic_model_instance.soil_moisture_fraction[timestep.index-1]
        pet_m3 = pet_demand * self._hru_area / 1e3
        return pet_m3

    def get_soilmoisture_fraction(self, timestep):
        soil_moist_frac = self._hydrologic_model_instance.soil_moisture_fraction[timestep.index-1]
        # print("soil moisture fraction", 1-soil_moist_frac)
        return np.max([1 - soil_moist_frac, 0])

    def get_discharge(self, timestep):
        # get basin runoff in millimeters
        discharge_mm = self._hydrologic_model_instance.flux_q[timestep.index-1]
        # convert to cubic meters
        discharge_m3 = discharge_mm * self._hru_area / 1e3
        # print("\n", self._name, timestep.index, "\ndischarge", discharge_m3,"\nrunoff", discharge_mm, "\narea", self._hru_area)
        return discharge_m3

    def get_upperstate(self, timestep):
        uztwc = self._hydrologic_model_instance.states[timestep.index-1, 0]
        uzfwc = self._hydrologic_model_instance.states[timestep.index-1, 1]
        return (uztwc + uzfwc) * self._hru_area / 1e3

    def get_lowerstate(self, timestep):
        lztwc = self._hydrologic_model_instance.states[timestep.index-1, 2]
        lzfsc = self._hydrologic_model_instance.states[timestep.index-1, 3]
        return (lztwc + lzfsc ) * self._hru_area / 1e3

    def get_evapotranspiration(self, timestep):
        # get basin et in millimeters
        et_mm = self._hydrologic_model_instance.flux_ea[timestep.index-1]
        # convert to cubic meters
        et_m3 = et_mm * self._hru_area / 1e3
        return et_m3

    def get_precip(self, timestep):
        # get basin et in millimeters
        pr_mm = self._hydrologic_model_instance.flux_pr[timestep.index-1]
        # convert to cubic meters
        pr_m3 = pr_mm * self._hru_area / 1e3
        return pr_m3

class GetHydrologicModelValue(Parameter):
    def __init__(self, model, hydrological_model_parameter, val):
        super().__init__(model)
        self._val = val
        self.hydrological_model_parameter = hydrological_model_parameter
        # self.hydrological_model_parameter.run_model()
        # self.parents.add(hsm_frac_irr_dist_45ydrological_model_parameter)

    def value(self, timestep, scenario_index):
        if self._val == 'flow':
            return self.hydrological_model_parameter.get_discharge(timestep)
        elif self._val == 'et':
            return self.hydrological_model_parameter.get_evapotranspiration(timestep)
        elif self._val == 'pet_demand':
            return self.hydrological_model_parameter.get_pet_demand(timestep)
        elif self._val == 'sm_fraction':
            return self.hydrological_model_parameter.get_soilmoisture_fraction(timestep)
        elif self._val == 'state_uz':
            return self.hydrological_model_parameter.get_upperstate(timestep)
        elif self._val == 'precip':
            return self.hydrological_model_parameter.get_precip(timestep)
        else:
            print("ERROR")

    @classmethod
    def load(cls, model, data):
        print("run the HydrologicModelClass")
        """ This is the function that actually reads the
             date defined in the node from which
             it is called """
        url = data["url"]["url"]
        name = data["url"]["name"]
        irr = data["url"]["irr"]
        model_type = model.metadata["hydrologic_model_type"]
        title = model.metadata["title"]
        val = data["url"]["value"]
        print("title", title, val)
        factor = data.pop("factor")
        print("catchment name", name)
        # call the hydrologic model with the data passed in this class
        hydrologic_model = HydrologicModelClass(model, url, title, name, factor, irr, model_type)
        print(hydrologic_model.get_url())
        c = cls(model, hydrologic_model, val)
        return c
