# IMPORT OTHER FUNCTIONS
import pandas as pd
import numpy as np
# IMPORT PYWR MODULE FUNCTIONS
from pywr.parameters import Parameter
import HBVFunctions.PET_MODEL as PET_MODEL

def heat_index(T):
    i  = (T/5)**1.514
    return i

def thornthwaite(T, dates):
    # print(T)
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

class ReservoirPETDemand(Parameter):
    """
    This will compute the irrigation requirements
    as the potential evaporation
    """
    def __init__(self, model, folder, title, scheme_name, hydrologic_model_type, pet_factor):
        print("running the IrrigationRequirements class")
        super().__init__(model)
        # self._value = factor
        self._pet_factor = pet_factor
        self._model_type = hydrologic_model_type
        self._model = model
        self._name = scheme_name
        self._title = title
        self._folder = folder
        self._input_data_folder = "{}/{}/reservoirs/{}".format(folder, title, scheme_name)
        print(self._input_data_folder)
        # READ THE TEMP DATA
        self._temp_data = pd.read_csv("{}/temp.csv".format(self._input_data_folder),
                                                header=[0], parse_dates=True, index_col=[0])
        # READ THE HRU DATA
        self._hru_data = pd.read_csv("{}/hru_info.csv".format(self._input_data_folder),
                                                header=[0])

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
        lat = self._hru_data.north.values[0]
        current_date = "{:04}-{:02}-{:02}".format(timestep.year, timestep.month, timestep.day)
        next_date = self.get_next_date(timestep, self._model)
        current_temp = self._temp_data.loc[current_date].values[0]
        T = current_temp
        # CURRENT_DATE = pd.to_datetime([current_date])
        PET = PET_MODEL.PET_Hammon(start_date=current_date,
                                        end_date=next_date,
                                        tavg=np.array(current_temp),
                                        latitude=lat,
                                        hammon_coeff=1)
        # PET = thornthwaite(T=np.array([T]), dates=CURRENT_DATE)
        # pet is returned in mm/day
        return PET*self._pet_factor

    @classmethod
    def load(cls, model, data):
        """ This is the function that actually reads the
             data defined in the node from which it is called """
        _url = data.pop("url")
        url = _url[0]
        name = _url[1]
        title = model.metadata['title']
        hydrologic_model_type = model.metadata["hydrologic_model_type"]
        print("irrigation title", title)
        factor = data.pop("factor"); print("reservoir pet factor", factor)
        # print(value)
        print("irrigation catchment name", name)
        # print(data)
        c = cls(model, url, title, name, hydrologic_model_type, factor) # here we call the
        return c
