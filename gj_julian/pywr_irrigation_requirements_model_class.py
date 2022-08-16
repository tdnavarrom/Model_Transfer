# IMPORT OTHER FUNCTIONS
import pandas as pd
import numpy as np
# IMPORT PYWR MODULE FUNCTIONS
from pywr.parameters import Parameter
import HBVFunctions.PET_MODEL as PET_MODEL

class CropCoefficients(Parameter):
    """
    This will read a list of crops and their
    constituent fractions and return a
    weighted sum of the monthly crop coefficients
    """
    def __init__(self, model, title, cropdata, name):
        super().__init__(model)
        # crop_coeff = pd.read_csv
        # TODO READ FROM FILE THAT IS DERIVED FROM THE FAO-DATABASE
        print("here")
        self._model = model
        self._name = name
        # print(cropdata)
        self._cropdata = cropdata
        self._start_date = self._model.timestepper.start
        self._end_date = self._model.timestepper.end
        self._freq = model.timestepper.freq
        self._dates = pd.date_range(start=self._start_date, end=self._end_date, freq=self._freq[-1])
        self.fao_data = pd.read_csv("../Data/FAO_KC_DATA.csv", index_col='Crop')

    def kc_timeseries(self, crop, start_month, current_date):
        dates = pd.date_range(start=self._start_date, end=self._end_date, freq='1D')
        start_month_ix = np.where(dates.month == start_month)[0]
        year_ix = start_month_ix[:-1][np.diff(start_month_ix)>300]
        kci, kc1, kc2, kce, di, d1, d2, de = self.fao_data.loc[crop, ['kc_init', 'kc_m1', 'kc_m2', 'kc_end', 'd_init', 'd_m1', 'd_m2', 'd_end']]
        kcs = [kci, kc1, kc2, kce]
        ds = [di, d1, d2, de]
        kc_array = np.zeros(dates.shape[0])
        for m in year_ix:
            d = m
            for i in range(len(kcs)):
                kc = kcs[i]
                d_ = ds[i]
                # print(kc, d, d_+d)
                kc_array[int(d):int(d_+d)] = kc
                d += d_
        return pd.Series(kc_array, index=dates).loc[current_date]

    def value(self, timestep, scenario_index):
        # print("value function")
        # # print("cropdata", self._cropdata)
        day = timestep.day
        month = timestep.month
        year = timestep.year
        unit = self._model.timestepper.freq
        current_date = pd.to_datetime("{:04}-{:02}-{:02}".format(year, month, day))
        crops = list(self._cropdata.keys())
        weights = list(self._cropdata.values())
        # print(weights)
        # kc_list = np.array([ self._crop_coeff[crop] for crop in crops ])
        # weight_list = np.array([ self._cropdata[crop] for crop in crops ])
        # kc = np.average(a=kc_list, axis=0, weights=weight_list)
        kc_list = []
        for crop in crops:
            start_month = self.fao_data.loc[crop, 'start_month']
            _kc = self.kc_timeseries(crop=crop, start_month=start_month, current_date=current_date)
            kc_list.append(_kc)
        # print("kc", kc)
        if np.all(np.array(kc_list) == 0):
            kc_array = np.array(kc_list)
            weight_array = np.array(weights)
        else:
            kc_array = np.array(kc_list)[np.array(kc_list) > 0]
            weight_array = np.array(weights)[np.array(kc_list) > 0]
        kc = np.average(a=kc_array, axis=0, weights=weight_array)
        return kc

    @classmethod
    def load(cls, model, data):
        print("running the CropCoefficients class")
        factor = data.pop("factor")
        cropdata = data.pop("cropdata")
        name = data.pop("url")
        print("crop coeff cropdata", cropdata)
        title = model.metadata['title']
        c = cls(model, title=title, cropdata=cropdata, name=name[1])
        return c

class PETDemand(Parameter):
    """
    This will compute the irrigation requirements
    as the potential evaporation
    """
    def __init__(self, model, folder, title, scheme_name, hydrologic_model_type):
        print("running the IrrigationRequirements class")
        super().__init__(model)
        # self._value = factor
        self._model_type = hydrologic_model_type
        self._model = model
        self._name = scheme_name
        self._title = title
        self._folder = folder
        self._input_data_folder = "{}/{}/catchments/{}".format(folder, title, scheme_name)
        # READ THE TEMP DATA
        self._temp_data = pd.read_csv("{}/temp.csv".format(self._input_data_folder),
                                                header=[0], parse_dates=True, index_col=[0])
        # READ THE PARAMETERS
        self._parameters = pd.read_csv("{}/{}_parameters.csv".format(self._input_data_folder, self._model_type),
                                                header=None)
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
        current_date = pd.to_datetime("{:04}-{:02}-{:02}".format(year, month, day))
        # print("current date", current_date)
        # print("timedelta value: {}\ntimedelta unit: {}".format(unit[0], unit[1]))
        end_date = current_date + pd.Timedelta(value=int(unit[0]), unit=str(unit[1]))
        return end_date

    def value(self, timestep, scenario_index):
        lat = self._hru_data.north.values[0]
        current_date = "{:04}-{:02}-{:02}".format(timestep.year, timestep.month, timestep.day)
        next_date = self.get_next_date(timestep, self._model)
        current_temp = self._temp_data.loc[current_date].values[0]
        hammon_coeff = self._parameters.values[0, :][0]
        result = PET_MODEL.PET_Hammon(start_date=current_date,
                                        end_date=next_date,
                                        tavg=np.array(current_temp),
                                        latitude=lat,
                                        hammon_coeff=hammon_coeff)
        self._area = self._hru_data.area.values[0]
        # area is in sq meters
        # pet is returned in mm/day, we convert this to m3/second
        val = (result[0]*self._area/(1000*24*3600))
        # print("running value function", val)
        return val

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
        factor = data.pop("factor")
        # print(value)
        print("irrigation catchment name", name)
        # print(data)
        c = cls(model, url, title, name, hydrologic_model_type) # here we call the
        return c
