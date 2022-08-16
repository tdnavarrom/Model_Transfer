import pandas as pd
import numpy as np
# IMPORT PYWR MODULE FUNCTIONS
from pywr.parameters import Parameter

class VolumeAreaCurve(Parameter):
    def __init__(self, model, folder, reservoir_name):
        super().__init__(model)
        self._model = model
        self._reservoir_name = reservoir_name
        self._parameters = pd.read_csv("{}/volume_area_parameters.csv".format(folder)).loc[:, self._reservoir_name.lower()]
        print("reservoir parameters", self._parameters)

    def power_function(self, x, a, b, c):
        y = a*x**b + c
        return np.maximum(0, y)

    def value(self, timestep, scenario_index):
        current_storage = self._model.nodes['{}_reservoir'.format(self._reservoir_name)].volume
        # print(current_storage)
        area = self.power_function(x=current_storage, a=self._parameters[0], b=self._parameters[1], c=self._parameters[2])
        # print(area)
        return area

    @classmethod
    def load(cls, model, data):
        """ This is the function that actually reads the
             data defined in the node from which it is called """
        _url = data.pop("url")
        url = _url[0]
        name = _url[1]
        # print(value)
        print("reservoir", name)
        # print(data)
        c = cls(model=model, folder=url, reservoir_name=name) # here we call the
        return c
