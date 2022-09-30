# -*- coding: utf-8 -*-
from typing import Dict, TypeVar, Union

from sympl import DataArray as SymplDataArray

from gt4py.storage import Storage


Array = Storage
ArrayDict = Dict[str, Storage]
DataArray = SymplDataArray
DataArrayDict = Dict[str, DataArray]
ParameterDict = Dict[str, Union[bool, float, int]]
StorageDict = dict[str, Storage]
Range = TypeVar("Range")
