# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import numpy as np
from typing import Dict, TypeVar, Union

from sympl import DataArray as SymplDataArray

try:
    import cupy as cp
except ImportError:
    cp = np


DataArray = SymplDataArray
DataArrayDict = Dict[str, DataArray]
ParameterDict = Dict[str, Union[bool, float, int]]
Storage = Union[np.ndarray, cp.ndarray]
StorageDict = Dict[str, Storage]
Range = TypeVar("Range")
