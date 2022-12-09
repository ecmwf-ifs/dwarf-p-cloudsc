# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from cloudsc4py.utils.numpyx import assign

if TYPE_CHECKING:
    from sympl._core.data_array import DataArray

    from cloudsc4py.utils.typingx import Storage


def initialize_storage_2d(storage: Storage, buffer: np.ndarray) -> None:
    ni = storage.shape[0]
    mi = buffer.size
    nb = ni // mi
    for b in range(nb):
        assign(storage[b * mi : (b + 1) * mi, 0:1], buffer[:, np.newaxis])
    assign(storage[nb * mi :, 0:1], buffer[: ni - nb * mi, np.newaxis])


def initialize_storage_3d(storage: Storage, buffer: np.ndarray) -> None:
    ni, _, nk = storage.shape
    mi, mk = buffer.shape
    lk = min(nk, mk)
    nb = ni // mi
    for b in range(nb):
        assign(storage[b * mi : (b + 1) * mi, 0:1, :lk], buffer[:, np.newaxis, :lk])
    assign(storage[nb * mi :, 0:1, :lk], buffer[: ni - nb * mi, np.newaxis, :lk])


def initialize_field(field: DataArray, buffer: np.ndarray) -> None:
    if field.ndim == 2:
        initialize_storage_2d(field.data, buffer)
    elif field.ndim == 3:
        initialize_storage_3d(field.data, buffer)
    else:
        raise ValueError("The field to initialize must be either 2-d or 3-d.")
