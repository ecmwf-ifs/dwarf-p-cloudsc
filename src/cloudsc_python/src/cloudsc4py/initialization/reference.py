# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
from datetime import datetime
from functools import partial
from typing import TYPE_CHECKING

from cloudsc4py.framework.grid import I, J, K
from cloudsc4py.framework.storage import allocate_data_array
from cloudsc4py.initialization.utils import initialize_field

if TYPE_CHECKING:
    from typing import Literal, Tuple

    from sympl._core.data_array import DataArray
    from sympl._core.typingx import DataArrayDict

    from cloudsc4py.framework.config import GT4PyConfig
    from cloudsc4py.framework.grid import ComputationalGrid, DimSymbol
    from cloudsc4py.utils.iox import HDF5Reader


def allocate_tendencies(
    computational_grid: ComputationalGrid, *, gt4py_config: GT4PyConfig
) -> DataArrayDict:
    def allocate(units: str = "") -> DataArray:
        return allocate_data_array(
            computational_grid, (I, J, K), units, gt4py_config=gt4py_config, dtype="float"
        )

    return {
        "time": datetime(year=2022, month=1, day=1),
        "f_a": allocate(),
        "f_qi": allocate(),
        "f_ql": allocate(),
        "f_qr": allocate(),
        "f_qs": allocate(),
        "f_qv": allocate(),
        "f_t": allocate(),
    }


def initialize_tendencies(tendencies: DataArrayDict, hdf5_reader: HDF5Reader) -> None:
    hdf5_reader_keys = {"f_a": "TENDENCY_LOC_A", "f_qv": "TENDENCY_LOC_Q", "f_t": "TENDENCY_LOC_T"}
    for name, hdf5_reader_key in hdf5_reader_keys.items():
        buffer = hdf5_reader.get_field(hdf5_reader_key)
        initialize_field(tendencies[name], buffer)

    cld = hdf5_reader.get_field("TENDENCY_LOC_CLD")
    for idx, name in enumerate(("f_ql", "f_qi", "f_qr", "f_qs")):
        initialize_field(tendencies[name], cld[..., idx])


def allocate_diagnostics(
    computational_grid: ComputationalGrid, *, gt4py_config: GT4PyConfig
) -> DataArrayDict:
    def _allocate(
        grid_id: Tuple[DimSymbol, ...], units: str, dtype: Literal["bool", "float", "int"]
    ) -> DataArray:
        return allocate_data_array(
            computational_grid, grid_id, units, gt4py_config=gt4py_config, dtype=dtype
        )

    allocate = partial(_allocate, grid_id=(I, J, K), units="", dtype="float")
    allocate_h = partial(_allocate, grid_id=(I, J, K - 1 / 2), units="", dtype="float")
    allocate_ij = partial(_allocate, grid_id=(I, J), units="", dtype="float")

    return {
        "time": datetime(year=2022, month=1, day=1),
        "f_covptot": allocate(),
        "f_fcqlng": allocate_h(),
        "f_fcqnng": allocate_h(),
        "f_fcqrng": allocate_h(),
        "f_fcqsng": allocate_h(),
        "f_fhpsl": allocate_h(),
        "f_fhpsn": allocate_h(),
        "f_fplsl": allocate_h(),
        "f_fplsn": allocate_h(),
        "f_fsqif": allocate_h(),
        "f_fsqitur": allocate_h(),
        "f_fsqlf": allocate_h(),
        "f_fsqltur": allocate_h(),
        "f_fsqrf": allocate_h(),
        "f_fsqsf": allocate_h(),
        "f_rainfrac_toprfz": allocate_ij(),
    }


def initialize_diagnostics(diagnostics: DataArrayDict, hdf5_reader: HDF5Reader) -> None:
    hdf5_reader_keys = {name: "P" + name[2:].upper() for name in diagnostics if name != "time"}
    for name, hdf5_reader_key in hdf5_reader_keys.items():
        buffer = hdf5_reader.get_field(hdf5_reader_key)
        initialize_field(diagnostics[name], buffer)


def get_reference_tendencies(
    computational_grid: ComputationalGrid, hdf5_reader: HDF5Reader, *, gt4py_config: GT4PyConfig
) -> DataArrayDict:
    tendencies = allocate_tendencies(computational_grid, gt4py_config=gt4py_config)
    initialize_tendencies(tendencies, hdf5_reader)
    return tendencies


def get_reference_diagnostics(
    computational_grid: ComputationalGrid, hdf5_reader: HDF5Reader, *, gt4py_config: GT4PyConfig
) -> DataArrayDict:
    diagnostics = allocate_diagnostics(computational_grid, gt4py_config=gt4py_config)
    initialize_diagnostics(diagnostics, hdf5_reader)
    return diagnostics
