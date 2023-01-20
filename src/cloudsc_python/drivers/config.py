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
from os.path import dirname, join, normpath, splitext
from pydantic import BaseModel, validator
import socket
from typing import Optional

from cloudsc4py.framework.config import DataTypes, GT4PyConfig


class IOConfig(BaseModel):
    """Gathers options for I/O."""

    output_csv_file: Optional[str]
    host_name: Optional[str]

    @validator("output_csv_file")
    @classmethod
    def check_extension(cls, v: Optional[str]) -> Optional[str]:
        if v is None:
            return v

        basename, extension = splitext(v)
        if extension == "":
            return v + ".csv"
        elif extension == ".csv":
            return v
        else:
            return basename + ".csv"

    @validator("host_name")
    @classmethod
    def set_host_name(cls, v: Optional[str]) -> str:
        return v or socket.gethostname()

    def with_host_name(self, host_name: str) -> IOConfig:
        args = self.dict()
        args["host_name"] = host_name
        return IOConfig(**args)

    def with_output_csv_file(self, output_csv_file: str) -> IOConfig:
        args = self.dict()
        args["output_csv_file"] = output_csv_file
        return IOConfig(**args)


default_io_config = IOConfig(output_file=None, host_name=None)


class PythonConfig(BaseModel):
    """Gathers options controlling execution of Python/GT4Py code."""

    # domain
    num_cols: Optional[int]

    # validation
    enable_validation: bool
    input_file: str
    reference_file: str

    # run
    num_runs: int

    # low-level and/or backend-related
    data_types: DataTypes
    gt4py_config: GT4PyConfig
    sympl_enable_checks: bool

    @validator("gt4py_config")
    @classmethod
    def add_dtypes(cls, v, values) -> GT4PyConfig:
        return v.with_dtypes(values["data_types"])

    def with_backend(self, backend: Optional[str]) -> PythonConfig:
        args = self.dict()
        args["gt4py_config"] = GT4PyConfig(**args["gt4py_config"]).with_backend(backend).dict()
        return PythonConfig(**args)

    def with_checks(self, enabled: bool) -> PythonConfig:
        args = self.dict()
        args["gt4py_config"] = (
            GT4PyConfig(**args["gt4py_config"]).with_validate_args(enabled).dict()
        )
        args["sympl_enable_checks"] = enabled
        return PythonConfig(**args)

    def with_num_cols(self, num_cols: Optional[int]) -> PythonConfig:
        args = self.dict()
        if num_cols is not None:
            args["num_cols"] = num_cols
        return PythonConfig(**args)

    def with_num_runs(self, num_runs: Optional[int]) -> PythonConfig:
        args = self.dict()
        if num_runs is not None:
            args["num_runs"] = num_runs
        return PythonConfig(**args)

    def with_validation(self, enabled: bool) -> PythonConfig:
        args = self.dict()
        args["enable_validation"] = enabled
        return PythonConfig(**args)


config_files_dir = normpath(join(dirname(__file__), "../../../config-files"))
default_python_config = PythonConfig(
    num_cols=1,
    enable_validation=True,
    input_file=join(config_files_dir, "input.h5"),
    reference_file=join(config_files_dir, "reference.h5"),
    num_runs=15,
    data_types=DataTypes(bool=bool, float=np.float64, int=int),
    gt4py_config=GT4PyConfig(backend="numpy", rebuild=False, validate_args=True, verbose=True),
    sympl_enable_checks=True,
)


class FortranConfig(BaseModel):
    """Gathers options controlling execution of FORTRAN code."""

    build_dir: str
    variant: str
    nproma: int
    num_cols: int
    num_runs: int
    num_threads: int

    def with_build_dir(self, build_dir: str) -> FortranConfig:
        args = self.dict()
        args["build_dir"] = build_dir
        return FortranConfig(**args)

    def with_nproma(self, nproma: int) -> FortranConfig:
        args = self.dict()
        args["nproma"] = nproma
        return FortranConfig(**args)

    def with_num_cols(self, num_cols: int) -> FortranConfig:
        args = self.dict()
        args["num_cols"] = num_cols
        return FortranConfig(**args)

    def with_num_runs(self, num_runs: int) -> FortranConfig:
        args = self.dict()
        args["num_runs"] = num_runs
        return FortranConfig(**args)

    def with_num_threads(self, num_threads: int) -> FortranConfig:
        args = self.dict()
        args["num_threads"] = num_threads
        return FortranConfig(**args)

    def with_variant(self, variant: str) -> FortranConfig:
        args = self.dict()
        args["variant"] = variant
        return FortranConfig(**args)


default_fortran_config = FortranConfig(
    build_dir=".", variant="fortran", nproma=32, num_cols=1, num_runs=1, num_threads=1
)
