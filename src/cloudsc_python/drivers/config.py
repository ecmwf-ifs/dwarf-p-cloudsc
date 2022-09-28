# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from os.path import dirname, join, normpath, splitext
from pydantic import BaseModel, validator
import socket
from typing import Optional

from cloudsc4py.framework.config import DataTypes, GT4PyConfig


class IOConfig(BaseModel):
    output_file: Optional[str]
    host_name: Optional[str]

    @validator("output_file")
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
    def set_host_name(cls, v: Optional[str]) -> str:
        return v or socket.gethostname()

    def with_host_name(self, host_name: str) -> IOConfig:
        args = self.dict()
        args["host_name"] = host_name
        return IOConfig(**args)

    def with_output_file(self, output_file: str) -> IOConfig:
        args = self.dict()
        args["output_file"] = output_file
        return IOConfig(**args)


default_io_config = IOConfig(output_file=None, host_name=None)


class PythonConfig(BaseModel):
    # domain
    nx: Optional[int]

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

    def with_num_runs(self, num_runs: Optional[int]) -> PythonConfig:
        args = self.dict()
        if num_runs is not None:
            args["num_runs"] = num_runs
        return PythonConfig(**args)

    def with_nx(self, nx: Optional[int]) -> PythonConfig:
        args = self.dict()
        if nx is not None:
            args["nx"] = nx
        return PythonConfig(**args)

    def with_validation(self, enabled: bool) -> PythonConfig:
        args = self.dict()
        args["enable_validation"] = enabled
        return PythonConfig(**args)


config_files_dir = normpath(join(dirname(__file__), "../../../config-files"))
default_python_config = PythonConfig(
    nx=None,
    enable_validation=True,
    input_file=join(config_files_dir, "input.h5"),
    reference_file=join(config_files_dir, "reference.h5"),
    num_runs=15,
    data_types=DataTypes(bool=bool, float=np.float64, int=int),
    gt4py_config=GT4PyConfig(backend="numpy", rebuild=False, validate_args=True, verbose=True),
    sympl_enable_checks=True,
)


class FortranConfig(BaseModel):
    mode: str
    num_runs: int
    num_threads: int
    nx: int

    def with_mode(self, mode: str) -> FortranConfig:
        args = self.dict()
        args["mode"] = mode
        return FortranConfig(**args)

    def with_num_runs(self, num_runs: int) -> FortranConfig:
        args = self.dict()
        args["num_runs"] = num_runs
        return FortranConfig(**args)

    def with_num_threads(self, num_threads: int) -> FortranConfig:
        args = self.dict()
        args["num_threads"] = num_threads
        return FortranConfig(**args)

    def with_nx(self, nx: int) -> FortranConfig:
        args = self.dict()
        args["nx"] = nx
        return FortranConfig(**args)


default_fortran_config = FortranConfig(mode="fortran", num_runs=1, num_threads=1, nx=1)
