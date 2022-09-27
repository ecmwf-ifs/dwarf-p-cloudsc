# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from os.path import dirname, join, normpath
from pydantic import BaseModel, validator
from typing import Optional

from cloudsc4py.framework.config import DataTypes, GT4PyConfig


class Config(BaseModel):
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

    def with_backend(self, backend: Optional[str]) -> Config:
        args = self.dict()
        args["gt4py_config"] = GT4PyConfig(**args["gt4py_config"]).with_backend(backend).dict()
        return Config(**args)

    def with_checks(self, enabled: bool) -> Config:
        args = self.dict()
        args["gt4py_config"] = (
            GT4PyConfig(**args["gt4py_config"]).with_validate_args(enabled).dict()
        )
        args["sympl_enable_checks"] = enabled
        return Config(**args)

    def with_num_runs(self, num_runs: Optional[int]) -> Config:
        args = self.dict()
        if num_runs is not None:
            args["num_runs"] = num_runs
        return Config(**args)

    def with_nx(self, nx: Optional[int]) -> Config:
        args = self.dict()
        if nx is not None:
            args["nx"] = nx
        return Config(**args)

    def with_validation(self, enabled: bool) -> Config:
        args = self.dict()
        args["enable_validation"] = enabled
        return Config(**args)


config_files_dir = normpath(join(dirname(__file__), "../../../config-files"))
default_config = Config(
    nx=None,
    enable_validation=True,
    input_file=join(config_files_dir, "input.h5"),
    reference_file=join(config_files_dir, "reference.h5"),
    num_runs=15,
    data_types=DataTypes(bool=bool, float=np.float64, int=int),
    gt4py_config=GT4PyConfig(backend="numpy", rebuild=False, validate_args=True, verbose=True),
    sympl_enable_checks=True,
)
