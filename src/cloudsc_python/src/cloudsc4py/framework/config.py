# -*- coding: utf-8 -*-
from __future__ import annotations
from pydantic import BaseModel
from typing import Any, Dict, Optional, Union, Type


class DataTypes(BaseModel):
    bool: Type
    float: Type
    int: Type


class GT4PyConfig(BaseModel):
    backend: str
    backend_opts: Dict[str, Any] = {}
    build_info: Optional[Dict[str, Any]] = None
    device_sync: bool = True
    dtypes: DataTypes = DataTypes(bool=bool, float=float, int=int)
    exec_info: Optional[Dict[str, Any]] = None
    managed: Union[bool, str] = "gt4py"
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True

    def with_backend(self, backend: Optional[str]) -> GT4PyConfig:
        args = self.dict()
        if backend is not None:
            args["backend"] = backend
        return GT4PyConfig(**args)

    def with_dtypes(self, dtypes: DataTypes) -> GT4PyConfig:
        args = self.dict()
        args["dtypes"] = dtypes
        return GT4PyConfig(**args)

    def with_validate_args(self, flag: bool) -> GT4PyConfig:
        args = self.dict()
        args["validate_args"] = flag
        return GT4PyConfig(**args)
