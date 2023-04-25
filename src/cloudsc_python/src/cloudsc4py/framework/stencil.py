# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
from typing import TYPE_CHECKING

from gt4py.cartesian import gtscript

if TYPE_CHECKING:
    from typing import Any, Dict

    from gt4py.cartesian import StencilObject

    from cloudsc4py.framework.config import GT4PyConfig


FUNCTION_COLLECTION = {}
STENCIL_COLLECTION = {}


def function_collection(name: str):
    """Decorator for GT4Py functions."""
    if name in FUNCTION_COLLECTION:
        raise RuntimeError(f"Another function called `{name}` found.")

    def core(definition):
        FUNCTION_COLLECTION[name] = {"definition": definition}
        return definition

    return core


def stencil_collection(name: str):
    """Decorator for GT4Py stencil definitions."""
    if name in STENCIL_COLLECTION:
        raise RuntimeError(f"Another stencil called `{name}` found.")

    def core(definition):
        STENCIL_COLLECTION[name] = {"definition": definition}
        return definition

    return core


def compile_stencil(
    name: str,
    gt4py_config: GT4PyConfig,
    externals: Dict[str, Any] = None,
) -> StencilObject:
    """Automate and customize the compilation of GT4Py stencils."""
    stencil_info = STENCIL_COLLECTION.get(name, None)
    if stencil_info is None:
        raise RuntimeError(f"Unknown stencil `{name}`.")
    definition = stencil_info["definition"]

    dtypes = gt4py_config.dtypes.dict()
    dtypes[float] = gt4py_config.dtypes.float
    dtypes[int] = gt4py_config.dtypes.int
    externals = externals or {}

    kwargs = gt4py_config.backend_opts.copy()
    if gt4py_config.backend not in ("debug", "numpy", "gtc:numpy"):
        kwargs["verbose"] = gt4py_config.verbose

    return gtscript.stencil(
        gt4py_config.backend,
        definition,
        name=name,
        build_info=gt4py_config.build_info,
        dtypes=dtypes,
        externals=externals,
        rebuild=gt4py_config.rebuild,
        **kwargs,
    )
