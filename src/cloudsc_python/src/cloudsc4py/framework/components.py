# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
from abc import abstractmethod
from functools import cached_property
from typing import Optional, TYPE_CHECKING

from sympl._core.core_components import (
    DiagnosticComponent as SymplDiagnosticComponent,
    ImplicitTendencyComponent as SymplImplicitTendencyComponent,
)

from cloudsc4py.framework.config import GT4PyConfig
from cloudsc4py.framework.stencil import compile_stencil
from cloudsc4py.framework.storage import get_data_shape_from_name, get_dtype_from_name, zeros

if TYPE_CHECKING:
    from typing import Any, Dict

    from gt4py import StencilObject
    from gt4py.storage import Storage
    from sympl._core.typingx import PropertyDict

    from cloudsc4py.framework.grid import ComputationalGrid


class ComputationalGridComponent:
    """Model component defined over a computational grid."""

    def __init__(self, computational_grid: ComputationalGrid, *, gt4py_config: GT4PyConfig) -> None:
        self.computational_grid = computational_grid
        self.gt4py_config = gt4py_config

    def compile_stencil(
        self, name: str, externals: Optional[Dict[str, Any]] = None
    ) -> StencilObject:
        return compile_stencil(name, self.gt4py_config, externals)

    def fill_properties_with_dims(self, properties: PropertyDict) -> PropertyDict:
        for field_name, field_prop in properties.items():
            field_prop["dims"] = self.computational_grid.grids[field_prop["grid"]].dims
        return properties

    def allocate(self, name: str, properties: PropertyDict) -> Storage:
        data_shape = get_data_shape_from_name(name)
        dtype = get_dtype_from_name(name)
        return zeros(
            self.computational_grid,
            properties[name]["grid"],
            data_shape,
            gt4py_config=self.gt4py_config,
            dtype=dtype,
        )


class DiagnosticComponent(ComputationalGridComponent, SymplDiagnosticComponent):
    """Grid-aware variant of Sympl's ``DiagnosticComponent``."""

    def __init__(
        self,
        computational_grid: ComputationalGrid,
        *,
        enable_checks: bool = True,
        gt4py_config: GT4PyConfig,
    ) -> None:
        super().__init__(computational_grid, gt4py_config=gt4py_config)
        super(ComputationalGridComponent, self).__init__(enable_checks=enable_checks)

    @cached_property
    def input_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._input_properties)

    @abstractmethod
    @cached_property
    def _input_properties(self) -> PropertyDict:
        """
        Dictionary where each key is the name of an input field, and the corresponding value is a
        dictionary specifying the units for that field ('units') and the identifier of the grid over
        which it is defined ('grid').
        """
        ...

    def allocate_diagnostic(self, name: str) -> Storage:
        return self.allocate(name, self.diagnostic_properties)

    @cached_property
    def diagnostic_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._diagnostic_properties)

    @abstractmethod
    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        """
        Dictionary where each key is the name of a field diagnosed by the component, and the
        corresponding value is a dictionary specifying the units for that field ('units') and the
        identifier of the grid over which it is defined ('grid').
        """
        ...


class ImplicitTendencyComponent(ComputationalGridComponent, SymplImplicitTendencyComponent):
    """Grid-aware variant of Sympl's ``ImplicitTendencyComponent``."""

    def __init__(
        self,
        computational_grid: ComputationalGrid,
        *,
        enable_checks: bool = True,
        gt4py_config: GT4PyConfig,
    ) -> None:
        super().__init__(computational_grid, gt4py_config=gt4py_config)
        super(ComputationalGridComponent, self).__init__(enable_checks=enable_checks)

    @cached_property
    def input_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._input_properties)

    @abstractmethod
    @cached_property
    def _input_properties(self) -> PropertyDict:
        """
        Dictionary where each key is the name of an input field, and the corresponding value is a
        dictionary specifying the units for that field ('units') and the identifier of the grid over
        which it is defined ('grid').
        """
        ...

    def allocate_tendency(self, name: str) -> Storage:
        return self.allocate(name, self.tendency_properties)

    @cached_property
    def tendency_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._tendency_properties)

    @abstractmethod
    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        """
        Dictionary where each key is the name of a tendency field computed by the component, and the
        corresponding value is a dictionary specifying the units for that field ('units') and the
        identifier of the grid over which it is defined ('grid').
        """
        ...

    def allocate_diagnostic(self, name: str) -> Storage:
        return self.allocate(name, self.diagnostic_properties)

    @cached_property
    def diagnostic_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._diagnostic_properties)

    @abstractmethod
    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        """
        Dictionary where each key is the name of a field diagnosed by the component, and the
        corresponding value is a dictionary specifying the units for that field ('units') and the
        identifier of the grid over which it is defined ('grid').
        """
        ...
