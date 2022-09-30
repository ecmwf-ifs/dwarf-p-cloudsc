# -*- coding: utf-8 -*-
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
    from typing import Any

    from gt4py import StencilObject
    from gt4py.storage import Storage
    from sympl._core.typingx import PropertyDict

    from cloudsc4py.framework.grid import ComputationalGrid


class ComputationalGridComponent:
    def __init__(self, computational_grid: ComputationalGrid, *, gt4py_config: GT4PyConfig) -> None:
        self.computational_grid = computational_grid
        self.gt4py_config = gt4py_config

    def compile_stencil(
        self, name: str, externals: Optional[dict[str, Any]] = None
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
        ...

    def allocate_diagnostic(self, name: str) -> Storage:
        return self.allocate(name, self.diagnostic_properties)

    @cached_property
    def diagnostic_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._diagnostic_properties)

    @abstractmethod
    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        ...


class ImplicitTendencyComponent(ComputationalGridComponent, SymplImplicitTendencyComponent):
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
        ...

    def allocate_tendency(self, name: str) -> Storage:
        return self.allocate(name, self.tendency_properties)

    @cached_property
    def tendency_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._tendency_properties)

    @abstractmethod
    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        ...

    def allocate_diagnostic(self, name: str) -> Storage:
        return self.allocate(name, self.diagnostic_properties)

    @cached_property
    def diagnostic_properties(self) -> PropertyDict:
        return self.fill_properties_with_dims(self._diagnostic_properties)

    @abstractmethod
    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        ...
