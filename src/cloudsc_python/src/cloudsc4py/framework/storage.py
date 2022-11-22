# -*- coding: utf-8 -*-
from __future__ import annotations
from contextlib import contextmanager
import numpy as np
from typing import TYPE_CHECKING

import gt4py
from sympl._core.data_array import DataArray

from cloudsc4py.framework.grid import get_mask

if TYPE_CHECKING:
    from typing import Literal, Optional

    from cloudsc4py.framework.config import GT4PyConfig
    from cloudsc4py.framework.grid import ComputationalGrid, DimSymbol


def zeros(
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    data_shape: Optional[tuple[int, ...]] = None,
    *,
    gt4py_config: GT4PyConfig,
    dtype: Literal["bool", "float", "int"],
) -> gt4py.storage.Storage:
    grid = computational_grid.grids[grid_id]
    data_shape = data_shape or ()
    shape = grid.storage_shape + data_shape
    dtype = gt4py_config.dtypes.dict()[dtype]
    return gt4py.storage.zeros(shape, dtype, backend=gt4py_config.backend)


def get_data_array(
    buffer: gt4py.storage.Storage,
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    units: str,
    data_dims: Optional[tuple[str, ...]] = None,
) -> DataArray:
    grid = computational_grid.grids[grid_id]
    data_dims = data_dims or ()
    dims = grid.dims + data_dims
    coords = grid.coords + tuple(
        np.arange(data_size) for data_size in buffer.shape[len(grid.dims) :]
    )
    return DataArray(buffer, dims=dims, coords=coords, attrs={"units": units})


def allocate_data_array(
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    units: str,
    data_shape: Optional[tuple[int, ...]] = None,
    data_dims: Optional[tuple[str, ...]] = None,
    *,
    gt4py_config: GT4PyConfig,
    dtype: Literal["bool", "float", "int"],
) -> DataArray:
    buffer = zeros(
        computational_grid, grid_id, data_shape=data_shape, gt4py_config=gt4py_config, dtype=dtype
    )
    return get_data_array(buffer, computational_grid, grid_id, units, data_dims=data_dims)


def get_dtype_from_name(field_name: str) -> str:
    if field_name.startswith("b"):
        return "bool"
    elif field_name.startswith("f"):
        return "float"
    elif field_name.startswith("i"):
        return "int"
    else:
        raise RuntimeError(f"Cannot retrieve dtype for field `{field_name}`.")


def get_data_shape_from_name(field_name: str) -> tuple[int]:
    data_dims = field_name.split("_", maxsplit=1)[0][1:]
    out = tuple(int(c) for c in data_dims)
    return out


TEMPORARY_STORAGE_POOL: dict[int, list[gt4py.storage.Storage]] = {}


@contextmanager
def managed_temporary_storage(
    computational_grid: ComputationalGrid,
    *args: tuple[tuple[DimSymbol, ...], Literal["bool", "float", "int"]],
    gt4py_config: GT4PyConfig,
):
    grid_hashes = []
    storages = []
    for grid_id, dtype in args:
        grid = computational_grid.grids[grid_id]
        grid_hash = hash((grid.shape + grid_id, dtype))
        pool = TEMPORARY_STORAGE_POOL.setdefault(grid_hash, [])
        if len(pool) > 0:
            storage = pool.pop()
        else:
            storage = zeros(computational_grid, grid_id, gt4py_config=gt4py_config, dtype=dtype)
        grid_hashes.append(grid_hash)
        storages.append(storage)

    try:
        if len(storages) == 1:
            yield storages[0]
        else:
            yield storages
    finally:
        for grid_hash, storage in zip(grid_hashes, storages):
            TEMPORARY_STORAGE_POOL[grid_hash].append(storage)


@contextmanager
def managed_temporary_storage_pool():
    """
    Context manager clearing the pool of temporary storages on entry and exit.

    Useful when running multiple simulations using different backends within the same session.
    All simulations using the same backend should be wrapped by this context manager.
    """
    try:
        TEMPORARY_STORAGE_POOL.clear()
        yield None
    finally:
        for grid_hash, storages in TEMPORARY_STORAGE_POOL.items():
            num_storages = len(storages)
            for _ in range(num_storages):
                storage = storages.pop()
                del storage
        TEMPORARY_STORAGE_POOL.clear()
