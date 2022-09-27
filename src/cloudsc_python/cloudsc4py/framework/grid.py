# -*- coding: utf-8 -*-
from __future__ import annotations
from functools import cached_property
import numpy as np


class DimSymbol:
    """Symbol identifying a dimension, e.g. I or I-1/2."""

    _instances: dict[int, DimSymbol] = {}

    name: str
    offset: float

    def __new__(cls, *args) -> DimSymbol:
        key = hash(args)
        if key not in cls._instances:
            cls._instances[key] = super().__new__(cls)
        return cls._instances[key]

    def __init__(self, name: str, offset: float) -> None:
        self.name = name
        self.offset = offset

    def __add__(self, other: float) -> DimSymbol:
        return DimSymbol(self.name, self.offset + other)

    def __sub__(self, other: float) -> DimSymbol:
        return self + (-other)

    def __repr__(self) -> str:
        if self.offset > 0:
            return f"{self.name} + {self.offset}"
        elif self.offset < 0:
            return f"{self.name} - {-self.offset}"
        else:
            return f"{self.name}"


I = DimSymbol("I", 0)
J = DimSymbol("J", 0)
K = DimSymbol("K", 0)


class Grid:
    """Grid of points."""

    def __init__(
        self, shape: tuple[int, ...], dims: tuple[str, ...], storage_shape: tuple[int, ...] = None
    ) -> None:
        assert len(shape) == len(dims)
        self.shape = shape
        self.dims = dims
        self.storage_shape = storage_shape or self.shape

    @cached_property
    def coords(self) -> tuple[np.ndarray, ...]:
        return tuple(np.arange(size) for size in self.storage_shape)


class ComputationalGrid:
    """A three-dimensional computational grid consisting of mass and staggered grid points."""

    grids: dict[tuple[DimSymbol, ...], Grid]

    def __init__(self, nx: int, ny: int, nz: int) -> None:
        self.grids = {
            (I, J, K): Grid((nx, ny, nz), ("x", "y", "z"), (nx, ny, nz + 1)),
            (I, J, K - 1 / 2): Grid((nx, ny, nz + 1), ("x", "y", "z_h")),
            (I, J): Grid((nx, ny), ("x", "y")),
            (K,): Grid((nz,), ("z",), (nz + 1,)),
        }


def get_mask(grid_id: tuple[DimSymbol, ...], data_shape: tuple[int, ...]) -> list[bool]:
    """Compute the mask for a storage defined over `grid_id`."""
    out = []
    for target_name in ("I", "J", "K"):
        for dim_symbol in grid_id:
            if dim_symbol.name == target_name:
                out.append(True)
                break
        else:
            out.append(False)
    return out + [True] * len(data_shape)
