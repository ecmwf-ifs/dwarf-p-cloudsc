# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

try:
    import cupy as cp
except ImportError:
    cp = np

if TYPE_CHECKING:
    from cloudsc4py.utils.typingx import Storage


def to_numpy(storage: Storage) -> np.ndarray:
    try:
        return storage.get()
    except AttributeError:
        return storage


def assign(lhs: Storage, rhs: Storage) -> None:
    if isinstance(lhs, cp.ndarray) and isinstance(rhs, np.ndarray):
        lhs[...] = cp.asarray(rhs)
    elif isinstance(lhs, np.ndarray) and isinstance(rhs, cp.ndarray):
        lhs[...] = rhs.get()
    else:
        lhs[...] = rhs
