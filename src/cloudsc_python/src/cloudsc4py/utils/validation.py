# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from cloudsc4py.utils.numpyx import to_numpy

if TYPE_CHECKING:
    from typing import Tuple

    from sympl._core.data_array import DataArray
    from sympl._core.typingx import DataArrayDict

    from cloudsc4py.utils.typingx import Storage


def validate_storage_2d(src: Storage, trg: Storage) -> bool:
    src_np = to_numpy(src)
    trg_np = to_numpy(trg)
    mi = min(src_np.shape[0], trg_np.shape[0])
    mj = min(src_np.shape[1], trg_np.shape[1])
    return np.allclose(src_np[:mi, :mj], trg_np[:mi, :mj], atol=1e-18, rtol=1e-12)


def validate_storage_3d(src: Storage, trg: Storage) -> bool:
    src_np = to_numpy(src)
    trg_np = to_numpy(trg)
    mi = min(src_np.shape[0], trg_np.shape[0])
    mj = min(src_np.shape[1], trg_np.shape[1])
    mk = min(src_np.shape[2], trg_np.shape[2])
    return np.allclose(src_np[:mi, :mj, :mk], trg_np[:mi, :mj, :mk], atol=1e-18, rtol=1e-12)


def validate_field(src: DataArray, trg: DataArray) -> bool:
    if src.ndim == 2:
        return validate_storage_2d(src.data, trg.data)
    elif src.ndim == 3:
        return validate_storage_3d(src.data, trg.data)
    else:
        raise ValueError("The field to validate must be either 2-d or 3-d.")


def validate(src: DataArrayDict, trg: DataArrayDict) -> Tuple[str]:
    return tuple(
        name
        for name in src
        if name in trg and name != "time" and not validate_field(src[name], trg[name])
    )
