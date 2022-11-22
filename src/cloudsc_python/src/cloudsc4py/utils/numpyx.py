# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cloudsc4py.utils.typingx import Storage


def to_numpy(storage: Storage) -> np.ndarray:
    try:
        return storage.asnumpy()
    except AttributeError:
        return storage
