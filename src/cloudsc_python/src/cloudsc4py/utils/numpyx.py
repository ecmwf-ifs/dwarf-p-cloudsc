# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from typing import TYPE_CHECKING

import gt4py

if TYPE_CHECKING:
    from cloudsc4py.utils.typingx import Storage


@dataclass(frozen=False)
class Status:
    enabled: bool = False


NUMPY_PATCH = Status()


class prepare_numpy:
    def __init__(self):
        self._patch_already_enabled = False

    def __enter__(self):
        if NUMPY_PATCH.enabled:
            self._patch_already_enabled = True
        else:
            NUMPY_PATCH.enabled = True
            self._patch_already_enabled = False
            gt4py.storage.prepare_numpy()
        return None

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self._patch_already_enabled:
            NUMPY_PATCH.enabled = False
            gt4py.storage.restore_numpy()


class restore_numpy:
    def __init__(self):
        self._patch_already_disabled = False

    def __enter__(self):
        if NUMPY_PATCH.enabled:
            NUMPY_PATCH.enabled = False
            self._patch_already_disabled = False
            gt4py.storage.restore_numpy()
        else:
            self._patch_already_disabled = True
        return None

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self._patch_already_disabled:
            NUMPY_PATCH.enabled = True
            gt4py.storage.prepare_numpy()


def to_numpy(storage: Storage) -> np.ndarray:
    try:
        storage.synchronize()
    except AttributeError:
        pass

    with restore_numpy():
        return np.asarray(storage)
