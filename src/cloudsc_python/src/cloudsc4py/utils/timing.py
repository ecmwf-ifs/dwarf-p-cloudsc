# -*- coding: utf-8 -*-
from __future__ import annotations
from typing import TYPE_CHECKING

from sympl._core.time import Timer

if TYPE_CHECKING:
    from typing import Type


class timing:
    def __init__(self, label: str) -> None:
        self.label = label

    def __enter__(self) -> Type[Timer]:
        Timer.start(self.label)
        return Timer

    def __exit__(self, exc_type, exc_value, exc_tb) -> None:
        Timer.stop()
