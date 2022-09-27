# -*- coding: utf-8 -*-
from __future__ import annotations

from sympl._core.time import Timer


class timing:
    def __init__(self, label: str) -> None:
        self.label = label

    def __enter__(self) -> type[Timer]:
        Timer.start(self.label)
        return Timer

    def __exit__(self, exc_type, exc_value, exc_tb) -> None:
        Timer.stop()
