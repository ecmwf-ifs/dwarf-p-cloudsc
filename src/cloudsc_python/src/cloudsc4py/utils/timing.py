# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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
