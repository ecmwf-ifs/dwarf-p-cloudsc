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

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Optional, Union


PORTED_OBJECTS = {}


def ported_object(
    handle: Optional[Callable] = None,
    from_file: Optional[Union[str, Sequence[str]]] = None,
    from_line: Optional[int] = None,
    to_line: Optional[int] = None,
) -> Callable:
    if from_line is not None and to_line is not None:
        assert from_line <= to_line

    def core(obj):
        PORTED_OBJECTS[obj.__name__] = obj
        setattr(obj, "from_file", from_file)
        setattr(obj, "from_line", from_line)
        setattr(obj, "to_line", to_line)
        return obj

    if handle is not None:
        return core(handle)
    else:
        return core


# convenient aliases to improve readability
ported_class = ported_object
ported_function = ported_object
ported_method = ported_object
