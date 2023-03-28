# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""
Utility routines to dynamically load generated Python modules
"""

import sys
from pathlib import Path
from importlib import import_module, invalidate_caches, reload


__all__ = ['load_module']


def load_module(module, modpath=None):
    """
    Utility routine to dynamically load the requested Python module.
    """

    modpath = Path.cwd() if modpath is None else modpath
    modpath = str(Path(modpath).absolute())
    if modpath not in sys.path:
        sys.path.insert(0, modpath)
    if module in sys.modules:
        reload(sys.modules[module])
        return sys.modules[module]

    # Trigger the actual module import
    try:
        return import_module(module)
    except ModuleNotFoundError:
        # If module caching interferes, try again with clean caches
        invalidate_caches()
        return import_module(module)
