# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from gt4py.cartesian import gtscript

from cloudsc4py.framework.stencil import function_collection
from cloudsc4py.physics._stencils.fcttre import f_foeeice, f_foeeliq
from cloudsc4py.utils.f2py import ported_function


@ported_function(from_file="common/include/fccld.func.h", from_line=26, to_line=27)
@function_collection("f_fokoop")
@gtscript.function
def f_fokoop(t):
    from __externals__ import RKOOP1, RKOOP2

    return min(RKOOP1 - RKOOP2 * t, f_foeeliq(t) / f_foeeice(t))
