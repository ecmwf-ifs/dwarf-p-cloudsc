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
from cloudsc4py.physics._stencils.fcttre import f_foedem, f_foeewm, f_foeldcpm
from cloudsc4py.utils.f2py import ported_function


@function_collection("f_cuadjtq_5")
@gtscript.function
def f_cuadjtq_5(qp, qsmix, t):
    from __externals__ import RETV

    qsat = min(f_foeewm(t) * qp, 0.5)
    cor = 1 / (1 - RETV * qsat)
    qsat *= cor
    cond = (qsmix - qsat) / (1 + qsat * cor * f_foedem(t))
    t += f_foeldcpm(t) * cond
    qsmix -= cond
    return qsmix, t


@ported_function(from_file="cloudsc_fortran/cloudsc2.F90", from_line=1297, to_line=1314)
@function_collection("f_cuadjtq")
@gtscript.function
def f_cuadjtq(ap, qsmix, t):
    qp = 1 / ap
    qsmix, t = f_cuadjtq_5(qp, qsmix, t)
    qsmix, t = f_cuadjtq_5(qp, qsmix, t)
    return qsmix, t
