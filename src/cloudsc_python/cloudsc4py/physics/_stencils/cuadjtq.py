# -*- coding: utf-8 -*-
from gt4py import gtscript

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
