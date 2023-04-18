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


@function_collection("f_helper_0")
@gtscript.function
def f_helper_0(
    order,
    index1_ql,
    index1_qi,
    index1_qr,
    index1_qs,
    index1_qv,
    ratio_ql,
    ratio_qi,
    ratio_qr,
    ratio_qs,
    ratio_qv,
):
    minimum = 1e32

    if index1_ql and ratio_ql < minimum:
        order = 1
        minimum = ratio_ql
    if index1_qi and ratio_qi < minimum:
        order = 2
        minimum = ratio_qi
    if index1_qr and ratio_qr < minimum:
        order = 3
        minimum = ratio_qr
    if index1_qs and ratio_qs < minimum:
        order = 4
        minimum = ratio_qs
    if index1_qv and ratio_qv < minimum:
        order = 0

    if order == 1:
        index1_ql = False
    if order == 2:
        index1_qi = False
    if order == 3:
        index1_qr = False
    if order == 4:
        index1_qs = False
    if order == 0:
        index1_qv = False

    return order, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv


@function_collection("f_helper_1")
@gtscript.function
def f_helper_1(
    order,
    index3_ql_ql,
    index3_ql_qi,
    index3_ql_qr,
    index3_ql_qs,
    index3_ql_qv,
    index3_qi_ql,
    index3_qi_qi,
    index3_qi_qr,
    index3_qi_qs,
    index3_qi_qv,
    index3_qr_ql,
    index3_qr_qi,
    index3_qr_qr,
    index3_qr_qs,
    index3_qr_qv,
    index3_qs_ql,
    index3_qs_qi,
    index3_qs_qr,
    index3_qs_qs,
    index3_qs_qv,
    index3_qv_ql,
    index3_qv_qi,
    index3_qv_qr,
    index3_qv_qs,
    index3_qv_qv,
    ql,
    qi,
    qr,
    qs,
    qv,
    ratio_ql,
    ratio_qi,
    ratio_qr,
    ratio_qs,
    ratio_qv,
    sinksum_ql,
    sinksum_qi,
    sinksum_qr,
    sinksum_qs,
    sinksum_qv,
    solqa_ql_ql,
    solqa_ql_qi,
    solqa_ql_qr,
    solqa_ql_qs,
    solqa_ql_qv,
    solqa_qi_ql,
    solqa_qi_qi,
    solqa_qi_qr,
    solqa_qi_qs,
    solqa_qi_qv,
    solqa_qr_ql,
    solqa_qr_qi,
    solqa_qr_qr,
    solqa_qr_qs,
    solqa_qr_qv,
    solqa_qs_ql,
    solqa_qs_qi,
    solqa_qs_qr,
    solqa_qs_qs,
    solqa_qs_qv,
    solqa_qv_ql,
    solqa_qv_qi,
    solqa_qv_qr,
    solqa_qv_qs,
    solqa_qv_qv,
):
    from __externals__ import EPSEC

    # recalculate sum and scaling factor
    if order == 1:
        index3_ql_ql = solqa_ql_ql < 0.0
        index3_ql_qi = solqa_ql_qi < 0.0
        index3_ql_qr = solqa_ql_qr < 0.0
        index3_ql_qs = solqa_ql_qs < 0.0
        index3_ql_qv = solqa_ql_qv < 0.0
        sinksum_ql -= solqa_ql_ql + solqa_ql_qi + solqa_ql_qr + solqa_ql_qs + solqa_ql_qv
        mm = max(ql, EPSEC)
        rr = max(sinksum_ql, mm)
        ratio_ql = mm / rr
    elif order == 2:
        index3_qi_ql = solqa_qi_ql < 0.0
        index3_qi_qi = solqa_qi_qi < 0.0
        index3_qi_qr = solqa_qi_qr < 0.0
        index3_qi_qs = solqa_qi_qs < 0.0
        index3_qi_qv = solqa_qi_qv < 0.0
        sinksum_qi -= solqa_qi_ql + solqa_qi_qi + solqa_qi_qr + solqa_qi_qs + solqa_qi_qv
        mm = max(qi, EPSEC)
        rr = max(sinksum_qi, mm)
        ratio_qi = mm / rr
    elif order == 3:
        index3_qr_ql = solqa_qr_ql < 0.0
        index3_qr_qi = solqa_qr_qi < 0.0
        index3_qr_qr = solqa_qr_qr < 0.0
        index3_qr_qs = solqa_qr_qs < 0.0
        index3_qr_qv = solqa_qr_qv < 0.0
        sinksum_qr -= solqa_qr_ql + solqa_qr_qi + solqa_qr_qr + solqa_qr_qs + solqa_qr_qv
        mm = max(qr, EPSEC)
        rr = max(sinksum_qr, mm)
        ratio_qr = mm / rr
    elif order == 4:
        index3_qs_ql = solqa_qs_ql < 0.0
        index3_qs_qi = solqa_qs_qi < 0.0
        index3_qs_qr = solqa_qs_qr < 0.0
        index3_qs_qs = solqa_qs_qs < 0.0
        index3_qs_qv = solqa_qs_qv < 0.0
        sinksum_qs -= solqa_qs_ql + solqa_qs_qi + solqa_qs_qr + solqa_qs_qs + solqa_qs_qv
        mm = max(qs, EPSEC)
        rr = max(sinksum_qs, mm)
        ratio_qs = mm / rr
    elif order == 0:
        index3_qv_ql = solqa_qv_ql < 0.0
        index3_qv_qi = solqa_qv_qi < 0.0
        index3_qv_qr = solqa_qv_qr < 0.0
        index3_qv_qs = solqa_qv_qs < 0.0
        index3_qv_qv = solqa_qv_qv < 0.0
        sinksum_qv -= solqa_qv_ql + solqa_qv_qi + solqa_qv_qr + solqa_qv_qs + solqa_qv_qv
        mm = max(qv, EPSEC)
        rr = max(sinksum_qv, mm)
        ratio_qv = mm / rr

    # scale
    if order == 1:
        if index3_ql_ql:
            solqa_ql_ql *= ratio_ql
            solqa_ql_ql *= ratio_ql
        if index3_ql_qi:
            solqa_ql_qi *= ratio_ql
            solqa_qi_ql *= ratio_ql
        if index3_ql_qr:
            solqa_ql_qr *= ratio_ql
            solqa_qr_ql *= ratio_ql
        if index3_ql_qs:
            solqa_ql_qs *= ratio_ql
            solqa_qs_ql *= ratio_ql
        if index3_ql_qv:
            solqa_ql_qv *= ratio_ql
            solqa_qv_ql *= ratio_ql
    elif order == 2:
        if index3_qi_ql:
            solqa_qi_ql *= ratio_qi
            solqa_ql_qi *= ratio_qi
        if index3_qi_qi:
            solqa_qi_qi *= ratio_qi
            solqa_qi_qi *= ratio_qi
        if index3_qi_qr:
            solqa_qi_qr *= ratio_qi
            solqa_qr_qi *= ratio_qi
        if index3_qi_qs:
            solqa_qi_qs *= ratio_qi
            solqa_qs_qi *= ratio_qi
        if index3_qi_qv:
            solqa_qi_qv *= ratio_qi
            solqa_qv_qi *= ratio_qi
    elif order == 3:
        if index3_qr_ql:
            solqa_qr_ql *= ratio_qr
            solqa_ql_qr *= ratio_qr
        if index3_qr_qi:
            solqa_qr_qi *= ratio_qr
            solqa_qi_qr *= ratio_qr
        if index3_qr_qr:
            solqa_qr_qr *= ratio_qr
            solqa_qr_qr *= ratio_qr
        if index3_qr_qs:
            solqa_qr_qs *= ratio_qr
            solqa_qs_qr *= ratio_qr
        if index3_qr_qv:
            solqa_qr_qv *= ratio_qr
            solqa_qv_qr *= ratio_qr
    elif order == 4:
        if index3_qs_ql:
            solqa_qs_ql *= ratio_qs
            solqa_ql_qs *= ratio_qs
        if index3_qs_qi:
            solqa_qs_qi *= ratio_qs
            solqa_qi_qs *= ratio_qs
        if index3_qs_qr:
            solqa_qs_qr *= ratio_qs
            solqa_qr_qs *= ratio_qs
        if index3_qs_qs:
            solqa_qs_qs *= ratio_qs
            solqa_qs_qs *= ratio_qs
        if index3_qs_qv:
            solqa_qs_qv *= ratio_qs
            solqa_qv_qs *= ratio_qs
    elif order == 0:
        if index3_qv_ql:
            solqa_qv_ql *= ratio_qv
            solqa_ql_qv *= ratio_qv
        if index3_qv_qi:
            solqa_qv_qi *= ratio_qv
            solqa_qi_qv *= ratio_qv
        if index3_qv_qr:
            solqa_qv_qr *= ratio_qv
            solqa_qr_qv *= ratio_qv
        if index3_qv_qs:
            solqa_qv_qs *= ratio_qv
            solqa_qs_qv *= ratio_qv
        if index3_qv_qv:
            solqa_qv_qv *= ratio_qv
            solqa_qv_qv *= ratio_qv

    return ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv
