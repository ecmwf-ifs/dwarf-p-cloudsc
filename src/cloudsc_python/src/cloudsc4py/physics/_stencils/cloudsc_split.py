# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations

from gt4py.cartesian.gtscript import Field, IJ, K

from cloudsc4py.framework.stencil import stencil_collection
from cloudsc4py.physics._stencils.cuadjtq import f_cuadjtq
from cloudsc4py.physics._stencils.fccld import f_fokoop
from cloudsc4py.physics._stencils.fcttre import (
    f_foealfa,
    f_foedelta,
    f_foedem,
    f_foeeice,
    f_foeeliq,
    f_foeewm,
    f_foeldcpm,
)
from cloudsc4py.physics._stencils.helpers import f_helper_0, f_helper_1


@stencil_collection("cloudsc_tendencies")
def cloudsc_tendencies(
    in_a: Field["float"],
    in_ap: Field["float"],
    in_aph: Field["float"],  # staggered
    in_ccn: Field["float"],
    in_convection_on: Field[IJ, "bool"],
    in_convection_type: Field[IJ, "int"],
    in_hrlw: Field["float"],
    in_hrsw: Field["float"],
    in_icrit_aer: Field["float"],
    in_lcrit_aer: Field["float"],
    in_lsm: Field[IJ, "float"],
    in_lu: Field["float"],
    in_lude: Field["float"],
    in_mfd: Field["float"],
    in_mfu: Field["float"],
    in_nice: Field["float"],
    in_qi: Field["float"],
    in_ql: Field["float"],
    in_qr: Field["float"],
    in_qs: Field["float"],
    in_qv: Field["float"],
    in_re_ice: Field["float"],
    in_snde: Field["float"],
    in_supsat: Field["float"],
    in_t: Field["float"],
    in_tnd_tmp_a: Field["float"],
    in_tnd_tmp_qi: Field["float"],
    in_tnd_tmp_ql: Field["float"],
    in_tnd_tmp_qr: Field["float"],
    in_tnd_tmp_qs: Field["float"],
    in_tnd_tmp_qv: Field["float"],
    in_tnd_tmp_t: Field["float"],
    in_w: Field["float"],
    out_covptot: Field["float"],
    out_foealfa: Field["float"],
    out_lneg_qi: Field["float"],
    out_lneg_ql: Field["float"],
    out_lneg_qr: Field["float"],
    out_lneg_qs: Field["float"],
    out_lude: Field["float"],
    out_pfplsi: Field["float"],
    out_pfplsl: Field["float"],
    out_pfplsr: Field["float"],
    out_pfplss: Field["float"],
    out_qi0: Field["float"],
    out_qin: Field["float"],
    out_ql0: Field["float"],
    out_qln: Field["float"],
    out_qr0: Field["float"],
    out_qrn: Field["float"],
    out_qs0: Field["float"],
    out_qsn: Field["float"],
    out_rainfrac_toprfz: Field[IJ, "float"],
    out_tnd_loc_a: Field["float"],
    out_tnd_loc_qi: Field["float"],
    out_tnd_loc_ql: Field["float"],
    out_tnd_loc_qr: Field["float"],
    out_tnd_loc_qs: Field["float"],
    out_tnd_loc_qv: Field["float"],
    out_tnd_loc_t: Field["float"],
    tmp_aph_s: Field[IJ, "float"],
    tmp_cldtopdist: Field[IJ, "float"],
    tmp_covpmax: Field[IJ, "float"],
    tmp_covptot: Field[IJ, "float"],
    tmp_klevel: Field[K, "int"],
    tmp_paphd: Field[IJ, "float"],
    tmp_rainliq: Field[IJ, "bool"],
    tmp_trpaus: Field[IJ, "float"],
    *,
    dt: "float",
):
    from __externals__ import (
        DEPICE,
        EPSEC,
        EPSILON,
        EVAPRAIN,
        EVAPSNOW,
        FALLQI,
        FALLQL,
        FALLQR,
        FALLQS,
        FALLQV,
        LAERICEAUTO,
        LAERICESED,
        LAERLIQAUTOLSP,
        LAERLIQCOLL,
        NCLDTOP,
        NLEV,
        NSSOPT,
        PHASEQI,
        PHASEQL,
        PHASEQR,
        PHASEQS,
        PHASEQV,
        R4IES,
        R4LES,
        R5IES,
        R5LES,
        RALFDCP,
        RALSDCP,
        RALVDCP,
        RAMID,
        RAMIN,
        RCCN,
        RCL_APB1,
        RCL_APB2,
        RCL_APB3,
        RCL_CDENOM1,
        RCL_CDENOM2,
        RCL_CDENOM3,
        RCL_CONST1I,
        RCL_CONST1R,
        RCL_CONST1S,
        RCL_CONST2I,
        RCL_CONST2R,
        RCL_CONST2S,
        RCL_CONST3I,
        RCL_CONST3R,
        RCL_CONST3S,
        RCL_CONST4I,
        RCL_CONST4R,
        RCL_CONST4S,
        RCL_CONST5I,
        RCL_CONST5R,
        RCL_CONST5S,
        RCL_CONST6I,
        RCL_CONST6R,
        RCL_CONST6S,
        RCL_CONST7S,
        RCL_CONST8S,
        RCL_FAC1,
        RCL_FAC2,
        RCL_FZRAB,
        RCL_KK_cloud_num_land,
        RCL_KK_cloud_num_sea,
        RCL_KKAac,
        RCL_KKAau,
        RCL_KKBac,
        RCL_KKBaun,
        RCL_KKBauq,
        RCLCRIT_LAND,
        RCLCRIT_SEA,
        RCLDIFF,
        RCLDIFF_CONVI,
        RCLDTOPCF,
        RCOVPMIN,
        RD,
        RDCP,
        RDENSREF,
        RDEPLIQREFDEPTH,
        RDEPLIQREFRATE,
        RETV,
        RG,
        RICEINIT,
        RKCONV,
        RKOOPTAU,
        RLCRITSNOW,
        RLDCP,
        RLMIN,
        RLSTT,
        RLVTT,
        RNICE,
        RPECONS,
        RPRC1,
        RPRECRHMAX,
        RSNOWLIN1,
        RSNOWLIN2,
        RTAUMEL,
        RTHOMO,
        RTT,
        RV,
        RVRFACTOR,
        TW1,
        TW2,
        TW3,
        TW4,
        TW5,
        VQI,
        VQL,
        VQR,
        VQS,
        VQV,
        WARMRAIN,
    )

    with computation(FORWARD), interval(0, 1):
        # zero arrays
        out_rainfrac_toprfz[0, 0] = 0.0
        tmp_cldtopdist[0, 0] = 0.0
        tmp_covpmax[0, 0] = 0.0
        tmp_covptot[0, 0] = 0.0
        tmp_paphd[0, 0] = 0.0
        tmp_rainliq[0, 0] = True
        tmp_trpaus[0, 0] = 0.0

    with computation(FORWARD), interval(...):
        # === 1: initial values for variables
        # --- initialization of output tendencies
        out_tnd_loc_t[0, 0, 0] = 0
        out_tnd_loc_a[0, 0, 0] = 0
        out_tnd_loc_ql[0, 0, 0] = 0
        out_tnd_loc_qr[0, 0, 0] = 0
        out_tnd_loc_qi[0, 0, 0] = 0
        out_tnd_loc_qs[0, 0, 0] = 0
        out_tnd_loc_qv[0, 0, 0] = 0

        # --- non CLV initialization
        t = in_t[0, 0, 0] + dt * in_tnd_tmp_t[0, 0, 0]
        a = in_a[0, 0, 0] + dt * in_tnd_tmp_a[0, 0, 0]
        a0 = a

        # --- initialization for CLV family
        ql = in_ql[0, 0, 0] + dt * in_tnd_tmp_ql[0, 0, 0]
        out_ql0[0, 0, 0] = ql
        qi = in_qi[0, 0, 0] + dt * in_tnd_tmp_qi[0, 0, 0]
        out_qi0[0, 0, 0] = qi
        qr = in_qr[0, 0, 0] + dt * in_tnd_tmp_qr[0, 0, 0]
        out_qr0[0, 0, 0] = qr
        qs = in_qs[0, 0, 0] + dt * in_tnd_tmp_qs[0, 0, 0]
        out_qs0[0, 0, 0] = qs
        qv = in_qv[0, 0, 0] + dt * in_tnd_tmp_qv[0, 0, 0]

        # --- zero arrays
        out_lneg_ql[0, 0, 0] = 0.0
        out_lneg_qi[0, 0, 0] = 0.0
        out_lneg_qr[0, 0, 0] = 0.0
        out_lneg_qs[0, 0, 0] = 0.0

        # --- tidy up very small cloud cover or total cloud water
        expr1 = ql + qi
        if expr1 < RLMIN or a < RAMIN:
            # evaporate small cloud liquid water amounts
            out_lneg_ql[0, 0, 0] += ql
            qadj = ql / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            out_tnd_loc_t[0, 0, 0] -= RALVDCP * qadj
            qv += ql
            ql = 0.0

            # evaporate small cloud ice water amounts
            out_lneg_qi[0, 0, 0] += qi
            qadj = qi / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            out_tnd_loc_t[0, 0, 0] -= RALSDCP * qadj
            qv += qi
            qi = 0.0

            # set cloud cover to zero
            a = 0.0

        # --- tidy up small CLV variables: ql
        if ql < RLMIN:
            out_lneg_ql[0, 0, 0] += ql
            qadj = ql / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            if __INLINED(PHASEQL == 1):
                out_tnd_loc_t[0, 0, 0] -= RALVDCP * qadj
            elif __INLINED(PHASEQL == 2):
                out_tnd_loc_t[0, 0, 0] -= RALSDCP * qadj
            qv += ql
            ql = 0.0

        # --- tidy up small CLV variables: qi
        if qi < RLMIN:
            out_lneg_qi[0, 0, 0] += qi
            qadj = qi / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            if __INLINED(PHASEQI == 1):
                out_tnd_loc_t[0, 0, 0] -= RALVDCP * qadj
            elif __INLINED(PHASEQI == 2):
                out_tnd_loc_t[0, 0, 0] -= RALSDCP * qadj
            qv += qi
            qi = 0.0

        # --- tidy up small CLV variables: qr
        if qr < RLMIN:
            out_lneg_qr[0, 0, 0] += qr
            qadj = qr / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            if __INLINED(PHASEQR == 1):
                out_tnd_loc_t[0, 0, 0] -= RALVDCP * qadj
            elif __INLINED(PHASEQR == 2):
                out_tnd_loc_t[0, 0, 0] -= RALSDCP * qadj
            qv += qr
            qr = 0.0

        # --- tidy up small CLV variables: qs
        if qs < RLMIN:
            out_lneg_qs[0, 0, 0] += qs
            qadj = qs / dt
            out_tnd_loc_qv[0, 0, 0] += qadj
            if __INLINED(PHASEQS == 1):
                out_tnd_loc_t[0, 0, 0] -= RALVDCP * qadj
            elif __INLINED(PHASEQS == 2):
                out_tnd_loc_t[0, 0, 0] -= RALSDCP * qadj
            qv += qs
            qs = 0.0

        # --- define saturation values
        # --- old *diagnostic* mixed phase saturation
        foealfa = f_foealfa(t)
        out_foealfa[0, 0, 0] = foealfa
        foeewmt = min(f_foeewm(t) / in_ap[0, 0, 0], 0.5)
        qsmix = foeewmt / (1 - RETV * foeewmt)

        # --- ice saturation T < 273K
        # --- liquid water saturation for T > 273K
        alfa = f_foedelta(t)
        foeew = min((alfa * f_foeeliq(t) + (1 - alfa) * f_foeeice(t)) / in_ap[0, 0, 0], 0.5)
        qsice = foeew / (1 - RETV * foeew)

        # --- liquid water saturation
        foeeliqt = min(f_foeeliq(t) / in_ap[0, 0, 0], 0.5)
        qsliq = foeeliqt / (1 - RETV * foeeliqt)

        # --- ensure cloud fraction is between 0 and 1
        a = max(0, min(1, a))

        # --- calculate liq/ice fractions (no longer a diagnostic relationship)
        li = ql + qi
        if li > RLMIN:
            liqfrac = ql / li
            icefrac = 1 - liqfrac
        else:
            liqfrac = 0.0
            icefrac = 0.0

    # === 2: constants and parameters
    # --- find tropopause level
    with computation(FORWARD), interval(0, 1):
        tmp_trpaus[0, 0] = 0.1
        tmp_paphd[0, 0] = 1 / tmp_aph_s[0, 0]
    with computation(FORWARD), interval(0, -1):
        sig = in_ap[0, 0, 0] * tmp_paphd[0, 0]
        if sig > 0.1 and sig < 0.4 and t[0, 0, 0] > t[0, 0, 1]:
            tmp_trpaus[0, 0] = sig

    # === 3: physics
    # --- main vertical loop
    with computation(FORWARD):
        with interval(0, NCLDTOP - 1):
            # --- initialize variables
            out_lude[0, 0, 0] = in_lude[0, 0, 0]
            out_pfplsl[0, 0, 0] = 0.0
            out_pfplsi[0, 0, 0] = 0.0
            out_pfplsr[0, 0, 0] = 0.0
            out_pfplss[0, 0, 0] = 0.0
            pfplsv = 0.0
            out_qln[0, 0, 0] = 0.0
            out_qin[0, 0, 0] = 0.0
            out_qrn[0, 0, 0] = 0.0
            out_qsn[0, 0, 0] = 0.0
            qvn = 0.0
            anew = 0.0
        with interval(NCLDTOP - 1, None):
            # *** 3.0: initialize variables
            # --- first guess microphysics
            qlfg = ql
            qifg = qi
            qrfg = qr
            qsfg = qs
            qvfg = qv

            convsink_ql = 0.0
            convsink_qi = 0.0
            convsink_qr = 0.0
            convsink_qs = 0.0
            convsrce_ql = 0.0
            convsrce_qi = 0.0
            convsrce_qr = 0.0
            convsrce_qs = 0.0
            convsrce_qv = 0.0
            fallsrce_ql = 0.0
            fallsrce_qi = 0.0
            fallsrce_qr = 0.0
            fallsrce_qs = 0.0
            index1_ql = True
            index1_qi = True
            index1_qr = True
            index1_qs = True
            index1_qv = True
            index3_ql_ql = False
            index3_ql_qi = False
            index3_ql_qr = False
            index3_ql_qs = False
            index3_ql_qv = False
            index3_qi_ql = False
            index3_qi_qi = False
            index3_qi_qr = False
            index3_qi_qs = False
            index3_qi_qv = False
            index3_qr_ql = False
            index3_qr_qi = False
            index3_qr_qr = False
            index3_qr_qs = False
            index3_qr_qv = False
            index3_qs_ql = False
            index3_qs_qi = False
            index3_qs_qr = False
            index3_qs_qs = False
            index3_qs_qv = False
            index3_qv_ql = False
            index3_qv_qi = False
            index3_qv_qr = False
            index3_qv_qs = False
            index3_qv_qv = False
            lcust_ql = 0.0
            lcust_qi = 0.0
            lcust_qr = 0.0
            lcust_qs = 0.0
            lcust_qv = 0.0
            ldefr = 0.0
            lfinalsum = 0.0
            order_ql = -999
            order_qi = -999
            order_qr = -999
            order_qs = -999
            order_qv = -999
            psupsatsrce_ql = 0.0
            psupsatsrce_qi = 0.0
            psupsatsrce_qr = 0.0
            psupsatsrce_qs = 0.0
            qpretot = 0.0
            solab = 0.0
            solac = 0.0
            solqa_ql_ql = 0.0
            solqa_ql_qi = 0.0
            solqa_ql_qr = 0.0
            solqa_ql_qs = 0.0
            solqa_ql_qv = 0.0
            solqa_qi_ql = 0.0
            solqa_qi_qi = 0.0
            solqa_qi_qr = 0.0
            solqa_qi_qs = 0.0
            solqa_qi_qv = 0.0
            solqa_qr_ql = 0.0
            solqa_qr_qi = 0.0
            solqa_qr_qr = 0.0
            solqa_qr_qs = 0.0
            solqa_qr_qv = 0.0
            solqa_qs_ql = 0.0
            solqa_qs_qi = 0.0
            solqa_qs_qr = 0.0
            solqa_qs_qs = 0.0
            solqa_qs_qv = 0.0
            solqa_qv_ql = 0.0
            solqa_qv_qi = 0.0
            solqa_qv_qr = 0.0
            solqa_qv_qs = 0.0
            solqa_qv_qv = 0.0
            solqb_ql_ql = 0.0
            solqb_ql_qi = 0.0
            solqb_ql_qr = 0.0
            solqb_ql_qs = 0.0
            solqb_ql_qv = 0.0
            solqb_qi_ql = 0.0
            solqb_qi_qi = 0.0
            solqb_qi_qr = 0.0
            solqb_qi_qs = 0.0
            solqb_qi_qv = 0.0
            solqb_qr_ql = 0.0
            solqb_qr_qi = 0.0
            solqb_qr_qr = 0.0
            solqb_qr_qs = 0.0
            solqb_qr_qv = 0.0
            solqb_qs_ql = 0.0
            solqb_qs_qi = 0.0
            solqb_qs_qr = 0.0
            solqb_qs_qs = 0.0
            solqb_qs_qv = 0.0
            solqb_qv_ql = 0.0
            solqb_qv_qi = 0.0
            solqb_qv_qr = 0.0
            solqb_qv_qs = 0.0
            solqb_qv_qv = 0.0

            # derived variables needed
            dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
            gdp = RG / dp
            rho = in_ap[0, 0, 0] / (RD * t)
            dtgdp = dt * gdp
            rdtgdp = dp / (RG * dt)

            # --- calculate dqs/dT correction factor
            # liquid
            facw = R5LES / (t - R4LES) ** 2
            cor = 1 / (1 - RETV * foeeliqt)
            dqsliqdt = facw * cor * qsliq
            corqsliq = 1 + RALVDCP * dqsliqdt

            # ice
            faci = R5IES / (t - R4IES) ** 2
            cor = 1 / (1 - RETV * foeew)
            dqsicedt = faci * cor * qsice
            corqsice = 1 + RALSDCP * dqsicedt

            # diagnostic mixed
            fac = out_foealfa[0, 0, 0] * facw + (1 - out_foealfa[0, 0, 0]) * faci
            cor = 1 / (1 - RETV * foeewmt)
            dqsmixdt = fac * cor * qsmix
            corqsmix = 1 + f_foeldcpm(t) * dqsmixdt

            # evaporation/sublimation limits
            evaplimmix = max((qsmix - qv) / corqsmix, 0.0)
            evaplimice = max((qsice - qv) / corqsice, 0.0)

            # --- in-cloud condensate amount
            tmpa = 1 / max(a, EPSEC)
            liqcld = ql * tmpa
            icecld = qi * tmpa
            licld = liqcld + icecld

            # --- evaporate very small amounts of liquid...
            if ql < RLMIN:
                solqa_qv_ql += ql
                solqa_ql_qv -= ql

            # --- ...and ice
            if qi < RLMIN:
                solqa_qv_qi += qi
                solqa_qi_qv -= qi

            # *** 3.1: ice supersaturation adjustment
            # --- supersaturation limit (from Koop)
            fokoop = f_fokoop(t)

            if t >= RTT or NSSOPT == 0:
                fac = 1.0
                faci = 1.0
            else:
                fac = a + fokoop * (1 - a)
                faci = dt / RKOOPTAU

            # calculate supersaturation to add to cloud
            if a > 1 - RAMIN:
                supsat = max((qv - fac * qsice) / corqsice, 0.0)
            else:
                # calculate environmental humidity supersaturation
                qp1env = (qv - a * qsice) / max(1 - a, EPSILON)
                supsat = max((1 - a) * (qp1env - fac * qsice) / corqsice, 0.0)

            # --- here the supersaturation is turned into liquid water
            if supsat > EPSEC:
                if t > RTHOMO:
                    # turn supersaturation into liquid water
                    solqa_ql_qv += supsat
                    solqa_qv_ql -= supsat
                    # include liquid in first guess
                    qlfg += supsat
                else:
                    # turn supersaturation into ice water
                    solqa_qi_qv += supsat
                    solqa_qv_qi -= supsat
                    # add ice to first guess for deposition term
                    qifg += supsat

                # increase cloud amount using RKOOPTAU timescale
                solac = (1 - a) * faci

            # --- include supersaturation from previous timestep
            if in_supsat[0, 0, 0] > EPSEC:
                if t > RTHOMO:
                    # turn supersaturation into liquid water
                    solqa_ql_ql += in_supsat[0, 0, 0]
                    psupsatsrce_ql = in_supsat[0, 0, 0]
                    # add liquid to first guess for deposition term
                    qlfg += in_supsat[0, 0, 0]
                else:
                    # turn supersaturation into ice water
                    solqa_qi_qi += in_supsat[0, 0, 0]
                    psupsatsrce_qi = in_supsat[0, 0, 0]
                    # add ice to first guess for deposition term
                    qifg += in_supsat[0, 0, 0]

                # increase cloud amount using RKOOPTAU timescale
                solac = (1 - a) * faci

            # *** 3.2: detrainment from convection
            if tmp_klevel[0] < NLEV - 1:
                out_lude[0, 0, 0] = in_lude[0, 0, 0] * dtgdp

                if in_convection_on[0, 0] and out_lude[0, 0, 0] > RLMIN and in_lu[0, 0, 1] > EPSEC:
                    solac += out_lude[0, 0, 0] / in_lu[0, 0, 1]
                    # diagnostic temperature split
                    convsrce_ql = out_foealfa[0, 0, 0] * out_lude[0, 0, 0]
                    convsrce_qi = (1 - out_foealfa[0, 0, 0]) * out_lude[0, 0, 0]
                    solqa_ql_ql += convsrce_ql
                    solqa_qi_qi += convsrce_qi
                else:
                    out_lude[0, 0, 0] = 0.0

                # convective snow detrainment source
                if in_convection_on[0, 0]:
                    solqa_qs_qs += in_snde[0, 0, 0] * dtgdp
            else:
                out_lude[0, 0, 0] = in_lude[0, 0, 0]

            # *** 3.3: subsidence compensating convective updraughts
            # --- subsidence source from layer above and evaporation of cloud within the layer
            if tmp_klevel[0] > NCLDTOP - 1:
                mf = max(0.0, (in_mfu + in_mfd) * dtgdp)
                acust = mf * anew[0, 0, -1]

                if __INLINED(not FALLQL and PHASEQL > 0):
                    lcust_ql = mf * out_qln[0, 0, -1]
                    # record total flux for enthalpy budget
                    convsrce_ql += lcust_ql

                if __INLINED(not FALLQI and PHASEQI > 0):
                    lcust_qi = mf * out_qin[0, 0, -1]
                    # record total flux for enthalpy budget
                    convsrce_qi += lcust_qi

                if __INLINED(not FALLQR and PHASEQR > 0):
                    lcust_qr = mf * out_qrn[0, 0, -1]
                    # record total flux for enthalpy budget
                    convsrce_qr += lcust_qr

                if __INLINED(not FALLQS and PHASEQS > 0):
                    lcust_qs = mf * out_qsn[0, 0, -1]
                    # record total flux for enthalpy budget
                    convsrce_qs += lcust_qs

                if __INLINED(not FALLQV and PHASEQV > 0):
                    lcust_qv = mf * qvn[0, 0, -1]
                    # record total flux for enthalpy budget
                    convsrce_qv += lcust_qv

                # work out how much liquid evaporates at arrival point
                dtdp = RDCP * 0.5 * (t[0, 0, -1] + t[0, 0, 0]) / in_aph[0, 0, 0]
                dtforc = dtdp[0, 0, 0] * (in_ap[0, 0, 0] - in_ap[0, 0, -1])
                dqs = anew[0, 0, -1] * dtforc * dqsmixdt

                if __INLINED(not FALLQL and PHASEQL > 0):
                    lfinal = max(0.0, lcust_ql - dqs)
                    evap = min(lcust_ql - lfinal, evaplimmix)
                    lfinal = lcust_ql - evap
                    lfinalsum += lfinal
                    solqa_ql_ql += lcust_ql
                    solqa_qv_ql += evap
                    solqa_ql_qv -= evap

                if __INLINED(not FALLQI and PHASEQI > 0):
                    lfinal = max(0.0, lcust_qi - dqs)
                    evap = min(lcust_qi - lfinal, evaplimmix)
                    lfinal = lcust_qi - evap
                    lfinalsum += lfinal
                    solqa_qi_qi += lcust_qi
                    solqa_qv_qi += evap
                    solqa_qi_qv -= evap

                if __INLINED(not FALLQR and PHASEQR > 0):
                    lfinal = max(0.0, lcust_qr - dqs)
                    evap = min(lcust_qr - lfinal, evaplimmix)
                    lfinal = lcust_qr - evap
                    lfinalsum += lfinal
                    solqa_qr_qr += lcust_qr
                    solqa_qv_qr += evap
                    solqa_qr_qv -= evap

                if __INLINED(not FALLQS and PHASEQS > 0):
                    lfinal = max(0.0, lcust_qs - dqs)
                    evap = min(lcust_qs - lfinal, evaplimmix)
                    lfinal = lcust_qs - evap
                    lfinalsum += lfinal
                    solqa_qs_qs += lcust_qs
                    solqa_qv_qs += evap
                    solqa_qs_qv -= evap

                if __INLINED(not FALLQV and PHASEQV > 0):
                    lfinal = max(0.0, lcust_qv - dqs)
                    evap = min(lcust_qv - lfinal, evaplimmix)
                    lfinal = lcust_qv - evap
                    lfinalsum += lfinal
                    solqa_qv_qv += lcust_qv

                # reset the cloud contribution if no cloud water survives to this level
                if lfinalsum < EPSEC:
                    acust = 0.0
                solac += acust

            # --- subsidence sink of cloud to the layer below
            if tmp_klevel[0] < NLEV - 1:
                mfdn = max(0.0, (in_mfu[0, 0, 1] + in_mfd[0, 0, 1]) * dtgdp)
                solab += mfdn
                solqb_ql_ql += mfdn
                solqb_qi_qi += mfdn

                # record sink for cloud budget and enthalpy budget diagnostics
                convsink_ql = mfdn
                convsink_qi = mfdn

            # *** 3.4: erosion of clouds by turbulent mixing
            # --- define turbulent erosion rate
            ldifdt = RCLDIFF * dt
            if in_convection_type[0, 0] > 0 and out_lude[0, 0, 0] > EPSEC:
                ldifdt *= RCLDIFF_CONVI

            if li > EPSEC:
                # calculate environmental humidity
                e = ldifdt * max(qsmix - qv, 0.0)
                leros = min(min(a * e, evaplimmix), li)
                aeros = leros / licld

                # erosion is -ve linear in L, A
                solac -= aeros
                solqa_qv_ql += liqfrac * leros
                solqa_ql_qv -= liqfrac * leros
                solqa_qv_qi += icefrac * leros
                solqa_qi_qv -= icefrac * leros

            # *** 3.5: condensation/evaporation due to dqsat/dT
            dtdp = RDCP * t / in_ap[0, 0, 0]
            dpmxdt = dp / dt
            mfdn = in_mfu[0, 0, 1] + in_mfd[0, 0, 1] if tmp_klevel[0] < NLEV - 1 else 0.0
            wtot = in_w[0, 0, 0] + 0.5 * RG * (in_mfu[0, 0, 0] + in_mfd[0, 0, 0] + mfdn)
            wtot = min(dpmxdt, max(-dpmxdt, wtot))
            zzdt = in_hrsw[0, 0, 0] + in_hrlw[0, 0, 0]
            dtdiab = min(dpmxdt * dtdp, max(-dpmxdt * dtdp, zzdt)) * dt + RALFDCP * ldefr
            dtforc = dtdp * wtot * dt + dtdiab
            qold = qsmix
            told = t
            t = max(t + dtforc, 160.0)

            qsmix, t = f_cuadjtq(in_ap, qsmix, t)

            dqs = qsmix - qold
            qsmix = qold
            t = told

            # ***: 3.5a: evaporation of clouds
            if dqs > 0:
                levap = min(min(a * min(dqs, licld), evaplimmix), max(qsmix - qv, 0.0))
                solqa_qv_ql += liqfrac * levap
                solqa_ql_qv -= liqfrac * levap
                solqa_qv_qi += icefrac * levap
                solqa_qi_qv -= icefrac * levap

            # *** 3.5b: formation of clouds
            # increase of cloud water in existing clouds
            if a > EPSEC and dqs <= -RLMIN:
                lcond1 = max(-dqs, 0.0)

                # old limiter
                if a > 0.99:
                    cor = 1 / (1 - RETV * qsmix)
                    cdmax = (qv - qsmix) / (1 + cor * qsmix * f_foedem(t))
                else:
                    cdmax = (qv - a * qsmix) / a

                lcond1 = a * max(min(lcond1, cdmax), 0.0)
                if lcond1 < RLMIN:
                    lcond1 = 0.0

                # --- all increase goes into liquid unless so cold cloud homogeneously freezes
                if t > RTHOMO:
                    solqa_ql_qv += lcond1
                    solqa_qv_ql -= lcond1
                    qlfg += lcond1
                else:
                    solqa_qi_qv += lcond1
                    solqa_qv_qi -= lcond1
                    qifg += lcond1

            # generation of new clouds (da/dt > 0)
            if dqs <= -RLMIN and a < 1 - EPSEC:
                # --- critical relative humidity
                rhc = RAMID
                sigk = in_ap[0, 0, 0] / tmp_aph_s[0, 0]
                if sigk > 0.8:
                    rhc += (1 - RAMID) * ((sigk - 0.8) / 0.2) ** 2

                # --- supersaturation options
                if __INLINED(NSSOPT == 0):
                    # no scheme
                    qe = max(0.0, (qv - a * qsice) / max(EPSEC, 1 - a))
                elif __INLINED(NSSOPT == 1):
                    # Tompkins
                    qe = max(0.0, (qv - a * qsice) / max(EPSEC, 1 - a))
                elif __INLINED(NSSOPT == 2):
                    # Lohmann and Karcher
                    qe = qv
                else:
                    # Gierens
                    qe = qv + li

                if t >= RTT or NSSOPT == 0:
                    # no ice supersaturation allowed
                    fac = 1.0
                else:
                    # ice supersaturation
                    fac = fokoop

                if qe >= rhc * qsice * fac and qe < qsice * fac:
                    acond = -(1 - a) * fac * dqs / max(2 * (fac * qsice - qe), EPSEC)
                    acond = min(acond, 1 - a)

                    # linear term
                    lcond2 = -fac * dqs * 0.5 * acond

                    # new limiter formulation
                    zdl = 2 * (fac * qsice - qe) / max(EPSEC, 1 - a)
                    expr2 = fac * dqs
                    if expr2 < -zdl:
                        lcondlim = (a - 1) * expr2 - fac * qsice + qv
                        lcond2 = min(lcond2, lcondlim)
                    lcond2 = max(lcond2, 0.0)

                    expr10 = 1 - a
                    if lcond2 < RLMIN or expr10 < EPSEC:
                        lcond2 = 0.0
                        acond = 0.0
                    if lcond2 == 0.0:
                        acond = 0.0

                    # large-scale generation is linear in A and linear in L
                    solac += acond

                    # --- all increase goes into liquid unless so cold cloud homogeneously freezes
                    if t > RTHOMO:
                        solqa_ql_qv += lcond2
                        solqa_qv_ql -= lcond2
                        qlfg += lcond2
                    else:  # homogeneous freezing
                        solqa_qi_qv += lcond2
                        solqa_qv_qi -= lcond2
                        qifg += lcond2

            # *** 3.6: growth of ice by vapour deposition
            if __INLINED(DEPICE == 1):  # --- ice deposition following Rotstayn et al. (2001)
                # --- calculate distance from cloud top
                if a[0, 0, -1] < RCLDTOPCF and a[0, 0, 0] >= RCLDTOPCF:
                    tmp_cldtopdist[0, 0] = 0.0
                else:
                    tmp_cldtopdist[0, 0] += dp / (rho * RG)

                # --- only treat depositional growth if liquid present
                if t < RTT and qlfg > RLMIN:
                    vpice = f_foeeice(t) * RV / RD
                    vpliq = vpice * fokoop
                    icenuclei = 1000 * exp(12.96 * (vpliq - vpice) / vpliq - 0.639)

                    # --- 0.024 is conductivity of air
                    # --- 8.8 = 700 ** (1/3) = density of ice to the third
                    add = RLSTT * (RLSTT / (RV * t) - 1) / (0.024 * t)
                    bdd = RV * t * in_ap[0, 0, 0] / (2.21 * vpice)
                    cvds = (
                        7.8
                        * (icenuclei / rho) ** 0.666
                        * (vpliq - vpice)
                        / (8.87 * (add + bdd) * vpice)
                    )

                    # --- RICEINIT = 1e-12 is initial mass of ice particle
                    ice0 = max(icecld, icenuclei * RICEINIT / rho)

                    # --- new value of ice
                    inew = (0.666 * cvds * dt + ice0**0.666) ** 1.5

                    # --- grid-mean deposition rate
                    depos = max(a * (inew - ice0), 0.0)

                    # --- limit deposition to liquid water amount
                    depos = min(depos, qlfg)

                    # --- at top of cloud, reduce deposition rate near cloud top
                    infactor = min(icenuclei / 15000, 1.0)
                    depos *= min(
                        infactor
                        + (1 - infactor)
                        * (RDEPLIQREFRATE + tmp_cldtopdist[0, 0] / RDEPLIQREFDEPTH),
                        1.0,
                    )

                    # --- add to matrix
                    solqa_qi_ql += depos
                    solqa_ql_qi -= depos
                    qifg += depos
                    qlfg -= depos
            elif __INLINED(DEPICE == 2):  # --- ice deposition assuming ice PSD
                # --- calculate distance from cloud top
                if a[0, 0, -1] < RCLDTOPCF and a[0, 0, 0] >= RCLDTOPCF:
                    tmp_cldtopdist = 0.0
                else:
                    tmp_cldtopdist += dp / (rho * RG)

                # --- only treat depositional growth if liquid present
                if t < RTT and qlfg > RLMIN:
                    vpice = f_foeeice(t) * RV / RD
                    vpliq = vpice * fokoop
                    icenuclei = 1000 * exp(12.96 * (vpliq - vpice) / vpliq - 0.639)

                    # --- RICEINIT=1e-12 is the initial mass of ice particle
                    ice0 = max(icecld, icenuclei * RICEINIT / rho)

                    # particle size distribution
                    tcg = 1
                    facx1i = 1
                    apb = RCL_APB1 * vpice - RCL_APB2 * vpice * t + in_ap * RCL_APB3 * t**3
                    corrfac = (1 / rho) ** 0.5
                    corrfac2 = ((t / 273) ** 1.5) * 393 / (t + 120)
                    pr02 = rho * ice0 * RCL_CONST1I / (tcg * facx1i)
                    term1 = (
                        (vpliq - vpice)
                        * t**2
                        * vpice
                        * corrfac2
                        * tcg
                        * RCL_CONST2I
                        * facx1i
                        / (rho * apb * vpice)
                    )
                    term2 = (
                        0.65 * RCL_CONST6I * pr02**RCL_CONST4I
                        + RCL_CONST3I
                        * corrfac**0.5
                        * rho**0.5
                        * pr02**RCL_CONST5I
                        / corrfac2**0.5
                    )
                    depos = max(a * term1 * term2 * dt, 0.0)

                    # --- limit deposition to liquid water amount
                    depos = min(depos, qlfg)

                    # --- at top of cloud, reduce deposition rate near cloud top to account for
                    # --- small scale turbulent processes
                    infactor = min(icenuclei / 15000, 1.0)
                    depos *= min(
                        infactor
                        + (1 - infactor) * (RDEPLIQREFRATE + tmp_cldtopdist / RDEPLIQREFDEPTH),
                        1.0,
                    )

                    # --- add to matrix
                    solqa_qi_ql += depos
                    solqa_ql_qi -= depos
                    qifg += depos
                    qlfg -= depos

            # === 4: precipitation processes
            # --- revise in-cloud condensate amount
            tmpa = 1 / max(a, EPSEC)
            liqcld = qlfg * tmpa
            icecld = qifg * tmpa

            # *** 4.1a: sedimentation/falling of ql
            if __INLINED(FALLQL):
                # --- source from layer above
                if tmp_klevel[0] > NCLDTOP - 1:
                    fallsrce_ql = out_pfplsl[0, 0, -1] * dtgdp
                    solqa_ql_ql += fallsrce_ql
                    qlfg += fallsrce_ql
                    # use first guess precip
                    qpretot += qlfg

                # --- sink to next layer, constant fall speed
                fallsink_ql = dtgdp * VQL * rho
            else:
                fallsink_ql = 0.0

            # *** 4.1b: sedimentation/falling of qi
            # --- source from layer above
            if tmp_klevel[0] > NCLDTOP - 1:
                fallsrce_qi = out_pfplsi[0, 0, -1] * dtgdp
                solqa_qi_qi += fallsrce_qi
                qifg += fallsrce_qi
                # use first guess precip
                qpretot += qifg

            # --- sink to next layer, constant fall speed
            if __INLINED(LAERICESED):
                vqi = 0.002 * in_re_ice[0, 0, 0]
            else:
                vqi = VQI
            fallsink_qi = dtgdp * vqi * rho

            # *** 4.1c: sedimentation/falling of qr
            if __INLINED(FALLQR):
                # --- source from layer above
                if tmp_klevel[0] > NCLDTOP - 1:
                    fallsrce_qr = out_pfplsr[0, 0, -1] * dtgdp
                    solqa_qr_qr += fallsrce_qr
                    qrfg += fallsrce_qr
                    # use first guess precip
                    qpretot += qrfg

                # --- sink to next layer, constant fall speed
                fallsink_qr = dtgdp * VQR * rho
            else:
                fallsink_qr = 0.0

            # *** 4.1d: sedimentation/falling of qs
            if __INLINED(FALLQS):
                # --- source from layer above
                if tmp_klevel[0] > NCLDTOP - 1:
                    fallsrce_qs = out_pfplss[0, 0, -1] * dtgdp
                    solqa_qs_qs += fallsrce_qs
                    qsfg += fallsrce_qs
                    # use first guess precip
                    qpretot += qsfg

                # --- sink to next layer, constant fall speed
                fallsink_qs = dtgdp * VQS * rho
            else:
                fallsink_qs = 0.0

            # *** 4.1e: sedimentation/falling of qv
            if __INLINED(FALLQV):
                # --- source from layer above
                if tmp_klevel[0] > NCLDTOP - 1:
                    fallsrce_qv = pfplsv[0, 0, -1] * dtgdp
                    solqa_qv_qv += fallsrce_qv
                    qvfg += fallsrce_qv
                    # use first guess precip
                    qpretot += qvfg

                # --- sink to next layer, constant fall speed
                fallsink_qv = dtgdp * VQV * rho
            else:
                fallsink_qv = 0.0

            # --- precip cover overlap using RAX-RAN Overlap
            if qpretot > EPSEC:
                tmp_covptot[0, 0] = 1 - (
                    (1 - tmp_covptot[0, 0])
                    * (1 - max(a[0, 0, 0], a[0, 0, -1]))
                    / (1 - min(a[0, 0, -1], 1 - 1e-6))
                )
                tmp_covptot[0, 0] = max(tmp_covptot[0, 0], RCOVPMIN)
                covpclr = max(0.0, tmp_covptot[0, 0] - a)
                raincld = qrfg / tmp_covptot[0, 0]
                snowcld = qsfg / tmp_covptot[0, 0]
                tmp_covpmax[0, 0] = max(tmp_covptot[0, 0], tmp_covpmax[0, 0])
            else:
                raincld = 0.0
                snowcld = 0.0
                tmp_covptot[0, 0] = 0.0
                covpclr = 0.0
                tmp_covpmax[0, 0] = 0.0

            # *** 4.2a: autoconversion to snow
            if t <= RTT:
                # --- snow autoconversion rate follow Lin et al. 1983
                if icecld > EPSEC:
                    co = dt * RSNOWLIN1 * exp(RSNOWLIN2 * (t - RTT))

                    if __INLINED(LAERICEAUTO):
                        lcrit = in_icrit_aer[0, 0, 0]
                        co *= (RNICE / in_nice[0, 0, 0]) ** 0.333
                    else:
                        lcrit = RLCRITSNOW

                    snowaut = co * (1 - exp(-((icecld / lcrit) ** 2)))
                    solqb_qs_qi += snowaut

            # *** 4.2b: autoconversion warm clouds
            if liqcld > EPSEC:
                if __INLINED(WARMRAIN == 1):  # --- warm-rain process follow Sundqvist (1989)
                    co = RKCONV * dt

                    if __INLINED(LAERLIQAUTOLSP):
                        lcrit = in_lcrit_aer[0, 0, 0]
                        co *= (RCCN / in_ccn[0, 0, 0]) ** 0.333
                    else:
                        lcrit = RCLCRIT_LAND if in_lsm[0, 0] > 0.5 else RCLCRIT_SEA

                    # --- parameters for cloud collection by rain and snow
                    precip = (out_pfplss[0, 0, -1] + out_pfplsr[0, 0, -1]) / max(
                        EPSEC, tmp_covptot[0, 0]
                    )
                    cfpr = 1 + RPRC1 * sqrt(max(precip, 0.0))
                    if __INLINED(LAERLIQCOLL):
                        cfpr *= (RCCN / in_ccn[0, 0, 0]) ** 0.333

                    co *= cfpr
                    lcrit /= max(cfpr, EPSEC)

                    rainaut = co
                    if liqcld / lcrit < 20:
                        rainaut *= 1 - exp(-((liqcld / lcrit) ** 2))

                    # rain freezes instantly
                    if t <= RTT:
                        solqb_qs_ql += rainaut
                    else:
                        solqb_qr_ql += rainaut
                elif __INLINED(
                    WARMRAIN == 2
                ):  # --- warm-rain process follow Khairoutdinov and Kogan (2000)
                    if in_lsm[0, 0] > 0.5:
                        const = RCL_KK_cloud_num_land
                        lcrit = RCLCRIT_LAND
                    else:
                        const = RCL_KK_cloud_num_sea
                        lcrit = RCLCRIT_SEA

                    if liqcld > lcrit:
                        rainaut = (
                            1.5 * a * dt * RCL_KKAau * liqcld**RCL_KKBauq * const**RCL_KKBaun
                        )
                        rainaut = min(rainaut, qlfg)
                        if rainaut < EPSEC:
                            rainaut = 0.0
                        rainacc = 2 * a * dt * RCL_KKAac * (liqcld * raincld) ** RCL_KKBac
                        rainacc = min(rainacc, qlfg)
                        if rainacc < EPSEC:
                            rainacc = 0.0
                    else:
                        rainaut = 0.0
                        rainacc = 0.0

                    expr3 = rainaut + rainacc
                    if t <= RTT:
                        solqa_qs_ql += expr3
                        solqa_ql_qs -= expr3
                    else:
                        solqa_qr_ql += expr3
                        solqa_ql_qr -= expr3

            # --- riming - collection of cloud liquid drops by snow and ice
            if __INLINED(WARMRAIN > 1):
                if t <= RTT and liqcld > EPSEC:
                    # fallspeed air density correction
                    fallcorr = (RDENSREF / rho) ** 0.4

                    # --- riming of snow by cloud water - implicit in lwc
                    if snowcld > EPSEC and tmp_covptot[0, 0] > 0.01:
                        # calculate riming term
                        snowrime = (
                            0.3
                            * tmp_covptot[0, 0]
                            * dt
                            * RCL_CONST7S
                            * fallcorr
                            * (rho * snowcld * RCL_CONST1S) ** RCL_CONST8S
                        )

                        # limit snow riming term
                        snowrime = min(snowrime, 1.0)

                        solqb_qs_ql += snowrime

            # *** 4.3a: melting of snow and ice
            icetot = qifg + qsfg
            meltmax = 0.0

            # if there are frozen hydrometeors present and dry-bulb temperature > 0degC
            if icetot > EPSEC and t > RTT:
                # calculate subsaturation
                subsat = max(qsice - qv, 0.0)

                # calculate difference between dry-bulb and the temperature at which the wet-buld=0degC
                # using and approx
                tdmtw0 = t - RTT - subsat * (TW1 + TW2 * (in_ap[0, 0, 0] - TW3) - TW4 * (t - TW5))

                # ensure cons1 is positive
                cons1 = abs(dt * (1 + 0.5 * tdmtw0) / RTAUMEL)
                meltmax = max(tdmtw0 * cons1 * RLDCP, 0.0)

            if meltmax > EPSEC and icetot > EPSEC:
                # apply melting in same proportion as frozen hydrometeor fractions
                alfa_qi = qifg / icetot
                melt_qi = min(qifg, alfa_qi * meltmax)
                alfa_qs = qsfg / icetot
                melt_qs = min(qsfg, alfa_qs * meltmax)

                # needed in first guess
                qifg -= melt_qi
                qrfg += melt_qi + melt_qs
                qsfg -= melt_qs
                solqa_qi_qr -= melt_qi
                solqa_qr_qi += melt_qi
                solqa_qr_qs += melt_qs
                solqa_qs_qr -= melt_qs

            # *** 4.3b: freezing of rain
            if qr > EPSEC:
                if t[0, 0, 0] <= RTT and t[0, 0, -1] > RTT:
                    # base of melting layer/top of refreezing layer so store rain/snow fraction for
                    # precip type diagnosis
                    qpretot = max(qs + qr, EPSEC)
                    out_rainfrac_toprfz[0, 0] = qr / qpretot
                    tmp_rainliq[0, 0] = out_rainfrac_toprfz[0, 0] > 0.8

                if t < RTT:
                    if tmp_rainliq[0, 0]:
                        # majority of raindrops completely melted
                        # slope of rain partical size distribution
                        lambda_ = (RCL_FAC1 / (rho * qr)) ** RCL_FAC2

                        # calculate freezing rate based on Bigg (1953) and Wisner (1972)
                        temp = RCL_FZRAB * (t - RTT)
                        frz = dt * (RCL_CONST5R / rho) * (exp(temp) - 1) * lambda_**RCL_CONST6R
                        frzmax = max(frz, 0.0)
                    else:
                        # majority of raindrops only partially melted
                        cons1 = abs(dt * (1 + 0.5 * (RTT - t)) / RTAUMEL)
                        frzmax = max((RTT - t) * cons1 * RLDCP, 0.0)

                    if frzmax > EPSEC:
                        frz = min(qr, frzmax)
                        solqa_qs_qr += frz
                        solqa_qr_qs -= frz

            # *** 4.3c: freezing of liquid
            frzmax = max((RTHOMO - t) * RLDCP, 0.0)
            if frzmax > EPSEC and qlfg > EPSEC:
                frz = min(qlfg, frzmax)
                solqa_qi_ql += frz
                solqa_ql_qi -= frz

            # *** 4.4: evaporation of rain/snow
            if __INLINED(EVAPRAIN == 1):  # --- rain evaporation scheme from Sundquist
                rh = RPRECRHMAX + (1 - RPRECRHMAX) * tmp_covpmax[0, 0] / max(EPSEC, 1 - a)
                rh = min(max(rh, RPRECRHMAX), 1.0)
                qe = (qv - a * qsliq) / max(EPSEC, 1 - a)

                # --- humidity in moistest covpclr part of domain
                qe = max(0.0, min(qe, qsliq))
                lo1 = covpclr > EPSEC and qrfg > EPSEC and qe < rh * qsliq
                if lo1:
                    # note: preclr is a rain flux
                    expr4 = tmp_covptot[0, 0] * dtgdp
                    expr5 = max(abs(expr4), EPSILON)
                    expr6 = expr5 if expr4 > 0 else -expr5
                    preclr = qrfg * covpclr / expr6

                    # --- actual microphysics formula in beta
                    beta1 = (
                        sqrt(in_ap[0, 0, 0] / tmp_aph_s[0, 0])
                        / RVRFACTOR
                        * preclr
                        / max(covpclr, EPSEC)
                    )
                    beta = RG * RPECONS * 0.5 * beta1**0.5777
                    denom = 1 + beta * dt * corqsliq
                    dpr = covpclr * beta * (qsliq - qe) / denom * dp / RG
                    dpevap = dpr * dtgdp

                    # --- add evaporation term to explicit sink
                    evap = min(dpevap, qrfg)
                    solqa_qv_qr += evap
                    solqa_qr_qv -= evap

                    # --- reduce the total precip coverage proportional to evaporation
                    tmp_covptot[0, 0] = max(
                        RCOVPMIN,
                        tmp_covptot[0, 0] - max(0.0, (tmp_covptot[0, 0] - a) * evap / qrfg),
                    )

                    # update fg field
                    qrfg -= evap
            elif __INLINED(
                EVAPRAIN == 2
            ):  # --- rain evaporation scheme based on Abel and Boutle (2013)
                # --- calculate relative humidity limit for rain evaporation
                # limit rh for rain evaporation dependent on precipitation fraction
                rh = RPRECRHMAX + (1 - RPRECRHMAX) * tmp_covpmax[0, 0] / max(EPSEC, 1 - a)
                rh = min(max(rh, RPRECRHMAX), 1.0)

                # further limit rh for rain evaporation to 80%
                rh = min(0.8, rh)

                qe = max(0.0, min(qv, qsliq))
                lo1 = covpclr > EPSEC and qrfg > EPSEC and qe < rh * qsliq
                if lo1:
                    # --- Abel and Boutle (2012) evaporation
                    # calculate local precipitation (kg/kg)
                    preclr = qrfg / tmp_covptot[0, 0]

                    # fallspeed air density correction
                    fallcorr = (RDENSREF / rho) ** 0.4

                    # saturation vapor pressure with respect to liquid phase
                    esatliq = RV / RD * f_foeeliq(t)

                    # slope of particle size distribution
                    lambda_ = (RCL_FAC1 / (rho * preclr)) ** RCL_FAC2

                    evap_denom = (
                        RCL_CDENOM1 * esatliq
                        - RCL_CDENOM2 * t * esatliq
                        + RCL_CDENOM3 * t**3 * in_ap[0, 0, 0]
                    )

                    # temperature dependent conductivity
                    corr2 = (t / 273) ** 1.5 * 393 / (t + 120)

                    subsat = max(rh * qsliq - qe, 0.0)
                    beta = (
                        0.5
                        / qsliq
                        * t**2
                        * esatliq
                        * RCL_CONST1R
                        * (corr2 / evap_denom)
                        * (
                            0.78 / lambda_**RCL_CONST4R
                            + RCL_CONST2R
                            * (rho * fallcorr) ** 0.5
                            / (corr2**0.5 * lambda_**RCL_CONST3R)
                        )
                    )
                    denom = 1 + beta * dt
                    dpevap = covpclr * beta * dt * subsat / denom

                    # --- add evaporation term to explicit sink
                    evap = min(dpevap, qrfg)
                    solqa_qv_qr += evap
                    solqa_qr_qv -= evap

                    # --- reduce the total precip coverage proportional to evaporation
                    tmp_covptot[0, 0] = max(
                        RCOVPMIN,
                        tmp_covptot[0, 0] - max(0.0, (tmp_covptot[0, 0] - a) * evap / qrfg),
                    )

                    # update fg field
                    qrfg -= evap

            # *** 4.5: evaporation of snow
            if __INLINED(EVAPSNOW == 1):
                rh = RPRECRHMAX + (1 - RPRECRHMAX) * tmp_covpmax[0, 0] / max(EPSEC, 1 - a)
                rh = min(max(rh, RPRECRHMAX), 1.0)
                qe = (qv - a * qsice) / max(EPSEC, 1 - a)

                # --- humidity in moistest covpclr part of domain
                qe = max(0.0, min(qe, qsice))
                lo1 = covpclr > EPSEC and qsfg > EPSEC and qe < rh * qsice
                if lo1:
                    expr7 = tmp_covptot[0, 0] * dtgdp
                    expr8 = max(abs(expr7), EPSILON)
                    expr9 = expr8 if expr7 > 0 else -expr8
                    preclr = qsfg * covpclr / expr9

                    # --- actual microphysics formula in beta
                    beta1 = (
                        sqrt(in_ap[0, 0, 0] / tmp_aph_s[0, 0])
                        / RVRFACTOR
                        * preclr
                        / max(covpclr, EPSEC)
                    )
                    beta = RG * RPECONS * beta1**0.5777
                    denom = 1 + beta * dt * corqsice
                    dpr = covpclr * beta * (qsice - qe) / denom * dp / RG
                    dpevap = dpr * dtgdp

                    # --- add evaporation term to explicit sink
                    evap = min(dpevap, qsfg)
                    solqa_qv_qs += evap
                    solqa_qs_qv -= evap

                    # --- reduce the total precip coverage proportional to evaporation
                    tmp_covptot[0, 0] = max(
                        RCOVPMIN,
                        tmp_covptot[0, 0] - max(0.0, (tmp_covptot[0, 0] - a) * evap / qsfg),
                    )

                    # update first guess field
                    qsfg -= evap
            elif __INLINED(EVAPSNOW == 2):
                # --- calculate relative humidity limit for snow evaporation
                rh = RPRECRHMAX + (1 - RPRECRHMAX) * tmp_covpmax[0, 0] / max(EPSEC, 1 - a)
                rh = min(max(rh, RPRECRHMAX), 1.0)
                qe = (qv - a * qsice) / max(EPSEC, 1 - a)

                # --- humidity in moistest covpclr part of domain
                qe = max(0.0, min(qe, qsice))
                lo1 = covpclr > EPSEC and qs > EPSEC and qe < rh * qsice
                if lo1:
                    # calculate local precipitation (kg/kg)
                    preclr = qsfg / tmp_covptot[0, 0]
                    vpice = f_foeeice(t) * RV / RD

                    # particle size distribution
                    tcg = 1.0
                    facx1s = 1.0
                    apb = (
                        RCL_APB1 * vpice - RCL_APB2 * vpice * t + in_ap[0, 0, 0] * RCL_APB3 * t**3
                    )
                    corrfac = (1 / rho) ** 0.5
                    corrfac2 = ((t / 273) ** 1.5) * 393 / (t + 120)
                    pr02 = rho * preclr * RCL_CONST1S / (tcg * facx1s)
                    term1 = (
                        (qsice - qe)
                        * t**2
                        * vpice
                        * corrfac2
                        * tcg
                        * RCL_CONST2S
                        * facx1s
                        / (rho * apb * qsice)
                    )
                    term2 = (
                        0.65 * RCL_CONST6S * pr02**RCL_CONST4S
                        + RCL_CONST3S
                        * corrfac**0.5
                        * rho**0.5
                        * pr02**RCL_CONST5S
                        / corrfac2**0.5
                    )
                    dpevap = max(covpclr * term1 * term2 * dt, 0.0)

                    # --- limit evaporation to snow amount
                    evap = min(min(dpevap, evaplimice), qs)
                    solqa_qv_qs += evap
                    solqa_qs_qv -= evap

                    # --- reduce the total precip coverage proportional to evaporation
                    tmp_covptot[0, 0] = max(
                        RCOVPMIN, tmp_covptot[0, 0] - max(0.0, (tmp_covptot[0, 0] - a) * evap / qs)
                    )

                    # update first guess field
                    qsfg -= evap

            # --- evaporate small precipitation amounts
            if __INLINED(FALLQL):
                if qlfg < RLMIN:
                    solqa_qv_ql += qlfg
                    solqa_ql_qv -= qlfg
            if __INLINED(FALLQI):
                if qifg < RLMIN:
                    solqa_qv_qi += qifg
                    solqa_qi_qv -= qifg
            if __INLINED(FALLQR):
                if qrfg < RLMIN:
                    solqa_qv_qr += qrfg
                    solqa_qr_qv -= qrfg
            if __INLINED(FALLQS):
                if qsfg < RLMIN:
                    solqa_qv_qs += qsfg
                    solqa_qs_qv -= qsfg

            # === 5: solvers for A and L
            # *** 5.1: solver for cloud cover
            anew = min((a + solac) / (1 + solab), 1.0)
            if anew < RAMIN:
                anew = 0.0
            da = anew - a0

            # *** 5.2: solver for the microphysics
            # --- collect sink terms and mark
            sinksum_ql = -(solqa_ql_ql + solqa_ql_qi + solqa_ql_qr + solqa_ql_qs + solqa_ql_qv)
            sinksum_qi = -(solqa_qi_ql + solqa_qi_qi + solqa_qi_qr + solqa_qi_qs + solqa_qi_qv)
            sinksum_qr = -(solqa_qr_ql + solqa_qr_qi + solqa_qr_qr + solqa_qr_qs + solqa_qr_qv)
            sinksum_qs = -(solqa_qs_ql + solqa_qs_qi + solqa_qs_qr + solqa_qs_qs + solqa_qs_qv)
            sinksum_qv = -(solqa_qv_ql + solqa_qv_qi + solqa_qv_qr + solqa_qv_qs + solqa_qv_qv)

            # --- calculate overshoot and scaling factor
            max_ql = max(ql, EPSEC)
            rat_ql = max(sinksum_ql, max_ql)
            ratio_ql = max_ql / rat_ql
            max_qi = max(qi, EPSEC)
            rat_qi = max(sinksum_qi, max_qi)
            ratio_qi = max_qi / rat_qi
            max_qr = max(qr, EPSEC)
            rat_qr = max(sinksum_qr, max_qr)
            ratio_qr = max_qr / rat_qr
            max_qs = max(qs, EPSEC)
            rat_qs = max(sinksum_qs, max_qs)
            ratio_qs = max_qs / rat_qs
            max_qv = max(qv, EPSEC)
            rat_qv = max(sinksum_qv, max_qv)
            ratio_qv = max_qv / rat_qv

            # --- now sort ratio to find out which species run out first
            order_ql, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv = f_helper_0(
                order_ql,
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
            )
            order_qi, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv = f_helper_0(
                order_qi,
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
            )
            order_qr, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv = f_helper_0(
                order_qr,
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
            )
            order_qs, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv = f_helper_0(
                order_qs,
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
            )
            order_qv, index1_ql, index1_qi, index1_qr, index1_qs, index1_qv = f_helper_0(
                order_qv,
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
            )

            # scale the sink terms, in the correct order, recalculating the scale factor each time
            sinksum_ql = 0.0
            sinksum_qi = 0.0
            sinksum_qr = 0.0
            sinksum_qs = 0.0
            sinksum_qv = 0.0

            # --- recalculate sum and scaling factor, and then scale
            ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv = f_helper_1(
                order_ql,
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
            )
            ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv = f_helper_1(
                order_qi,
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
            )
            ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv = f_helper_1(
                order_qr,
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
            )
            ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv = f_helper_1(
                order_qs,
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
            )
            ratio_ql, ratio_qi, ratio_qr, ratio_qs, ratio_qv = f_helper_1(
                order_qv,
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
            )

            # *** 5.2.2: solver
            # --- set the lhs of equation
            # --- diagonals: microphysical sink terms + transport
            lhs_ql_ql = (
                1
                + fallsink_ql
                + solqb_qv_ql
                + solqb_ql_ql
                + solqb_qi_ql
                + solqb_qr_ql
                + solqb_qs_ql
            )
            lhs_qi_qi = (
                1
                + fallsink_qi
                + solqb_qv_qi
                + solqb_ql_qi
                + solqb_qi_qi
                + solqb_qr_qi
                + solqb_qs_qi
            )
            lhs_qr_qr = (
                1
                + fallsink_qr
                + solqb_qv_qr
                + solqb_ql_qr
                + solqb_qi_qr
                + solqb_qr_qr
                + solqb_qs_qr
            )
            lhs_qs_qs = (
                1
                + fallsink_qs
                + solqb_qv_qs
                + solqb_ql_qs
                + solqb_qi_qs
                + solqb_qr_qs
                + solqb_qs_qs
            )
            lhs_qv_qv = (
                1
                + fallsink_qv
                + solqb_qv_qv
                + solqb_ql_qv
                + solqb_qi_qv
                + solqb_qr_qv
                + solqb_qs_qv
            )

            # --- non-diagonals: microphysical source terms
            lhs_ql_qi = -solqb_ql_qi
            lhs_ql_qr = -solqb_ql_qr
            lhs_ql_qs = -solqb_ql_qs
            lhs_ql_qv = -solqb_ql_qv
            lhs_qi_ql = -solqb_qi_ql
            lhs_qi_qr = -solqb_qi_qr
            lhs_qi_qs = -solqb_qi_qs
            lhs_qi_qv = -solqb_qi_qv
            lhs_qr_ql = -solqb_qr_ql
            lhs_qr_qi = -solqb_qr_qi
            lhs_qr_qs = -solqb_qr_qs
            lhs_qr_qv = -solqb_qr_qv
            lhs_qs_ql = -solqb_qs_ql
            lhs_qs_qi = -solqb_qs_qi
            lhs_qs_qr = -solqb_qs_qr
            lhs_qs_qv = -solqb_qs_qv
            lhs_qv_ql = -solqb_qv_ql
            lhs_qv_qi = -solqb_qv_qi
            lhs_qv_qr = -solqb_qv_qr
            lhs_qv_qs = -solqb_qv_qs

            # --- set the rhs of equation
            # --- sum the explicit source and sink
            out_qln[0, 0, 0] = (
                ql + solqa_ql_ql + solqa_ql_qi + solqa_ql_qr + solqa_ql_qs + solqa_ql_qv
            )
            out_qin[0, 0, 0] = (
                qi + solqa_qi_ql + solqa_qi_qi + solqa_qi_qr + solqa_qi_qs + solqa_qi_qv
            )
            out_qrn[0, 0, 0] = (
                qr + solqa_qr_ql + solqa_qr_qi + solqa_qr_qr + solqa_qr_qs + solqa_qr_qv
            )
            out_qsn[0, 0, 0] = (
                qs + solqa_qs_ql + solqa_qs_qi + solqa_qs_qr + solqa_qs_qs + solqa_qs_qv
            )
            qvn = qv + solqa_qv_ql + solqa_qv_qi + solqa_qv_qr + solqa_qv_qs + solqa_qv_qv

            # --- solve by LU decomposition
            # non pivoting recursive factorization
            lhs_qi_ql /= lhs_ql_ql  #             JN=1, JM=2
            lhs_qi_qi -= lhs_qi_ql * lhs_ql_qi  # JN=1, JM=2, IK=2
            lhs_qi_qr -= lhs_qi_ql * lhs_ql_qr  # JN=1, JM=2, IK=3
            lhs_qi_qs -= lhs_qi_ql * lhs_ql_qs  # JN=1, JM=2, IK=4
            lhs_qi_qv -= lhs_qi_ql * lhs_ql_qv  # JN=1, JM=2, IK=0
            lhs_qr_ql /= lhs_ql_ql  #             JN=1, JM=3
            lhs_qr_qi -= lhs_qr_ql * lhs_ql_qi  # JN=1, JM=3, IK=2
            lhs_qr_qr -= lhs_qr_ql * lhs_ql_qr  # JN=1, JM=3, IK=3
            lhs_qr_qs -= lhs_qr_ql * lhs_ql_qs  # JN=1, JM=3, IK=4
            lhs_qr_qv -= lhs_qr_ql * lhs_ql_qv  # JN=1, JM=3, IK=0
            lhs_qs_ql /= lhs_ql_ql  #             JN=1, JM=4
            lhs_qs_qi -= lhs_qs_ql * lhs_ql_qi  # JN=1, JM=4, IK=2
            lhs_qs_qr -= lhs_qs_ql * lhs_ql_qr  # JN=1, JM=4, IK=3
            lhs_qs_qs -= lhs_qs_ql * lhs_ql_qs  # JN=1, JM=4, IK=4
            lhs_qs_qv -= lhs_qs_ql * lhs_ql_qv  # JN=1, JM=4, IK=0
            lhs_qv_ql /= lhs_ql_ql  #             JN=1, JM=0
            lhs_qv_qi -= lhs_qv_ql * lhs_ql_qi  # JN=1, JM=0, IK=2
            lhs_qv_qr -= lhs_qv_ql * lhs_ql_qr  # JN=1, JM=0, IK=3
            lhs_qv_qs -= lhs_qv_ql * lhs_ql_qs  # JN=1, JM=0, IK=4
            lhs_qv_qv -= lhs_qv_ql * lhs_ql_qv  # JN=1, JM=0, IK=0
            lhs_qr_qi /= lhs_qi_qi  #             JN=2, JM=3
            lhs_qr_qr -= lhs_qr_qi * lhs_qi_qr  # JN=2, JM=3, IK=3
            lhs_qr_qs -= lhs_qr_qi * lhs_qi_qs  # JN=2, JM=3, IK=4
            lhs_qr_qv -= lhs_qr_qi * lhs_qi_qv  # JN=2, JM=3, IK=0
            lhs_qs_qi /= lhs_qi_qi  #             JN=2, JM=4
            lhs_qs_qr -= lhs_qs_qi * lhs_qi_qr  # JN=2, JM=4, IK=3
            lhs_qs_qs -= lhs_qs_qi * lhs_qi_qs  # JN=2, JM=4, IK=4
            lhs_qs_qv -= lhs_qs_qi * lhs_qi_qv  # JN=2, JM=4, IK=0
            lhs_qv_qi /= lhs_qi_qi  #             JN=2, JM=0
            lhs_qv_qr -= lhs_qv_qi * lhs_qi_qr  # JN=2, JM=0, IK=3
            lhs_qv_qs -= lhs_qv_qi * lhs_qi_qs  # JN=2, JM=0, IK=4
            lhs_qv_qv -= lhs_qv_qi * lhs_qi_qv  # JN=2, JM=0, IK=0
            lhs_qs_qr /= lhs_qr_qr  #             JN=3, JM=4
            lhs_qs_qs -= lhs_qs_qr * lhs_qr_qs  # JN=3, JM=4, IK=4
            lhs_qs_qv -= lhs_qs_qr * lhs_qr_qv  # JN=3, JM=4, IK=0
            lhs_qv_qr /= lhs_qr_qr  #             JN=3, JM=0
            lhs_qv_qs -= lhs_qv_qr * lhs_qr_qs  # JN=3, JM=0, IK=4
            lhs_qv_qv -= lhs_qv_qr * lhs_qr_qv  # JN=3, JM=0, IK=0
            lhs_qv_qs /= lhs_qs_qs  #             JN=4, JM=0
            lhs_qv_qv -= lhs_qv_qs * lhs_qs_qv  # JN=4, JM=0, IK=0

            # backsubstitution: step 1
            out_qin[0, 0, 0] -= lhs_qi_ql * out_qln[0, 0, 0]
            out_qrn[0, 0, 0] -= lhs_qr_ql * out_qln[0, 0, 0] + lhs_qr_qi * out_qin[0, 0, 0]
            out_qsn[0, 0, 0] -= (
                lhs_qs_ql * out_qln[0, 0, 0]
                + lhs_qs_qi * out_qin[0, 0, 0]
                + lhs_qs_qr * out_qrn[0, 0, 0]
            )
            qvn -= (
                lhs_qv_ql * out_qln[0, 0, 0]
                + lhs_qv_qi * out_qin[0, 0, 0]
                + lhs_qv_qr * out_qrn[0, 0, 0]
                + lhs_qv_qs * out_qsn[0, 0, 0]
            )

            # backsubstitution: step 2
            qvn /= lhs_qv_qv
            out_qsn[0, 0, 0] -= lhs_qs_qv * qvn
            out_qsn[0, 0, 0] /= lhs_qs_qs
            out_qrn[0, 0, 0] -= lhs_qr_qs * out_qsn[0, 0, 0] + lhs_qr_qv * qvn
            out_qrn[0, 0, 0] /= lhs_qr_qr
            out_qin[0, 0, 0] -= (
                lhs_qi_qr * out_qrn[0, 0, 0] + lhs_qi_qs * out_qsn[0, 0, 0] + lhs_qi_qv * qvn
            )
            out_qin[0, 0, 0] /= lhs_qi_qi
            out_qln[0, 0, 0] -= (
                lhs_ql_qi * out_qin[0, 0, 0]
                + lhs_ql_qr * out_qrn[0, 0, 0]
                + lhs_ql_qs * out_qsn[0, 0, 0]
                + lhs_ql_qv * qvn
            )
            out_qln[0, 0, 0] /= lhs_ql_ql

            # ensure no small values (including negatives) remain in cloud variables
            # nor precipitation rates
            if out_qln[0, 0, 0] < EPSEC:
                qvn += out_qln[0, 0, 0]
                out_qln[0, 0, 0] = 0.0
            if out_qin[0, 0, 0] < EPSEC:
                qvn += out_qin[0, 0, 0]
                out_qin[0, 0, 0] = 0.0
            if out_qrn[0, 0, 0] < EPSEC:
                qvn += out_qrn[0, 0, 0]
                out_qrn[0, 0, 0] = 0.0
            if out_qsn[0, 0, 0] < EPSEC:
                qvn += out_qsn[0, 0, 0]
                out_qsn[0, 0, 0] = 0.0

            # *** 5.3: precipitation/sedimentation fluxes to next level diagnostic precipitation fluxes
            out_pfplsl[0, 0, 0] = fallsink_ql * out_qln[0, 0, 0] * rdtgdp
            out_pfplsi[0, 0, 0] = fallsink_qi * out_qin[0, 0, 0] * rdtgdp
            out_pfplsr[0, 0, 0] = fallsink_qr * out_qrn[0, 0, 0] * rdtgdp
            out_pfplss[0, 0, 0] = fallsink_qs * out_qsn[0, 0, 0] * rdtgdp
            pfplsv = fallsink_qv * qvn * rdtgdp

            # ensure precipitation fraction is zero if no precipitation
            qpretot = out_pfplss[0, 0, 0] + out_pfplsr[0, 0, 0]
            if qpretot < EPSEC:
                tmp_covptot[0, 0] = 0.0

            # === 6: update tendencies
            # *** 6.1: temperature and CLV budgets
            flux_ql = (
                psupsatsrce_ql
                + convsrce_ql
                + fallsrce_ql
                - (fallsink_ql + convsink_ql) * out_qln[0, 0, 0]
            )
            if __INLINED(PHASEQL == 1):
                out_tnd_loc_t[0, 0, 0] += RALVDCP * (out_qln[0, 0, 0] - ql - flux_ql) / dt
            if __INLINED(PHASEQL == 2):
                out_tnd_loc_t[0, 0, 0] += RALSDCP * (out_qln[0, 0, 0] - ql - flux_ql) / dt
            out_tnd_loc_ql[0, 0, 0] += (out_qln[0, 0, 0] - out_ql0[0, 0, 0]) / dt

            flux_qi = (
                psupsatsrce_qi
                + convsrce_qi
                + fallsrce_qi
                - (fallsink_qi + convsink_qi) * out_qin[0, 0, 0]
            )
            if __INLINED(PHASEQI == 1):
                out_tnd_loc_t[0, 0, 0] += RALVDCP * (out_qin[0, 0, 0] - qi - flux_qi) / dt
            if __INLINED(PHASEQI == 2):
                out_tnd_loc_t[0, 0, 0] += RALSDCP * (out_qin[0, 0, 0] - qi - flux_qi) / dt
            out_tnd_loc_qi[0, 0, 0] += (out_qin[0, 0, 0] - out_qi0[0, 0, 0]) / dt

            flux_qr = (
                psupsatsrce_qr
                + convsrce_qr
                + fallsrce_qr
                - (fallsink_qr + convsink_qr) * out_qrn[0, 0, 0]
            )
            if __INLINED(PHASEQR == 1):
                out_tnd_loc_t[0, 0, 0] += RALVDCP * (out_qrn[0, 0, 0] - qr - flux_qr) / dt
            if __INLINED(PHASEQR == 2):
                out_tnd_loc_t[0, 0, 0] += RALSDCP * (out_qrn[0, 0, 0] - qr - flux_qr) / dt
            out_tnd_loc_qr[0, 0, 0] += (out_qrn[0, 0, 0] - out_qr0[0, 0, 0]) / dt

            flux_qs = (
                psupsatsrce_qs
                + convsrce_qs
                + fallsrce_qs
                - (fallsink_qs + convsink_qs) * out_qsn[0, 0, 0]
            )
            if __INLINED(PHASEQS == 1):
                out_tnd_loc_t[0, 0, 0] += RALVDCP * (out_qsn[0, 0, 0] - qs - flux_qs) / dt
            if __INLINED(PHASEQS == 2):
                out_tnd_loc_t[0, 0, 0] += RALSDCP * (out_qsn[0, 0, 0] - qs - flux_qs) / dt
            out_tnd_loc_qs[0, 0, 0] += (out_qsn[0, 0, 0] - out_qs0[0, 0, 0]) / dt

            # *** 6.2: humidity budget
            out_tnd_loc_qv[0, 0, 0] += (qvn - qv) / dt

            # *** 6.3: cloud cover
            out_tnd_loc_a[0, 0, 0] += da / dt

            # --- copy precipitation fraction into output variable
            out_covptot[0, 0, 0] = tmp_covptot[0, 0]


@stencil_collection("cloudsc_fluxes")
def cloudsc_fluxes(
    in_aph: Field["float"],  # staggered
    in_foealfa: Field["float"],
    in_lneg_qi: Field["float"],
    in_lneg_ql: Field["float"],
    in_lneg_qr: Field["float"],
    in_lneg_qs: Field["float"],
    in_lude: Field["float"],
    in_pfplsi: Field["float"],
    in_pfplsl: Field["float"],
    in_pfplsr: Field["float"],
    in_pfplss: Field["float"],
    in_qi0: Field["float"],
    in_qin: Field["float"],
    in_ql0: Field["float"],
    in_qln: Field["float"],
    in_qr0: Field["float"],
    in_qrn: Field["float"],
    in_qs0: Field["float"],
    in_qsn: Field["float"],
    in_vfi: Field["float"],
    in_vfl: Field["float"],
    out_fcqlng: Field["float"],  # staggered
    out_fcqnng: Field["float"],  # staggered
    out_fcqrng: Field["float"],  # staggered
    out_fcqsng: Field["float"],  # staggered
    out_fhpsl: Field["float"],  # staggered
    out_fhpsn: Field["float"],  # staggered
    out_fplsl: Field["float"],  # staggered
    out_fplsn: Field["float"],  # staggered
    out_fsqif: Field["float"],  # staggered
    out_fsqitur: Field["float"],  # staggered
    out_fsqlf: Field["float"],  # staggered
    out_fsqltur: Field["float"],  # staggered
    out_fsqrf: Field["float"],  # staggered
    out_fsqsf: Field["float"],  # staggered
    *,
    dt: "float",
):
    from __externals__ import RG, RLSTT, RLVTT

    # === 7: flux/diagnostics computations
    with computation(FORWARD):
        with interval(0, 1):
            out_fplsl[0, 0, 0] = 0.0
            out_fplsn[0, 0, 0] = 0.0
            out_fhpsl[0, 0, 0] = 0.0
            out_fhpsn[0, 0, 0] = 0.0
            out_fsqlf[0, 0, 0] = 0.0
            out_fsqif[0, 0, 0] = 0.0
            out_fsqrf[0, 0, 0] = 0.0
            out_fsqsf[0, 0, 0] = 0.0
            out_fcqlng[0, 0, 0] = 0.0
            out_fcqnng[0, 0, 0] = 0.0
            out_fcqrng[0, 0, 0] = 0.0
            out_fcqsng[0, 0, 0] = 0.0
            out_fsqltur[0, 0, 0] = 0.0
            out_fsqitur[0, 0, 0] = 0.0

        with interval(1, None):
            # --- copy general precip arrays back info PFP arrays for GRIB archiving
            out_fplsl[0, 0, 0] = in_pfplsr[0, 0, -1] + in_pfplsl[0, 0, -1]
            out_fplsn[0, 0, 0] = in_pfplss[0, 0, -1] + in_pfplsi[0, 0, -1]

            # --- enthalpy flux due to precipitation
            out_fhpsl[0, 0, 0] = -RLVTT * out_fplsl[0, 0, 0]
            out_fhpsn[0, 0, 0] = -RLSTT * out_fplsn[0, 0, 0]

            gdph_r = -(in_aph[0, 0, 0] - in_aph[0, 0, -1]) / (RG * dt)
            out_fsqlf[0, 0, 0] = out_fsqlf[0, 0, -1]
            out_fsqif[0, 0, 0] = out_fsqif[0, 0, -1]
            out_fsqrf[0, 0, 0] = out_fsqlf[0, 0, -1]
            out_fsqsf[0, 0, 0] = out_fsqif[0, 0, -1]
            out_fcqlng[0, 0, 0] = out_fcqlng[0, 0, -1]
            out_fcqnng[0, 0, 0] = out_fcqnng[0, 0, -1]
            out_fcqrng[0, 0, 0] = out_fcqlng[0, 0, -1]
            out_fcqsng[0, 0, 0] = out_fcqnng[0, 0, -1]
            out_fsqltur[0, 0, 0] = out_fsqltur[0, 0, -1]
            out_fsqitur[0, 0, 0] = out_fsqitur[0, 0, -1]

            # liquid, LS scheme minus detrainment
            out_fsqlf[0, 0, 0] += (
                in_qln[0, 0, -1]
                - in_ql0[0, 0, -1]
                + in_vfl[0, 0, -1] * dt
                - in_foealfa[0, 0, -1] * in_lude[0, 0, -1]
            ) * gdph_r
            # liquid, negative numbers
            out_fcqlng[0, 0, 0] += in_lneg_ql[0, 0, -1] * gdph_r
            # liquid, vertical diffusion
            out_fsqltur[0, 0, 0] += in_vfl[0, 0, -1] * dt * gdph_r

            # rain, LS scheme
            out_fsqrf[0, 0, 0] += (in_qrn[0, 0, -1] - in_qr0[0, 0, -1]) * gdph_r
            # rain, negative numbers
            out_fcqrng[0, 0, 0] += in_lneg_qr[0, 0, -1] * gdph_r

            # ice, LS scheme minus detrainment
            out_fsqif[0, 0, 0] += (
                in_qin[0, 0, -1]
                - in_qi0[0, 0, -1]
                + in_vfi[0, 0, -1] * dt
                - (1 - in_foealfa[0, 0, -1]) * in_lude[0, 0, -1]
            ) * gdph_r
            # ice, negative numbers
            out_fcqnng[0, 0, 0] += in_lneg_qi[0, 0, -1] * gdph_r
            # ice, vertical diffusion
            out_fsqitur[0, 0, 0] += in_vfi[0, 0, -1] * dt * gdph_r

            # snow, LS scheme
            out_fsqsf[0, 0, 0] += (in_qsn[0, 0, -1] - in_qs0[0, 0, -1]) * gdph_r
            # snow, negative numbers
            out_fcqsng[0, 0, 0] += in_lneg_qs[0, 0, -1] * gdph_r
