/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "cuda.h"
#include "yoecldp_c.h"
//struct TECLDP ;
#include <stdio.h>

__global__ void cloudsc_c(int kidia, int kfdia, int klon, double ptsphy, double *  v_pt, double *  v_pq,
	      double *  v_tendency_tmp_t, double *  v_tendency_tmp_q, double *  v_tendency_tmp_a, double *  v_tendency_tmp_cld,
	      double *  v_tendency_loc_t, double *  v_tendency_loc_q, double *  v_tendency_loc_a, double *  v_tendency_loc_cld,
	      double *  v_pvfa, double *  v_pvfl, double *  v_pvfi, double *  v_pdyna, double *  v_pdynl, double *  v_pdyni,
	      double *  v_phrsw, double *  v_phrlw, double *  v_pvervel, double *  v_pap, double *  v_paph, double *  v_plsm,
	      int *  v_ldcum,
	      int *  v_ktype, double *  v_plu, double *  v_plude, double *  v_psnde, double *  v_pmfu,
	      double *  v_pmfd, double *  v_pa, double *  v_pclv, double *  v_psupsat, double *  v_plcrit_aer, double *  v_picrit_aer,
	      double *  v_pre_ice, double *  v_pccn, double *  v_pnice, double *  v_pcovptot, double *  v_prainfrac_toprfz, double *  v_pfsqlf,
	      double *  v_pfsqif, double *  v_pfcqnng, double *  v_pfcqlng, double *  v_pfsqrf, double *  v_pfsqsf, double *  v_pfcqrng,
	      double *  v_pfcqsng, double *  v_pfsqltur, double *  v_pfsqitur, double *  v_pfplsl, double *  v_pfplsn, double *  v_pfhpsl,
	      double *  v_pfhpsn,
	      struct TECLDP *yrecldp, int nblocks,
              double rg, double rd, double rcpd, double retv, double rlvtt, double rlstt, double rlmlt, double rtt,
              double rv, double r2es, double r3les, double r3ies, double r4les, double r4ies, double r5les,
              double r5ies, double r5alvcp, double r5alscp, double ralvdcp, double ralsdcp, double ralfdcp,
              double rtwat, double rtice, double rticecu, double rtwat_rtice_r, double rtwat_rticecu_r,
              double rkoop1, double rkoop2,
	      double *zfoealfa, double *ztp1, double *zli,
              double *za, double *zaorig, double *zliqfrac,
              double *zicefrac, double *zqx, double *zqx0,
              double *zpfplsx, double *zlneg, double *zqxn2d,
              double *zqsmix, double *zqsliq, double *zqsice,
              double *zfoeewmt, double *zfoeew, double *zfoeeliqt);
