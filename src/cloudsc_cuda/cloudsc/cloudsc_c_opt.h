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
#include <stdio.h>
#include "dtype.h"

__global__ void cloudsc_c_opt(int kidia, int kfdia, int klon, dtype ptsphy,
  const dtype * __restrict__  pt,
  const dtype * __restrict__  pq, const dtype * __restrict__  tendency_tmp_t,
  const dtype * __restrict__  tendency_tmp_q, const dtype * __restrict__  tendency_tmp_a,
  const dtype * __restrict__  tendency_tmp_cld, dtype * __restrict__  tendency_loc_t,
  dtype * __restrict__  tendency_loc_q, dtype * __restrict__  tendency_loc_a,
  dtype * __restrict__  tendency_loc_cld, const dtype * __restrict__  pvfa,
  const dtype * __restrict__  pvfl, const dtype * __restrict__  pvfi, const dtype * __restrict__  pdyna,
  const dtype * __restrict__  pdynl, const dtype * __restrict__  pdyni, const dtype * __restrict__  phrsw,
  dtype * __restrict__  phrlw, const dtype * __restrict__  pvervel, const dtype * __restrict__  pap,
  const dtype * __restrict__  paph, const dtype * __restrict__  plsm,
  const int *  ktype, const dtype * __restrict__  plu, dtype * __restrict__  plude,
  const dtype * __restrict__  psnde, const dtype * __restrict__  pmfu, const dtype * __restrict__  pmfd,
  const dtype * __restrict__  pa, const dtype * __restrict__  pclv, const dtype * __restrict__  psupsat,
  const dtype * __restrict__  plcrit_aer, const dtype * __restrict__  picrit_aer,
  const dtype * __restrict__  pre_ice, const dtype * __restrict__  pccn, const dtype * __restrict__  pnice,
  dtype * __restrict__  pcovptot, dtype * __restrict__  prainfrac_toprfz,
  dtype * __restrict__  pfsqlf, dtype * __restrict__  pfsqif, dtype * __restrict__  pfcqnng,
  dtype * __restrict__  pfcqlng, dtype * __restrict__  pfsqrf, dtype * __restrict__  pfsqsf,
  dtype * __restrict__  pfcqrng, dtype * __restrict__  pfcqsng,
  dtype * __restrict__  pfsqltur, dtype * __restrict__  pfsqitur,
  dtype * __restrict__  pfplsl, dtype * __restrict__  pfplsn, dtype * __restrict__  pfhpsl,
  dtype * __restrict__  pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  dtype rg, dtype rd, dtype rcpd, dtype retv, dtype rlvtt, dtype rlstt, dtype rlmlt, dtype rtt,
  dtype rv, dtype r2es, dtype r3les, dtype r3ies, dtype r4les, dtype r4ies, dtype r5les,
  dtype r5ies, dtype r5alvcp, dtype r5alscp, dtype ralvdcp, dtype ralsdcp, dtype ralfdcp,
  dtype rtwat, dtype rtice, dtype rticecu, dtype rtwat_rtice_r, dtype rtwat_rticecu_r,
  dtype rkoop1, dtype rkoop2);

