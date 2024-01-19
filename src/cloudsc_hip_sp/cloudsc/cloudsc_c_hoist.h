#include "hip/hip_runtime.h"
/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "hip/hip_runtime.h"
#include "yoecldp_c.h"
#include <stdio.h>

__global__ void __launch_bounds__(128, 1) cloudsc_c(int kidia, int kfdia, int klon, float ptsphy,
  const float * __restrict__  pt,
  const float * __restrict__  pq, const float * __restrict__  tendency_tmp_t,
  const float * __restrict__  tendency_tmp_q, const float * __restrict__  tendency_tmp_a,
  const float * __restrict__  tendency_tmp_cld, float * __restrict__  tendency_loc_t,
  float * __restrict__  tendency_loc_q, float * __restrict__  tendency_loc_a,
  float * __restrict__  tendency_loc_cld, const float * __restrict__  pvfa,
  const float * __restrict__  pvfl, const float * __restrict__  pvfi, const float * __restrict__  pdyna,
  const float * __restrict__  pdynl, const float * __restrict__  pdyni, const float * __restrict__  phrsw,
  float * __restrict__  phrlw, const float * __restrict__  pvervel, const float * __restrict__  pap,
  const float * __restrict__  paph, const float * __restrict__  plsm,
  const int *  ktype, const float * __restrict__  plu, float * __restrict__  plude,
  const float * __restrict__  psnde, const float * __restrict__  pmfu, const float * __restrict__  pmfd,
  const float * __restrict__  pa, const float * __restrict__  pclv, const float * __restrict__  psupsat,
  const float * __restrict__  plcrit_aer, const float * __restrict__  picrit_aer,
  const float * __restrict__  pre_ice, const float * __restrict__  pccn, const float * __restrict__  pnice,
  float * __restrict__  pcovptot, float * __restrict__  prainfrac_toprfz,
  float * __restrict__  pfsqlf, float * __restrict__  pfsqif, float * __restrict__  pfcqnng,
  float * __restrict__  pfcqlng, float * __restrict__  pfsqrf, float * __restrict__  pfsqsf,
  float * __restrict__  pfcqrng, float * __restrict__  pfcqsng,
  float * __restrict__  pfsqltur, float * __restrict__  pfsqitur,
  float * __restrict__  pfplsl, float * __restrict__  pfplsn, float * __restrict__  pfhpsl,
  float * __restrict__  pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  float rg, float rd, float rcpd, float retv, float rlvtt, float rlstt, float rlmlt, float rtt,
  float rv, float r2es, float r3les, float r3ies, float r4les, float r4ies, float r5les,
  float r5ies, float r5alvcp, float r5alscp, float ralvdcp, float ralsdcp, float ralfdcp,
  float rtwat, float rtice, float rticecu, float rtwat_rtice_r, float rtwat_rticecu_r,
  float rkoop1, float rkoop2,
  float * __restrict__ zfoealfa, float * __restrict__ ztp1, float * __restrict__ zli,
  float * __restrict__ za, float * __restrict__ zaorig, float * __restrict__ zliqfrac,
  float * __restrict__ zicefrac, float * __restrict__ zqx, float * __restrict__ zqx0,
  float * __restrict__ zpfplsx, float * __restrict__ zlneg, float * __restrict__ zqxn2d,
  float * __restrict__ zqsmix, float * __restrict__ zqsliq, float * __restrict__ zqsice,
  float * __restrict__ zfoeewmt, float * __restrict__ zfoeew, float * __restrict__ zfoeeliqt);

