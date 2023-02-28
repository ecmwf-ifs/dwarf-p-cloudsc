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

__global__ void cloudsc_c(int kidia, int kfdia, int klon, int klev, double ptsphy,
  const double * __restrict__  pt,
  const double * __restrict__  pq, const double * __restrict__  tendency_tmp_t,
  const double * __restrict__  tendency_tmp_q, const double * __restrict__  tendency_tmp_a,
  const double * __restrict__  tendency_tmp_cld, double * __restrict__  tendency_loc_t,
  double * __restrict__  tendency_loc_q, double * __restrict__  tendency_loc_a,
  double * __restrict__  tendency_loc_cld, const double * __restrict__  pvfa,
  const double * __restrict__  pvfl, const double * __restrict__  pvfi, const double * __restrict__  pdyna,
  const double * __restrict__  pdynl, const double * __restrict__  pdyni, const double * __restrict__  phrsw,
  double * __restrict__  phrlw, const double * __restrict__  pvervel, const double * __restrict__  pap,
  const double * __restrict__  paph, const double * __restrict__  plsm,
  const int *  ktype, const double * __restrict__  plu, double * __restrict__  plude,
  const double * __restrict__  psnde, const double * __restrict__  pmfu, const double * __restrict__  pmfd,
  const double * __restrict__  pa, const double * __restrict__  pclv, const double * __restrict__  psupsat,
  const double * __restrict__  plcrit_aer, const double * __restrict__  picrit_aer,
  const double * __restrict__  pre_ice, const double * __restrict__  pccn, const double * __restrict__  pnice,
  double * __restrict__  pcovptot, double * __restrict__  prainfrac_toprfz,
  double * __restrict__  pfsqlf, double * __restrict__  pfsqif, double * __restrict__  pfcqnng,
  double * __restrict__  pfcqlng, double * __restrict__  pfsqrf, double * __restrict__  pfsqsf,
  double * __restrict__  pfcqrng, double * __restrict__  pfcqsng,
  double * __restrict__  pfsqltur, double * __restrict__  pfsqitur,
  double * __restrict__  pfplsl, double * __restrict__  pfplsn, double * __restrict__  pfhpsl,
  double * __restrict__  pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  double rg, double rd, double rcpd, double retv, double rlvtt, double rlstt, double rlmlt, double rtt,
  double rv, double r2es, double r3les, double r3ies, double r4les, double r4ies, double r5les,
  double r5ies, double r5alvcp, double r5alscp, double ralvdcp, double ralsdcp, double ralfdcp,
  double rtwat, double rtice, double rticecu, double rtwat_rtice_r, double rtwat_rticecu_r,
  double rkoop1, double rkoop2,
  double * __restrict__ zfoealfa, double * __restrict__ ztp1, double * __restrict__ zli,
  double * __restrict__ za, double * __restrict__ zaorig, double * __restrict__ zliqfrac,
  double * __restrict__ zicefrac, double * __restrict__ zqx, double * __restrict__ zqx0,
  double * __restrict__ zpfplsx, double * __restrict__ zlneg, double * __restrict__ zqxn2d,
  double * __restrict__ zqsmix, double * __restrict__ zqsliq, double * __restrict__ zqsice,
  double * __restrict__ zfoeewmt, double * __restrict__ zfoeew, double * __restrict__ zfoeeliqt);

