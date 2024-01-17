/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include "yoecldp_c.h"

struct TECLDP ;

void query_state(int *klon, int *klev);

void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		float* ptsphy, float* plcrit_aer, float* picrit_aer,
		float* pre_ice, float* pccn, float* pnice, float* pt, float* pq, 
		float* tend_cml_t, float* tend_cml_q, float* tend_cml_a, float* tend_cml_cld,
		float* tend_tmp_t, float* tend_tmp_q, float* tend_tmp_a, float* tend_tmp_cld,
		float* pvfa, float* pvfl, float* pvfi, float* pdyna, float* pdynl, float* pdyni, 
		float* phrsw, float* phrlw, float* pvervel, float* pap, float* paph, float* plsm,
		int* ktype, float* plu, float* plude, float* psnde, float* pmfu,
		float* pmfd, float* pa, float* pclv, float* psupsat, struct TECLDP* yrecldp,
		float* rg, float* rd, float* rcpd, float* retv, float* rlvtt, float* rlstt,
                float* rlmlt, float* rtt, float* rv, float* r2es, float* r3les, float* r3ies,
                float* r4les, float* r4ies, float* r5les, float* r5ies, float* r5alvcp, float* r5alscp,
                float* ralvdcp, float* ralsdcp, float* ralfdcp, float* rtwat,
                float* rtice, float* rticecu, float* rtwat_rtice_r, float *rtwat_rticecu_r,
                float* rkoop1, float* rkoop2);


void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    float *plude, float *pcovptot, float *prainfrac_toprfz, float *pfsqlf, float *pfsqif,
		    float *pfcqlng, float *pfcqnng, float *pfsqrf, float *pfsqsf, float *pfcqrng, float *pfcqsng,
		    float *pfsqltur, float *pfsqitur, float *pfplsl, float *pfplsn, float *pfhpsl, float *pfhpsn,
		    float *tend_loc_a, float *tend_loc_q, float *tend_loc_t, float *tend_loc_cld);
