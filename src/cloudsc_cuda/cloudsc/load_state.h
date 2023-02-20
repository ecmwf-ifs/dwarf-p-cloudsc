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
		double* ptsphy, double* plcrit_aer, double* picrit_aer,
		double* pre_ice, double* pccn, double* pnice, double* pt, double* pq, 
		double* tend_cml_t, double* tend_cml_q, double* tend_cml_a, double* tend_cml_cld,
		double* tend_tmp_t, double* tend_tmp_q, double* tend_tmp_a, double* tend_tmp_cld,
		double* pvfa, double* pvfl, double* pvfi, double* pdyna, double* pdynl, double* pdyni, 
		double* phrsw, double* phrlw, double* pvervel, double* pap, double* paph, double* plsm,
		int* ktype, double* plu, double* plude, double* psnde, double* pmfu,
		double* pmfd, double* pa, double* pclv, double* psupsat, struct TECLDP* yrecldp,
		double* rg, double* rd, double* rcpd, double* retv, double* rlvtt, double* rlstt,
                double* rlmlt, double* rtt, double* rv, double* r2es, double* r3les, double* r3ies,
                double* r4les, double* r4ies, double* r5les, double* r5ies, double* r5alvcp, double* r5alscp,
                double* ralvdcp, double* ralsdcp, double* ralfdcp, double* rtwat,
                double* rtice, double* rticecu, double* rtwat_rtice_r, double *rtwat_rticecu_r,
                double* rkoop1, double* rkoop2);


void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    double *plude, double *pcovptot, double *prainfrac_toprfz, double *pfsqlf, double *pfsqif,
		    double *pfcqlng, double *pfcqnng, double *pfsqrf, double *pfsqsf, double *pfcqrng, double *pfcqsng,
		    double *pfsqltur, double *pfsqitur, double *pfplsl, double *pfplsn, double *pfhpsl, double *pfhpsn,
		    double *tend_loc_a, double *tend_loc_q, double *tend_loc_t, double *tend_loc_cld);
