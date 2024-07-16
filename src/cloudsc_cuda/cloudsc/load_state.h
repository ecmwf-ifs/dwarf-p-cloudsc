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
#include "dtype.h"

struct TECLDP ;

void query_state(int *klon, int *klev);

void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		dtype* ptsphy, dtype* plcrit_aer, dtype* picrit_aer,
		dtype* pre_ice, dtype* pccn, dtype* pnice, dtype* pt, dtype* pq, 
		dtype* tend_cml_t, dtype* tend_cml_q, dtype* tend_cml_a, dtype* tend_cml_cld,
		dtype* tend_tmp_t, dtype* tend_tmp_q, dtype* tend_tmp_a, dtype* tend_tmp_cld,
		dtype* pvfa, dtype* pvfl, dtype* pvfi, dtype* pdyna, dtype* pdynl, dtype* pdyni, 
		dtype* phrsw, dtype* phrlw, dtype* pvervel, dtype* pap, dtype* paph, dtype* plsm,
		int* ktype, dtype* plu, dtype* plude, dtype* psnde, dtype* pmfu,
		dtype* pmfd, dtype* pa, dtype* pclv, dtype* psupsat, struct TECLDP* yrecldp,
		dtype* rg, dtype* rd, dtype* rcpd, dtype* retv, dtype* rlvtt, dtype* rlstt,
                dtype* rlmlt, dtype* rtt, dtype* rv, dtype* r2es, dtype* r3les, dtype* r3ies,
                dtype* r4les, dtype* r4ies, dtype* r5les, dtype* r5ies, dtype* r5alvcp, dtype* r5alscp,
                dtype* ralvdcp, dtype* ralsdcp, dtype* ralfdcp, dtype* rtwat,
                dtype* rtice, dtype* rticecu, dtype* rtwat_rtice_r, dtype *rtwat_rticecu_r,
                dtype* rkoop1, dtype* rkoop2);


void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    dtype *plude, dtype *pcovptot, dtype *prainfrac_toprfz, dtype *pfsqlf, dtype *pfsqif,
		    dtype *pfcqlng, dtype *pfcqnng, dtype *pfsqrf, dtype *pfsqsf, dtype *pfcqrng, dtype *pfcqsng,
		    dtype *pfsqltur, dtype *pfsqitur, dtype *pfplsl, dtype *pfplsn, dtype *pfhpsl, dtype *pfhpsn,
		    dtype *tend_loc_a, dtype *tend_loc_q, dtype *tend_loc_t, dtype *tend_loc_cld);
