/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "load_state.h"
//#include <limits>

int cloudsc_validate(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		     float *plude, float *pcovptot, float *prainfrac_toprfz, float *pfsqlf, float *pfsqif,
		     float *pfcqlng, float *pfcqnng, float *pfsqrf, float *pfsqsf, float *pfcqrng, float *pfcqsng,
		     float *pfsqltur, float *pfsqitur, float *pfplsl, float *pfplsn, float *pfhpsl, float *pfhpsn,
		     float *tend_loc_a, float *tend_loc_q, float *tend_loc_t, float *tend_loc_cld);
