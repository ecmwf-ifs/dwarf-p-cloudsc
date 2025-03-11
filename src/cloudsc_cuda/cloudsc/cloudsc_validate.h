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
#include "dtype.h"

int cloudsc_validate(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		     dtype *plude, dtype *pcovptot, dtype *prainfrac_toprfz, dtype *pfsqlf, dtype *pfsqif,
		     dtype *pfcqlng, dtype *pfcqnng, dtype *pfsqrf, dtype *pfsqsf, dtype *pfcqrng, dtype *pfcqsng,
		     dtype *pfsqltur, dtype *pfsqitur, dtype *pfplsl, dtype *pfplsn, dtype *pfhpsl, dtype *pfhpsn,
		     dtype *tend_loc_a, dtype *tend_loc_q, dtype *tend_loc_t, dtype *tend_loc_cld);
