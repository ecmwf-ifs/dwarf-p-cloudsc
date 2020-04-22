#include "load_state.h"

int cloudsc_validate(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		     double *plude, double *pcovptot, double *prainfrac_toprfz, double *pfsqlf, double *pfsqif,
		     double *pfcqlng, double *pfcqnng, double *pfsqrf, double *pfsqsf, double *pfcqrng, double *pfcqsng,
		     double *pfsqltur, double *pfsqitur, double *pfplsl, double *pfplsn, double *pfhpsl, double *pfhpsn,
		     double *tend_loc_a, double *tend_loc_q, double *tend_loc_t, double *tend_loc_cld);
