#include <stdio.h>
#include <stdlib.h>

#include "yomcst_c.h"
#include "yoethf_c.h"
#include "yoecldp_c.h"


void query_state(int *klon, int *klev);

void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		double* ptsphy, double* plcrit_aer, double* picrit_aer,
		double* pre_ice, double* pccn, double* pnice, double* pt, double* pq, 
		double* tend_cml_t, double* tend_cml_q, double* tend_cml_a, double* tend_cml_cld,
		double* tend_tmp_t, double* tend_tmp_q, double* tend_tmp_a, double* tend_tmp_cld,
		double* pvfa, double* pvfl, double* pvfi, double* pdyna, double* pdynl, double* pdyni, 
		double* phrsw, double* phrlw, double* pvervel, double* pap, double* paph, double* plsm,
		int* ktype, double* plu, double* plude, double* psnde, double* pmfu,
		double* pmfd, double* pa, double* pclv, double* psupsat);


void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    double *plude, double *pcovptot, double *prainfrac_toprfz, double *pfsqlf, double *pfsqif,
		    double *pfcqlng, double *pfcqnng, double *pfsqrf, double *pfsqsf, double *pfcqrng, double *pfcqsng,
		    double *pfsqltur, double *pfsqitur, double *pfplsl, double *pfplsn, double *pfhpsl, double *pfhpsn,
		    double *tend_loc_a, double *tend_loc_q, double *tend_loc_t, double *tend_loc_cld);
