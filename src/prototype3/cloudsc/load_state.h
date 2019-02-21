#include <stdio.h>
#include <stdlib.h>

void query_state(int *klon, int *klev);

void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma, 
		double* ptsphy, double* plcrit_aer, double* picrit_aer,
		double* pre_ice, double* pccn, double* pnice, double* pt, double* pq, 
		double* tend_cml_t, double* tend_cml_q, double* tend_cml_a, double* tend_cml_cld,
		double* tend_tmp_t, double* tend_tmp_q, double* tend_tmp_a, double* tend_tmp_cld,
		double* pvfa, double* pvfl, double* pvfi, double* pdyna, double* pdynl, double* pdyni, 
		double* phrsw, double* phrlw, double* pvervel, double* pap, double* paph, double* plsm,
		int* ldcum, int* ktype, double* plu, double* plude, double* psnde, double* pmfu,
		double* pmfd, double* pa, double* pclv, double* psupsat);

