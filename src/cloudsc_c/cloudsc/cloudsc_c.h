/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef CLOUDSC_C_H
#define CLOUDSC_C_H

#include "yomcst_c.h"
#include "yoethf_c.h"
#include "yoecldp_c.h"

int cloudsc_c(int kidia, int kfdia, int klon, int klev, double ptsphy, double * restrict v_pt, double * restrict v_pq,
	      double * restrict v_tendency_cml_t, double * restrict v_tendency_cml_q, double * restrict v_tendency_cml_a, double * restrict v_tendency_cml_cld,
	      double * restrict v_tendency_tmp_t, double * restrict v_tendency_tmp_q, double * restrict v_tendency_tmp_a, double * restrict v_tendency_tmp_cld,
	      double * restrict v_tendency_loc_t, double * restrict v_tendency_loc_q, double * restrict v_tendency_loc_a, double * restrict v_tendency_loc_cld,
	      double * restrict v_pvfa, double * restrict v_pvfl, double * restrict v_pvfi, double * restrict v_pdyna, double * restrict v_pdynl, double * restrict v_pdyni,
	      double * restrict v_phrsw, double * restrict v_phrlw, double * restrict v_pvervel, double * restrict v_pap, double * restrict v_paph, double * restrict v_plsm,
	      int * restrict v_ktype, double * restrict v_plu, double * restrict v_plude, double * restrict v_psnde, double * restrict v_pmfu,
	      double * restrict v_pmfd, double * restrict v_pa, double * restrict v_pclv, double * restrict v_psupsat, double * restrict v_plcrit_aer, double * restrict v_picrit_aer,
	      double * restrict v_pre_ice, double * restrict v_pccn, double * restrict v_pnice, double * restrict v_pcovptot, double * restrict v_prainfrac_toprfz, double * restrict v_pfsqlf,
	      double * restrict v_pfsqif, double * restrict v_pfcqnng, double * restrict v_pfcqlng, double * restrict v_pfsqrf, double * restrict v_pfsqsf, double * restrict v_pfcqrng,
	      double * restrict v_pfcqsng, double * restrict v_pfsqltur, double * restrict v_pfsqitur, double * restrict v_pfplsl, double * restrict v_pfplsn, double * restrict v_pfhpsl,
	      double * restrict v_pfhpsn);

#endif
