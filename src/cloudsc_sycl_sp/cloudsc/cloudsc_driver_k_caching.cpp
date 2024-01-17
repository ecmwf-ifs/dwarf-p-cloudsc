/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "cloudsc_driver.h"
#include "cloudsc_c_k_caching.kernel"

#include <omp.h>
#include "mycpu.h"
#include <CL/sycl.hpp>

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

void cloudsc_driver(int numthreads, int numcols, int nproma) {

  cl::sycl::default_selector device_select;
  cl::sycl::queue q( device_select );

  printf("Running on %s\n", q.get_device().get_info<cl::sycl::info::device::name>().c_str());
	
  float *tend_tmp_u;
  float *tend_tmp_v;
  float *tend_tmp_t;
  float *tend_tmp_q;
  float *tend_tmp_o3;
  float *tend_tmp_a;
  float *tend_tmp_cld;

  float *tend_loc_u;
  float *tend_loc_v;
  float *tend_loc_t;
  float *tend_loc_q;
  float *tend_loc_o3;
  float *tend_loc_a;
  float *tend_loc_cld;

  float *tend_cml_u;
  float *tend_cml_v;
  float *tend_cml_t;
  float *tend_cml_q;
  float *tend_cml_o3;
  float *tend_cml_a;
  float *tend_cml_cld;

  float ptsphy;            //! Physics timestep

  float *plcrit_aer;
  float *picrit_aer;
  float *pre_ice;
  float *pccn;       //! liquid cloud condensation nuclei
  float *pnice;     //! ice number concentration (cf. CCN)
  float *pt;           //! T at start of callpar
  float *pq;           //! Q at start of callpar
  float *pvfa;       //! CC from VDF scheme
  float *pvfl;       //! Liq from VDF scheme
  float *pvfi;       //! Ice from VDF scheme
  float *pdyna;     //! CC from Dynamics
  float *pdynl;     //! Liq from Dynamics
  float *pdyni;     //! Liq from Dynamics
  float *phrsw;     //! Short-wave heating rate
  float *phrlw;     //! Long-wave heating rate
  float *pvervel; //! Vertical velocity
  float *pap;         //! Pressure on full levels
  float *paph;       //! Pressure on half levels
  float *plsm;           //! Land fraction (0-1)
  int *ldcum;
  int    *ktype;          //! Convection type 0,1,2
  float *plu;         //! Conv. condensate
  float *plude;     //! Conv. detrained water
  float *plude_tmp;
  float *psnde;     //! Conv. detrained snow
  float *pmfu;       //! Conv. mass flux up
  float *pmfd;       //! Conv. mass flux down
  float *pa;           //! Original Cloud fraction (t)

  float *pclv;
  float *psupsat;

  float *pcovptot;   //! Precip fraction
  float *prainfrac_toprfz;
  float *pfsqlf;       //! Flux of liquid
  float *pfsqif;       //! Flux of ice
  float *pfcqlng;     //! -ve corr for liq
  float *pfcqnng;     //! -ve corr for ice
  float *pfsqrf;       //! Flux diagnostics
  float *pfsqsf;       //!    for DDH, generic
  float *pfcqrng;     //! rain
  float *pfcqsng;     //! snow
  float *pfsqltur;   //! liquid flux due to VDF
  float *pfsqitur;   //! ice flux due to VDF
  float *pfplsl;       //! liq+rain sedim flux
  float *pfplsn;       //! ice+snow sedim flux
  float *pfhpsl;       //! Enthalpy flux for liq
  float *pfhpsn;       //! Enthalpy flux for ice

  /* Define or query data dimensions from input file */
  int klon, nlev;
  int nblocks = (numcols / nproma) + min(numcols % nproma, 1);

  float zinfo[4][numthreads];
  const float zhpm = 12482329.0;  // IBM P7 HPM flop count for 100 points at L137

  int nclv;      // number of microphysics variables
  int ncldql;    // liquid cloud water
  int ncldqi;    // ice cloud water
  int ncldqr;    // rain water
  int ncldqs;    // snow
  int ncldqv;    // vapour

  nclv = 5;      // number of microphysics variables
  ncldql = 1;    // liquid cloud water
  ncldqi = 2;    // ice cloud water
  ncldqr = 3;    // rain water
  ncldqs = 4;    // snow
  ncldqv = 5;    // vapour

  struct TECLDP *yrecldp = (struct TECLDP*)malloc(sizeof(struct TECLDP));

  query_state(&klon, &nlev);

  tend_loc_t   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_loc_q   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_loc_a   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_loc_cld = (float*) malloc( sizeof(float) * nblocks*nlev*nproma*nclv );
  tend_cml_t   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_cml_q   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_cml_a   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_cml_cld = (float*) malloc( sizeof(float) * nblocks*nlev*nproma*nclv );
  tend_tmp_t   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_tmp_q   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_tmp_a   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  tend_tmp_cld = (float*) malloc( sizeof(float) * nblocks*nlev*nproma*nclv );

  plcrit_aer = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  picrit_aer = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pre_ice    = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pccn       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pnice      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pt         = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pq         = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pvfa       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pvfl       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pvfi       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pdyna      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pdynl      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pdyni      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  phrsw      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  phrlw      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pvervel    = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pap        = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  paph       = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  plsm       = (float*) malloc( sizeof(float) * nblocks*nproma );
  ldcum      = (int*) malloc( sizeof(int) * nblocks*nproma );
  ktype      = (int*) malloc( sizeof(int) * nblocks*nproma );
  plu        = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  plude      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  psnde      = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pmfu       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pmfd       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pa         = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pclv       = (float*) malloc( sizeof(float) * nblocks*nlev*nproma*nclv );
  psupsat    = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  pcovptot   = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  prainfrac_toprfz = (float*) malloc( sizeof(float) * nblocks*nproma );
  pfsqlf    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfsqif    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfcqnng   = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfcqlng   = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfsqrf    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfsqsf    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfcqrng   = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfcqsng   = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfsqltur  = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfsqitur  = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfplsl    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfplsn    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfhpsl    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  pfhpsn    = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );


  float rg;
  float rd;
  float rcpd;
  float retv;
  float rlvtt;
  float rlstt;
  float rlmlt;
  float rtt;
  float rv;
  float r2es;
  float r3les;
  float r3ies;
  float r4les;
  float r4ies;
  float r5les;
  float r5ies;
  float r5alvcp;
  float r5alscp;
  float ralvdcp;
  float ralsdcp;
  float ralfdcp;
  float rtwat;
  float rtice;
  float rticecu;
  float rtwat_rtice_r;
  float rtwat_rticecu_r;
  float rkoop1;
  float rkoop2;

  // device declarations
  float *d_plcrit_aer;
  float *d_picrit_aer;
  float *d_pre_ice;
  float *d_pccn;
  float *d_pnice;
  float *d_pt;
  float *d_pq;
  float *d_tend_loc_t;
  float *d_tend_loc_q;
  float *d_tend_loc_a;
  float *d_tend_loc_cld;
  float *d_tend_tmp_t;
  float *d_tend_tmp_q;
  float *d_tend_tmp_a;
  float *d_tend_tmp_cld;
  float *d_pvfa;
  float *d_pvfl;
  float *d_pvfi;
  float *d_pdyna;
  float *d_pdynl;
  float *d_pdyni;
  float *d_phrsw;
  float *d_phrlw;
  float *d_pvervel;
  float *d_pap;
  float *d_paph;
  float *d_plsm;
  int *d_ktype;
  float *d_plu;
  float *d_plude;
  float *d_psnde;
  float *d_pmfu;
  float *d_pmfd;
  float *d_pa;
  float *d_pclv;
  float *d_psupsat;
  struct TECLDP *d_yrecldp;
  float *d_pcovptot;
  float *d_prainfrac_toprfz;
  float *d_pfsqlf;
  float *d_pfsqif;
  float *d_pfcqnng;
  float *d_pfcqlng;
  float *d_pfsqrf;
  float *d_pfsqsf;
  float *d_pfcqrng;
  float *d_pfcqsng;
  float *d_pfsqltur;
  float *d_pfsqitur;
  float *d_pfplsl;
  float *d_pfplsn;
  float *d_pfhpsl;
  float *d_pfhpsn;
  // end device declarations


  d_plcrit_aer = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_picrit_aer = cl::sycl::malloc_device<float>( nblocks*nlev*nproma, q);
  d_pre_ice = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pccn = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pnice = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pt = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pq = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_loc_t = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_loc_q = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_loc_a = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_loc_cld = cl::sycl::malloc_device<float>(nblocks*nlev*nproma*nclv, q);
  d_tend_tmp_t = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_tmp_q = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_tmp_a = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_tend_tmp_cld = cl::sycl::malloc_device<float>(nblocks*nlev*nproma*nclv, q);
  d_pvfa = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pvfl = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pvfi = cl::sycl::malloc_device<float>(nblocks*nlev*nproma,q );
  d_pdyna = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pdynl = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pdyni = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_phrsw = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_phrlw = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pvervel = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pap = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_paph = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_plsm = cl::sycl::malloc_device<float>(nblocks*nproma, q);
  d_ktype = cl::sycl::malloc_device<int>(nblocks*nproma, q);
  d_plu = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_plude = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_psnde = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pmfu = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pmfd = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pa = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_pclv = cl::sycl::malloc_device<float>(nblocks*nlev*nproma*nclv, q);
  d_psupsat = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_yrecldp = (TECLDP*) cl::sycl::malloc_device( sizeof(TECLDP), q);
  d_pcovptot = cl::sycl::malloc_device<float>(nblocks*nlev*nproma, q);
  d_prainfrac_toprfz = cl::sycl::malloc_device<float>(nblocks*nproma, q);
  d_pfsqlf = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfsqif = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfcqnng = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfcqlng = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfsqrf = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfsqsf = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfcqrng = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfcqsng = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfsqltur = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfsqitur = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfplsl = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfplsn = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfhpsl = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);
  d_pfhpsn = cl::sycl::malloc_device<float>(nblocks*(nlev+1)*nproma, q);

  int n_runs = 2;
  for (int i_runs=0; i_runs<n_runs; i_runs++) {
  load_state(klon, nlev, nclv, numcols, nproma, &ptsphy, plcrit_aer, picrit_aer,
	     pre_ice, pccn, pnice, pt, pq,
	     tend_cml_t, tend_cml_q, tend_cml_a, tend_cml_cld,
	     tend_tmp_t, tend_tmp_q, tend_tmp_a, tend_tmp_cld,
	     pvfa, pvfl, pvfi, pdyna, pdynl, pdyni,
	     phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu,
	     plude, psnde, pmfu, pmfd, pa, pclv, psupsat, yrecldp,
	     &rg, &rd, &rcpd, &retv, &rlvtt, &rlstt,
             &rlmlt, &rtt, &rv, &r2es, &r3les, &r3ies,
             &r4les, &r4ies, &r5les, &r5ies, &r5alvcp, &r5alscp,
             &ralvdcp, &ralsdcp, &ralfdcp, &rtwat,
             &rtice, &rticecu, &rtwat_rtice_r, &rtwat_rticecu_r,
             &rkoop1, &rkoop2 );

  // host to device
  q.memcpy(d_plcrit_aer, plcrit_aer, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_picrit_aer, picrit_aer, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pre_ice, pre_ice, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pccn, pccn, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pnice, pnice, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pt, pt, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pq, pq, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_loc_t, tend_loc_t, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_loc_q, tend_loc_q, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_loc_a, tend_loc_a, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_loc_cld, tend_loc_cld, sizeof(float) * nblocks*nlev*nproma*nclv);
  q.memcpy(d_tend_tmp_t, tend_tmp_t, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_tmp_q, tend_tmp_q, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_tmp_a, tend_tmp_a, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_tend_tmp_cld, tend_tmp_cld, sizeof(float) * nblocks*nlev*nproma*nclv);
  q.memcpy(d_pvfa, pvfa, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pvfl, pvfl, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pvfi, pvfi, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pdyna, pdyna, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pdynl, pdynl, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pdyni, pdyni, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_phrsw, phrsw, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_phrlw, phrlw, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pvervel, pvervel, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pap, pap, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_paph, paph, sizeof(float) * nblocks*(nlev+1)*nproma);
  q.memcpy(d_plsm, plsm, sizeof(float) * nblocks*nproma);
  q.memcpy(d_ktype, ktype, sizeof(int) * nblocks*nproma);
  q.memcpy(d_plu, plu, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_plude, plude, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_psnde, psnde, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pmfu, pmfu, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pmfd, pmfd, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pa, pa, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_pclv, pclv, sizeof(float) * nblocks*nlev*nproma*nclv);
  q.memcpy(d_psupsat, psupsat, sizeof(float) * nblocks*nlev*nproma);
  q.memcpy(d_yrecldp, yrecldp, sizeof(TECLDP));
  q.wait();
  // end host to device

  float t1 = omp_get_wtime();

    int b, bsize, icalls=0, igpc=numcols;
    int coreid = mycpu();
    int tid = omp_get_thread_num();
    float start = omp_get_wtime();

    int jkglo = 0;
    int ibl = (jkglo - 1) / nproma + 1;
    int icend = min(nproma, numcols - jkglo + 1);


    cl::sycl::range<1> global(numcols);
    cl::sycl::range<1> local(nproma);

    q.submit([&](cl::sycl::handler &h) {
        /*cl::sycl::stream out_stream(16384, 512, h);*/
        h.parallel_for( cl::sycl::nd_range<1>( global, local), [=] (cl::sycl::nd_item<1> item_ct1) /* TODO: insert: [intel::reqd_sub_group_size(16)]*/ {			 

    cloudsc_c(1, icend, nproma, ptsphy, d_pt, d_pq,
    		d_tend_tmp_t, d_tend_tmp_q, d_tend_tmp_a, d_tend_tmp_cld,
    		d_tend_loc_t, d_tend_loc_q, d_tend_loc_a, d_tend_loc_cld,
    		d_pvfa, d_pvfl, d_pvfi,
    		d_pdyna, d_pdynl, d_pdyni,
    		d_phrsw, d_phrlw, d_pvervel,
    		d_pap, d_paph, d_plsm, d_ktype,
    		d_plu, d_plude, d_psnde, d_pmfu, d_pmfd,
    		d_pa, d_pclv, d_psupsat,
    		d_plcrit_aer, d_picrit_aer, d_pre_ice, d_pccn, d_pnice,
    		d_pcovptot, d_prainfrac_toprfz, d_pfsqlf,
    		d_pfsqif, d_pfcqnng, d_pfcqlng,
    		d_pfsqrf, d_pfsqsf, d_pfcqrng,
    		d_pfcqsng, d_pfsqltur, d_pfsqitur,
    		d_pfplsl, d_pfplsn, d_pfhpsl, d_pfhpsn, d_yrecldp,
    		nblocks, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt,
                rv, r2es, r3les, r3ies, r4les, r4ies, r5les,
                r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp,
                rtwat, rtice, rticecu, rtwat_rtice_r, rtwat_rticecu_r,
                rkoop1, rkoop2, /*out_stream,*/ item_ct1);

	});
    });

    q.wait();

    float end = omp_get_wtime();

    // device to host
    q.memcpy(tend_loc_t, d_tend_loc_t, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(tend_loc_q, d_tend_loc_q, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(tend_loc_a, d_tend_loc_a, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(tend_loc_cld, d_tend_loc_cld, sizeof(float) * nblocks*nlev*nproma*nclv);
    q.memcpy(phrlw, d_phrlw, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(plude, d_plude, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(yrecldp, d_yrecldp, sizeof(TECLDP));
    q.memcpy(pcovptot, d_pcovptot, sizeof(float) * nblocks*nlev*nproma);
    q.memcpy(prainfrac_toprfz, d_prainfrac_toprfz, sizeof(float) * nblocks*nproma);
    q.memcpy(pfsqlf, d_pfsqlf, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfsqif, d_pfsqif, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfcqnng, d_pfcqnng, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfcqlng, d_pfcqlng, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfsqrf, d_pfsqrf, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfsqsf, d_pfsqsf, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfcqrng, d_pfcqrng, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfcqsng, d_pfcqsng, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfsqltur, d_pfsqltur, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfsqitur, d_pfsqitur, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfplsl, d_pfplsl, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfplsn, d_pfplsn, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfhpsl, d_pfhpsl, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.memcpy(pfhpsn, d_pfhpsn, sizeof(float) * nblocks*(nlev+1)*nproma);
    q.wait();
    // end device to host

    /* int msec = diff * 1000 / CLOCKS_PER_SEC; */
    zinfo[0][tid] = end - start;
    zinfo[1][tid] = (float) coreid;
    zinfo[2][tid] = (float) icalls;
    zinfo[3][tid] = (float) igpc;

  float t2 = omp_get_wtime();

  printf("     NUMOMP=%d, NGPTOT=%d, NPROMA=%d, NGPBLKS=%d\n", numthreads, numcols, nproma, nblocks);
  printf(" %+10s%+10s%+10s%+10s%+10s %+4s : %+10s%+10s%+10s\n",
    "NUMOMP", "NGPTOT", "#GP-cols", "#BLKS", "NPROMA", "tid#", "Time(msec)", "MFlops/s", "col/s");
  float zfrac, zmflops, zthrput;
  for (int t = 0; t < numthreads; t++) {
    const float tloc = zinfo[0][t];
    const int coreid = (int) zinfo[1][t];
    const int icalls = (int) zinfo[2][t];
    const int igpc = (int) zinfo[3][t];
    zfrac = (float)igpc / (float)numcols;
    if (tloc > 0.0) {
      zmflops = 1.0e-06 * zfrac * zhpm * ((float)numcols / 100.) / tloc;
      zthrput = (float)numcols/tloc;
    } else {
      zmflops = 0.;
      zthrput = 0.0;
    }
    printf(" %10d%10d%10d%10d%10d %4d : %10d%10d%10d @ core#\n",
           numthreads, numcols, igpc, icalls, nproma, t, (int)(tloc * 1000.), (int)zmflops, (int)zthrput);
  }
  float tdiff = t2 - t1;
  zfrac = 1.0;
  if (tdiff > 0.0) {
    zmflops = 1.0e-06 * zfrac * zhpm * ((float)numcols / 100.) / tdiff;
    zthrput = (float)numcols/tdiff;
  } else {
    zmflops = 0.0;
    zthrput = 0.0;
  }
  printf(" %10d%10d%10d%10d%10d %4d: %10d%10d%10d TOTAL\n",
         numthreads, numcols, numcols, nblocks, nproma, -1, (int)(tdiff * 1000.), (int)zmflops, (int)zthrput);
 
  } // n_runs
  cloudsc_validate(klon, nlev, nclv, numcols, nproma,
		   plude, pcovptot, prainfrac_toprfz, pfsqlf, pfsqif,
		   pfcqlng, pfcqnng, pfsqrf, pfsqsf, pfcqrng, pfcqsng,
		   pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn,
		   tend_loc_a, tend_loc_q, tend_loc_t, tend_loc_cld);

  free(plcrit_aer);
  free(picrit_aer);
  free(pre_ice);
  free(pccn);
  free(pnice);
  free(pt);
  free(pq);
  free(pvfa);
  free(pvfl);
  free(pvfi);
  free(pdyna);
  free(pdynl);
  free(pdyni);
  free(phrsw);
  free(phrlw);
  free(pvervel);
  free(pap);
  free(paph);
  free(plsm);
  free(ktype);
  free(plu);
  free(plude);
  free(psnde);
  free(pmfu);
  free(pmfd);
  free(pa);
  free(pclv);
  free(psupsat);
  free(pcovptot);
  free(tend_loc_t);
  free(tend_loc_q);
  free(tend_loc_a);
  free(tend_loc_cld);
  free(tend_tmp_t);
  free(tend_tmp_q);
  free(tend_tmp_a);
  free(tend_tmp_cld);
  free(tend_cml_t);
  free(tend_cml_q);
  free(tend_cml_a);
  free(tend_cml_cld);
  free(prainfrac_toprfz);
  free(pfsqlf);
  free(pfsqif);
  free(pfcqnng);
  free(pfcqlng);
  free(pfsqrf);
  free(pfsqsf);
  free(pfcqrng);
  free(pfcqsng);
  free(pfsqltur);
  free(pfsqitur);
  free(pfplsl);
  free(pfplsn);
  free(pfhpsl);
  free(pfhpsn);
  free(yrecldp);

  // free device
  cl::sycl::free(d_plcrit_aer, q);
  cl::sycl::free(d_picrit_aer, q);
  cl::sycl::free(d_pre_ice, q);
  cl::sycl::free(d_pccn, q);
  cl::sycl::free(d_pnice, q);
  cl::sycl::free(d_pt, q);
  cl::sycl::free(d_pq, q);
  cl::sycl::free(d_tend_loc_t, q);
  cl::sycl::free(d_tend_loc_q, q);
  cl::sycl::free(d_tend_loc_a, q);
  cl::sycl::free(d_tend_loc_cld, q);
  cl::sycl::free(d_tend_tmp_t, q);
  cl::sycl::free(d_tend_tmp_q, q);
  cl::sycl::free(d_tend_tmp_a, q);
  cl::sycl::free(d_tend_tmp_cld, q);
  cl::sycl::free(d_pvfa, q);
  cl::sycl::free(d_pvfl, q);
  cl::sycl::free(d_pvfi, q);
  cl::sycl::free(d_pdyna, q);
  cl::sycl::free(d_pdynl, q);
  cl::sycl::free(d_pdyni, q);
  cl::sycl::free(d_phrsw, q);
  cl::sycl::free(d_phrlw, q);
  cl::sycl::free(d_pvervel, q);
  cl::sycl::free(d_pap, q);
  cl::sycl::free(d_paph, q);
  cl::sycl::free(d_plsm, q);
  cl::sycl::free(d_ktype, q);
  cl::sycl::free(d_plu, q);
  cl::sycl::free(d_plude, q);
  cl::sycl::free(d_psnde, q);
  cl::sycl::free(d_pmfu, q);
  cl::sycl::free(d_pmfd, q);
  cl::sycl::free(d_pa, q);
  cl::sycl::free(d_pclv, q);
  cl::sycl::free(d_psupsat, q);
  cl::sycl::free(d_yrecldp, q);
  cl::sycl::free(d_pcovptot, q);
  cl::sycl::free(d_prainfrac_toprfz, q);
  cl::sycl::free(d_pfsqlf, q);
  cl::sycl::free(d_pfsqif, q);
  cl::sycl::free(d_pfcqnng, q);
  cl::sycl::free(d_pfcqlng, q);
  cl::sycl::free(d_pfsqrf, q);
  cl::sycl::free(d_pfsqsf, q);
  cl::sycl::free(d_pfcqrng, q);
  cl::sycl::free(d_pfcqsng, q);
  cl::sycl::free(d_pfsqltur, q);
  cl::sycl::free(d_pfsqitur, q);
  cl::sycl::free(d_pfplsl, q);
  cl::sycl::free(d_pfplsn, q);
  cl::sycl::free(d_pfhpsl, q);
  cl::sycl::free(d_pfhpsn, q);
  // end free device
}

