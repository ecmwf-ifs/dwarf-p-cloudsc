/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "cloudsc_driver_hoist.h"

#include <omp.h>
#include "mycpu.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

void cloudsc_driver(int numthreads, int numcols, int nproma) {

  dtype *tend_tmp_u;
  dtype *tend_tmp_v;
  dtype *tend_tmp_t;
  dtype *tend_tmp_q;
  dtype *tend_tmp_o3;
  dtype *tend_tmp_a;
  dtype *tend_tmp_cld;

  dtype *tend_loc_u;
  dtype *tend_loc_v;
  dtype *tend_loc_t;
  dtype *tend_loc_q;
  dtype *tend_loc_o3;
  dtype *tend_loc_a;
  dtype *tend_loc_cld;

  dtype *tend_cml_u;
  dtype *tend_cml_v;
  dtype *tend_cml_t;
  dtype *tend_cml_q;
  dtype *tend_cml_o3;
  dtype *tend_cml_a;
  dtype *tend_cml_cld;

  dtype ptsphy;            //! Physics timestep

  dtype *plcrit_aer;
  dtype *picrit_aer;
  dtype *pre_ice;
  dtype *pccn;       //! liquid cloud condensation nuclei
  dtype *pnice;     //! ice number concentration (cf. CCN)
  dtype *pt;           //! T at start of callpar
  dtype *pq;           //! Q at start of callpar
  dtype *pvfa;       //! CC from VDF scheme
  dtype *pvfl;       //! Liq from VDF scheme
  dtype *pvfi;       //! Ice from VDF scheme
  dtype *pdyna;     //! CC from Dynamics
  dtype *pdynl;     //! Liq from Dynamics
  dtype *pdyni;     //! Liq from Dynamics
  dtype *phrsw;     //! Short-wave heating rate
  dtype *phrlw;     //! Long-wave heating rate
  dtype *pvervel; //! Vertical velocity
  dtype *pap;         //! Pressure on full levels
  dtype *paph;       //! Pressure on half levels
  dtype *plsm;           //! Land fraction (0-1)
  int *ldcum;
  int    *ktype;          //! Convection type 0,1,2
  dtype *plu;         //! Conv. condensate
  dtype *plude;     //! Conv. detrained water
  dtype *plude_tmp;
  dtype *psnde;     //! Conv. detrained snow
  dtype *pmfu;       //! Conv. mass flux up
  dtype *pmfd;       //! Conv. mass flux down
  dtype *pa;           //! Original Cloud fraction (t)

  dtype *pclv;
  dtype *psupsat;

  dtype *pcovptot;   //! Precip fraction
  dtype *prainfrac_toprfz;
  dtype *pfsqlf;       //! Flux of liquid
  dtype *pfsqif;       //! Flux of ice
  dtype *pfcqlng;     //! -ve corr for liq
  dtype *pfcqnng;     //! -ve corr for ice
  dtype *pfsqrf;       //! Flux diagnostics
  dtype *pfsqsf;       //!    for DDH, generic
  dtype *pfcqrng;     //! rain
  dtype *pfcqsng;     //! snow
  dtype *pfsqltur;   //! liquid flux due to VDF
  dtype *pfsqitur;   //! ice flux due to VDF
  dtype *pfplsl;       //! liq+rain sedim flux
  dtype *pfplsn;       //! ice+snow sedim flux
  dtype *pfhpsl;       //! Enthalpy flux for liq
  dtype *pfhpsn;       //! Enthalpy flux for ice

  /* Define or query data dimensions from input file */
  int klon, nlev;
  int nblocks = (numcols / nproma) + min(numcols % nproma, 1);

  dtype zinfo[4][numthreads];
  const dtype zhpm = 12482329.0;  // IBM P7 HPM flop count for 100 points at L137

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

  tend_loc_t   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_loc_q   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_loc_a   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_loc_cld = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma*nclv );
  tend_cml_t   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_cml_q   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_cml_a   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_cml_cld = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma*nclv );
  tend_tmp_t   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_tmp_q   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_tmp_a   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  tend_tmp_cld = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma*nclv );
  plcrit_aer = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  picrit_aer = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pre_ice    = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pccn       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pnice      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pt         = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pq         = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pvfa       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pvfl       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pvfi       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pdyna      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pdynl      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pdyni      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  phrsw      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  phrlw      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pvervel    = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pap        = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  paph       = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  plsm       = (dtype*) malloc( sizeof(dtype) * nblocks*nproma );
  ldcum      = (int*) malloc( sizeof(int) * nblocks*nproma );
  ktype      = (int*) malloc( sizeof(int) * nblocks*nproma );
  plu        = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  plude      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  psnde      = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pmfu       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pmfd       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pa         = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pclv       = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma*nclv );
  psupsat    = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  pcovptot   = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  prainfrac_toprfz = (dtype*) malloc( sizeof(dtype) * nblocks*nproma );
  pfsqlf    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfsqif    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfcqnng   = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfcqlng   = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfsqrf    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfsqsf    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfcqrng   = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfcqsng   = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfsqltur  = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfsqitur  = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfplsl    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfplsn    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfhpsl    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  pfhpsn    = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );


  dtype rg;
  dtype rd;
  dtype rcpd;
  dtype retv;
  dtype rlvtt;
  dtype rlstt;
  dtype rlmlt;
  dtype rtt;
  dtype rv;
  dtype r2es;
  dtype r3les;
  dtype r3ies;
  dtype r4les;
  dtype r4ies;
  dtype r5les;
  dtype r5ies;
  dtype r5alvcp;
  dtype r5alscp;
  dtype ralvdcp;
  dtype ralsdcp;
  dtype ralfdcp;
  dtype rtwat;
  dtype rtice;
  dtype rticecu;
  dtype rtwat_rtice_r;
  dtype rtwat_rticecu_r;
  dtype rkoop1;
  dtype rkoop2;

  // device declarations
  dtype *d_plcrit_aer;
  dtype *d_picrit_aer;
  dtype *d_pre_ice;
  dtype *d_pccn;
  dtype *d_pnice;
  dtype *d_pt;
  dtype *d_pq;
  dtype *d_tend_loc_t;
  dtype *d_tend_loc_q;
  dtype *d_tend_loc_a;
  dtype *d_tend_loc_cld;
  dtype *d_tend_tmp_t;
  dtype *d_tend_tmp_q;
  dtype *d_tend_tmp_a;
  dtype *d_tend_tmp_cld;
  dtype *d_pvfa;
  dtype *d_pvfl;
  dtype *d_pvfi;
  dtype *d_pdyna;
  dtype *d_pdynl;
  dtype *d_pdyni;
  dtype *d_phrsw;
  dtype *d_phrlw;
  dtype *d_pvervel;
  dtype *d_pap;
  dtype *d_paph;
  dtype *d_plsm;
  int *d_ktype;
  dtype *d_plu;
  dtype *d_plude;
  dtype *d_psnde;
  dtype *d_pmfu;
  dtype *d_pmfd;
  dtype *d_pa;
  dtype *d_pclv;
  dtype *d_psupsat;
  struct TECLDP *d_yrecldp;
  dtype *d_pcovptot;
  dtype *d_prainfrac_toprfz;
  dtype *d_pfsqlf;
  dtype *d_pfsqif;
  dtype *d_pfcqnng;
  dtype *d_pfcqlng;
  dtype *d_pfsqrf;
  dtype *d_pfsqsf;
  dtype *d_pfcqrng;
  dtype *d_pfcqsng;
  dtype *d_pfsqltur;
  dtype *d_pfsqitur;
  dtype *d_pfplsl;
  dtype *d_pfplsn;
  dtype *d_pfhpsl;
  dtype *d_pfhpsn;
  dtype *d_zfoealfa;
  dtype *d_ztp1;
  dtype *d_zli;
  dtype *d_za;
  dtype *d_zaorig;
  dtype *d_zliqfrac;
  dtype *d_zicefrac;
  dtype *d_zqx;
  dtype *d_zqx0;
  dtype *d_zpfplsx;
  dtype *d_zlneg;
  dtype *d_zqxn2d;
  dtype *d_zqsmix;
  dtype *d_zqsliq;
  dtype *d_zqsice;
  dtype *d_zfoeewmt;
  dtype *d_zfoeew;
  dtype *d_zfoeeliqt;
  // end device declarations

  //
  cudaMalloc(&d_plcrit_aer, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_picrit_aer, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pre_ice, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pccn, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pnice, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pt, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pq, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_loc_t, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_loc_q, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_loc_a, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_loc_cld, sizeof(dtype) * nblocks*nlev*nproma*nclv);
  cudaMalloc(&d_tend_tmp_t, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_tmp_q, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_tmp_a, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_tend_tmp_cld, sizeof(dtype) * nblocks*nlev*nproma*nclv);
  cudaMalloc(&d_pvfa, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pvfl, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pvfi, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pdyna, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pdynl, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pdyni, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_phrsw, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_phrlw, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pvervel, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pap, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_paph, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_plsm, sizeof(dtype) * nblocks*nproma);
  cudaMalloc(&d_ktype, sizeof(int) * nblocks*nproma);
  cudaMalloc(&d_plu, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_plude, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_psnde, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pmfu, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pmfd, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pa, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_pclv, sizeof(dtype) * nblocks*nlev*nproma*nclv);
  cudaMalloc(&d_psupsat, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_yrecldp, sizeof(struct TECLDP));
  cudaMalloc(&d_pcovptot, sizeof(dtype) * nblocks*nlev*nproma);
  cudaMalloc(&d_prainfrac_toprfz, sizeof(dtype) * nblocks*nproma);
  cudaMalloc(&d_pfsqlf, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfsqif, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfcqnng, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfcqlng, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfsqrf, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfsqsf, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfcqrng, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfcqsng, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfsqltur, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfsqitur, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfplsl, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfplsn, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfhpsl, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_pfhpsn, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zfoealfa, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_ztp1, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zli, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_za, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zaorig, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zliqfrac, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zicefrac, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zqx, sizeof(dtype) * nblocks*(nlev+1)*nproma*nclv);
  cudaMalloc(&d_zqx0, sizeof(dtype) * nblocks*(nlev+1)*nproma*nclv);
  cudaMalloc(&d_zpfplsx, sizeof(dtype) * nblocks*(nlev+1)*nproma*nclv);
  cudaMalloc(&d_zlneg, sizeof(dtype) * nblocks*(nlev+1)*nproma*nclv);
  cudaMalloc(&d_zqxn2d, sizeof(dtype) * nblocks*(nlev+1)*nproma*nclv);
  cudaMalloc(&d_zqsmix, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zqsliq, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zqsice, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zfoeewmt, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zfoeew, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  cudaMalloc(&d_zfoeeliqt, sizeof(dtype) * nblocks*(nlev+1)*nproma);
  //

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

  double t1 = omp_get_wtime();

  // host to device
  cudaMemcpy(d_plcrit_aer, plcrit_aer, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_picrit_aer, picrit_aer, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pre_ice, pre_ice, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pccn, pccn, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pnice, pnice, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pt, pt, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pq, pq, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_loc_t, tend_loc_t, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_loc_q, tend_loc_q, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_loc_a, tend_loc_a, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_loc_cld, tend_loc_cld, sizeof(dtype) * nblocks*nlev*nproma*nclv, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_tmp_t, tend_tmp_t, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_tmp_q, tend_tmp_q, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_tmp_a, tend_tmp_a, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_tend_tmp_cld, tend_tmp_cld, sizeof(dtype) * nblocks*nlev*nproma*nclv, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pvfa, pvfa, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pvfl, pvfl, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pvfi, pvfi, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pdyna, pdyna, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pdynl, pdynl, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pdyni, pdyni, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_phrsw, phrsw, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_phrlw, phrlw, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pvervel, pvervel, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pap, pap, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_paph, paph, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_plsm, plsm, sizeof(dtype) * nblocks*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_ktype, ktype, sizeof(int) * nblocks*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_plu, plu, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_plude, plude, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_psnde, psnde, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pmfu, pmfu, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pmfd, pmfd, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pa, pa, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pclv, pclv, sizeof(dtype) * nblocks*nlev*nproma*nclv, cudaMemcpyHostToDevice);
  cudaMemcpy(d_psupsat, psupsat, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyHostToDevice);
  cudaMemcpy(d_yrecldp, yrecldp, sizeof(TECLDP), cudaMemcpyHostToDevice);
  // end host to device

  int b, bsize, icalls=0, igpc=numcols;
  int coreid = mycpu();
  int tid = omp_get_thread_num();
  double start = omp_get_wtime();

  dim3 blockdim(nproma, 1, 1);
  dim3 griddim(ceil(((dtype)numcols) / ((dtype)nproma)), 1, 1);
  int jkglo = 0;
  int ibl = (jkglo - 1) / nproma + 1;
  int icend = min(nproma, numcols - jkglo + 1);

  cloudsc_c<<<griddim, blockdim>>>(1, icend, nproma, nlev, ptsphy, d_pt, d_pq,
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
                rkoop1, rkoop2,
                d_zfoealfa, d_ztp1, d_zli,
                d_za, d_zaorig, d_zliqfrac,
                d_zicefrac, d_zqx, d_zqx0,
                d_zpfplsx, d_zlneg, d_zqxn2d,
                d_zqsmix, d_zqsliq, d_zqsice,
                d_zfoeewmt, d_zfoeew, d_zfoeeliqt);

  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );

  double end = omp_get_wtime();

  // device to host
  cudaMemcpy(tend_loc_t, d_tend_loc_t, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(tend_loc_q, d_tend_loc_q, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(tend_loc_a, d_tend_loc_a, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(tend_loc_cld, d_tend_loc_cld, sizeof(dtype) * nblocks*nlev*nproma*nclv, cudaMemcpyDeviceToHost);
  cudaMemcpy(phrlw, d_phrlw, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(plude, d_plude, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(yrecldp, d_yrecldp, sizeof(TECLDP), cudaMemcpyDeviceToHost);
  cudaMemcpy(pcovptot, d_pcovptot, sizeof(dtype) * nblocks*nlev*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(prainfrac_toprfz, d_prainfrac_toprfz, sizeof(dtype) * nblocks*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqlf, d_pfsqlf, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqif, d_pfsqif, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfcqnng, d_pfcqnng, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfcqlng, d_pfcqlng, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqrf, d_pfsqrf, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqsf, d_pfsqsf, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfcqrng, d_pfcqrng, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfcqsng, d_pfcqsng, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqltur, d_pfsqltur, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfsqitur, d_pfsqitur, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfplsl, d_pfplsl, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfplsn, d_pfplsn, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfhpsl, d_pfhpsl, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);
  cudaMemcpy(pfhpsn, d_pfhpsn, sizeof(dtype) * nblocks*(nlev+1)*nproma, cudaMemcpyDeviceToHost);

  /* int msec = diff * 1000 / CLOCKS_PER_SEC; */
  zinfo[0][tid] = end - start;
  zinfo[1][tid] = (dtype) coreid;
  zinfo[2][tid] = (dtype) icalls;
  zinfo[3][tid] = (dtype) igpc;

  double t2 = omp_get_wtime();

  printf("     NUMOMP=%d, NGPTOT=%d, NPROMA=%d, NGPBLKS=%d\n", numthreads, numcols, nproma, nblocks);
  printf(" %10s%10s%10s%10s%10s %+4s : %10s%10s%10s\n",
    "NUMOMP", "NGPTOT", "#GP-cols", "#BLKS", "NPROMA", "tid#", "Time(msec)", "MFlops/s", "col/s");
  double zfrac, zmflops, zthrput;
  for (int t = 0; t < numthreads; t++) {
    const dtype tloc = zinfo[0][t];
    const int coreid = (int) zinfo[1][t];
    const int icalls = (int) zinfo[2][t];
    const int igpc = (int) zinfo[3][t];
    zfrac = (double)igpc / (double)numcols;
    if (tloc > 0.0) {
      zmflops = 1.0e-06 * zfrac * zhpm * ((double)numcols / 100.) / tloc;
      zthrput = (dtype)numcols/tloc;
    } else {
      zmflops = 0.;
      zthrput = 0.0;
    }
    printf(" %10d%10d%10d%10d%10d %4d : %10d%10d%10d @ core#\n",
           numthreads, numcols, igpc, icalls, nproma, t, (int)(tloc * 1000.), (int)zmflops, (int)zthrput);
  }
  double tdiff = t2 - t1;
  zfrac = 1.0;
  if (tdiff > 0.0) {
    zmflops = 1.0e-06 * zfrac * zhpm * ((double)numcols / 100.) / tdiff;
    zthrput = (dtype)numcols/tdiff;
  } else {
    zmflops = 0.0;
    zthrput = 0.0;
  }
  printf(" %10d%10d%10d%10d%10d %4d : %10d%10d%10d TOTAL\n",
         numthreads, numcols, numcols, nblocks, nproma, -1, (int)(tdiff * 1000.), (int)zmflops, (int)zthrput);

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
  cudaFree(d_plcrit_aer);
  cudaFree(d_picrit_aer);
  cudaFree(d_pre_ice);
  cudaFree(d_pccn);
  cudaFree(d_pnice);
  cudaFree(d_pt);
  cudaFree(d_pq);
  cudaFree(d_tend_loc_t);
  cudaFree(d_tend_loc_q);
  cudaFree(d_tend_loc_a);
  cudaFree(d_tend_loc_cld);
  cudaFree(d_tend_tmp_t);
  cudaFree(d_tend_tmp_q);
  cudaFree(d_tend_tmp_a);
  cudaFree(d_tend_tmp_cld);
  cudaFree(d_pvfa);
  cudaFree(d_pvfl);
  cudaFree(d_pvfi);
  cudaFree(d_pdyna);
  cudaFree(d_pdynl);
  cudaFree(d_pdyni);
  cudaFree(d_phrsw);
  cudaFree(d_phrlw);
  cudaFree(d_pvervel);
  cudaFree(d_pap);
  cudaFree(d_paph);
  cudaFree(d_plsm);
  cudaFree(d_ktype);
  cudaFree(d_plu);
  cudaFree(d_plude);
  cudaFree(d_psnde);
  cudaFree(d_pmfu);
  cudaFree(d_pmfd);
  cudaFree(d_pa);
  cudaFree(d_pclv);
  cudaFree(d_psupsat);
  cudaFree(d_yrecldp);
  cudaFree(d_pcovptot);
  cudaFree(d_prainfrac_toprfz);
  cudaFree(d_pfsqlf);
  cudaFree(d_pfsqif);
  cudaFree(d_pfcqnng);
  cudaFree(d_pfcqlng);
  cudaFree(d_pfsqrf);
  cudaFree(d_pfsqsf);
  cudaFree(d_pfcqrng);
  cudaFree(d_pfcqsng);
  cudaFree(d_pfsqltur);
  cudaFree(d_pfsqitur);
  cudaFree(d_pfplsl);
  cudaFree(d_pfplsn);
  cudaFree(d_pfhpsl);
  cudaFree(d_pfhpsn);
  cudaFree(d_zfoealfa);
  cudaFree(d_ztp1);
  cudaFree(d_zli);
  cudaFree(d_za);
  cudaFree(d_zaorig);
  cudaFree(d_zliqfrac);
  cudaFree(d_zicefrac);
  cudaFree(d_zqx);
  cudaFree(d_zqx0);
  cudaFree(d_zpfplsx);
  cudaFree(d_zlneg);
  cudaFree(d_zqxn2d);
  cudaFree(d_zqsmix);
  cudaFree(d_zqsliq);
  cudaFree(d_zqsice);
  cudaFree(d_zfoeewmt);
  cudaFree(d_zfoeew);
  cudaFree(d_zfoeeliqt);
  // end free device
}
