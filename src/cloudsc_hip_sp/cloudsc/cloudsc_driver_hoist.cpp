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
inline void gpuAssert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

void cloudsc_driver(int numthreads, int numcols, int nproma) {

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
  float *d_zfoealfa;
  float *d_ztp1; 
  float *d_zli;
  float *d_za; 
  float *d_zaorig; 
  float *d_zliqfrac;
  float *d_zicefrac; 
  float *d_zqx; 
  float *d_zqx0;
  float *d_zpfplsx; 
  float *d_zlneg;
  float *d_zqxn2d;
  float *d_zqsmix; 
  float *d_zqsliq; 
  float *d_zqsice;
  float *d_zfoeewmt; 
  float *d_zfoeew; 
  float *d_zfoeeliqt;
  // end device declarations

  //
  hipMalloc(&d_plcrit_aer, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_picrit_aer, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pre_ice, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pccn, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pnice, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pt, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pq, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_t, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_q, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_a, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_cld, sizeof(float) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_tend_tmp_t, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_q, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_a, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_cld, sizeof(float) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_pvfa, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pvfl, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pvfi, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pdyna, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pdynl, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pdyni, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_phrsw, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_phrlw, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pvervel, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pap, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_paph, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_plsm, sizeof(float) * nblocks*nproma);
  hipMalloc(&d_ktype, sizeof(int) * nblocks*nproma);
  hipMalloc(&d_plu, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_plude, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_psnde, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pmfu, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pmfd, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pa, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_pclv, sizeof(float) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_psupsat, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_yrecldp, sizeof(struct TECLDP));
  hipMalloc(&d_pcovptot, sizeof(float) * nblocks*nlev*nproma);
  hipMalloc(&d_prainfrac_toprfz, sizeof(float) * nblocks*nproma);
  hipMalloc(&d_pfsqlf, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqif, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqnng, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqlng, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqrf, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqsf, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqrng, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqsng, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqltur, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqitur, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfplsl, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfplsn, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfhpsl, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfhpsn, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zfoealfa, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_ztp1, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zli, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_za, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zaorig, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zliqfrac, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zicefrac, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zqx, sizeof(float) * nblocks*(nlev+1)*nproma*nclv);
  hipMalloc(&d_zqx0, sizeof(float) * nblocks*(nlev+1)*nproma*nclv);
  hipMalloc(&d_zpfplsx, sizeof(float) * nblocks*(nlev+1)*nproma*nclv);
  hipMalloc(&d_zlneg, sizeof(float) * nblocks*(nlev+1)*nproma*nclv);
  hipMalloc(&d_zqxn2d, sizeof(float) * nblocks*(nlev+1)*nproma*nclv);
  hipMalloc(&d_zqsmix, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zqsliq, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zqsice, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zfoeewmt, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zfoeew, sizeof(float) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_zfoeeliqt, sizeof(float) * nblocks*(nlev+1)*nproma);
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

  // host to device
  hipMemcpy(d_plcrit_aer, plcrit_aer, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_picrit_aer, picrit_aer, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pre_ice, pre_ice, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pccn, pccn, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pnice, pnice, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pt, pt, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pq, pq, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_t, tend_loc_t, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_q, tend_loc_q, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_a, tend_loc_a, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_cld, tend_loc_cld, sizeof(float) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_t, tend_tmp_t, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_q, tend_tmp_q, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_a, tend_tmp_a, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_cld, tend_tmp_cld, sizeof(float) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfa, pvfa, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfl, pvfl, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfi, pvfi, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdyna, pdyna, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdynl, pdynl, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdyni, pdyni, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_phrsw, phrsw, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_phrlw, phrlw, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvervel, pvervel, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pap, pap, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_paph, paph, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plsm, plsm, sizeof(float) * nblocks*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_ktype, ktype, sizeof(int) * nblocks*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plu, plu, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plude, plude, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_psnde, psnde, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pmfu, pmfu, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pmfd, pmfd, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pa, pa, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pclv, pclv, sizeof(float) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_psupsat, psupsat, sizeof(float) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_yrecldp, yrecldp, sizeof(TECLDP), hipMemcpyHostToDevice);
  // end host to device

  double t1 = omp_get_wtime();

    int b, bsize, icalls=0, igpc=numcols;
    int coreid = mycpu();
    int tid = omp_get_thread_num();
    double start = omp_get_wtime();

    dim3 blockdim(nproma, 1, 1);
    //dim3 griddim(1, 1, ceil(((float)numcols) / ((float)nproma)));
    dim3 griddim(ceil(((float)numcols) / ((float)nproma)), 1, 1);
    int jkglo = 0;
    int ibl = (jkglo - 1) / nproma + 1;
    int icend = min(nproma, numcols - jkglo + 1);

    cloudsc_c<<<griddim, blockdim>>>(1, icend/*bsize*/, nproma/*, nlev*/, ptsphy, d_pt, d_pq,
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

    gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );

    double end = omp_get_wtime();

    // device to host
    hipMemcpy(tend_loc_t, d_tend_loc_t, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_q, d_tend_loc_q, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_a, d_tend_loc_a, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_cld, d_tend_loc_cld, sizeof(float) * nblocks*nlev*nproma*nclv, hipMemcpyDeviceToHost);
    hipMemcpy(phrlw, d_phrlw, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(plude, d_plude, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(yrecldp, d_yrecldp, sizeof(TECLDP), hipMemcpyDeviceToHost);
    hipMemcpy(pcovptot, d_pcovptot, sizeof(float) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(prainfrac_toprfz, d_prainfrac_toprfz, sizeof(float) * nblocks*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqlf, d_pfsqlf, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqif, d_pfsqif, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqnng, d_pfcqnng, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqlng, d_pfcqlng, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqrf, d_pfsqrf, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqsf, d_pfsqsf, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqrng, d_pfcqrng, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqsng, d_pfcqsng, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqltur, d_pfsqltur, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqitur, d_pfsqitur, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfplsl, d_pfplsl, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfplsn, d_pfplsn, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfhpsl, d_pfhpsl, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfhpsn, d_pfhpsn, sizeof(float) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);

    /* int msec = diff * 1000 / CLOCKS_PER_SEC; */
    zinfo[0][tid] = end - start;
    zinfo[1][tid] = (float) coreid;
    zinfo[2][tid] = (float) icalls;
    zinfo[3][tid] = (float) igpc;

  double t2 = omp_get_wtime();

  printf("     NUMOMP=%d, NGPTOT=%d, NPROMA=%d, NGPBLKS=%d\n", numthreads, numcols, nproma, nblocks);
  printf(" %+10s%+10s%+10s%+10s%+10s %+4s : %+10s%+10s%+10s\n",
    "NUMOMP", "NGPTOT", "#GP-cols", "#BLKS", "NPROMA", "tid#", "Time(msec)", "MFlops/s", "col/s");
  double zfrac, zmflops, zthrput;
  for (int t = 0; t < numthreads; t++) {
    const float tloc = zinfo[0][t];
    const int coreid = (int) zinfo[1][t];
    const int icalls = (int) zinfo[2][t];
    const int igpc = (int) zinfo[3][t];
    zfrac = (double)igpc / (double)numcols;
    if (tloc > 0.0) {
      zmflops = 1.0e-06 * zfrac * zhpm * ((double)numcols / 100.) / tloc;
      zthrput = (double)numcols/tloc;
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
    zthrput = (double)numcols/tdiff;
  } else {
    zmflops = 0.0;
    zthrput = 0.0;
  }
  printf(" %10d%10d%10d%10d%10d %4d: %10d%10d%10d TOTAL\n",
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
  hipFree(d_plcrit_aer);
  hipFree(d_picrit_aer);
  hipFree(d_pre_ice);
  hipFree(d_pccn);
  hipFree(d_pnice);
  hipFree(d_pt);
  hipFree(d_pq);
  hipFree(d_tend_loc_t);
  hipFree(d_tend_loc_q);
  hipFree(d_tend_loc_a);
  hipFree(d_tend_loc_cld);
  hipFree(d_tend_tmp_t);
  hipFree(d_tend_tmp_q);
  hipFree(d_tend_tmp_a);
  hipFree(d_tend_tmp_cld);
  hipFree(d_pvfa);
  hipFree(d_pvfl);
  hipFree(d_pvfi);
  hipFree(d_pdyna);
  hipFree(d_pdynl);
  hipFree(d_pdyni);
  hipFree(d_phrsw);
  hipFree(d_phrlw);
  hipFree(d_pvervel);
  hipFree(d_pap);
  hipFree(d_paph);
  hipFree(d_plsm);
  hipFree(d_ktype);
  hipFree(d_plu);
  hipFree(d_plude);
  hipFree(d_psnde);
  hipFree(d_pmfu);
  hipFree(d_pmfd);
  hipFree(d_pa);
  hipFree(d_pclv);
  hipFree(d_psupsat);
  hipFree(d_yrecldp);
  hipFree(d_pcovptot);
  hipFree(d_prainfrac_toprfz);
  hipFree(d_pfsqlf);
  hipFree(d_pfsqif);
  hipFree(d_pfcqnng);
  hipFree(d_pfcqlng);
  hipFree(d_pfsqrf);
  hipFree(d_pfsqsf);
  hipFree(d_pfcqrng);
  hipFree(d_pfcqsng);
  hipFree(d_pfsqltur);
  hipFree(d_pfsqitur);
  hipFree(d_pfplsl);
  hipFree(d_pfplsn);
  hipFree(d_pfhpsl);
  hipFree(d_pfhpsn);
  hipFree(d_zfoealfa);
  hipFree(d_ztp1);
  hipFree(d_zli);
  hipFree(d_za);
  hipFree(d_zaorig);
  hipFree(d_zliqfrac);
  hipFree(d_zicefrac);
  hipFree(d_zqx);
  hipFree(d_zqx0);
  hipFree(d_zpfplsx);
  hipFree(d_zlneg);
  hipFree(d_zqxn2d);
  hipFree(d_zqsmix);
  hipFree(d_zqsliq);
  hipFree(d_zqsice);
  hipFree(d_zfoeewmt);
  hipFree(d_zfoeew);
  hipFree(d_zfoeeliqt);
  // end free device
}

