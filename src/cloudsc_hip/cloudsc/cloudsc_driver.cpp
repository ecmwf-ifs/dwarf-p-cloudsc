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

  double *tend_tmp_u;
  double *tend_tmp_v;
  double *tend_tmp_t;
  double *tend_tmp_q;
  double *tend_tmp_o3;
  double *tend_tmp_a;
  double *tend_tmp_cld;

  double *tend_loc_u;
  double *tend_loc_v;
  double *tend_loc_t;
  double *tend_loc_q;
  double *tend_loc_o3;
  double *tend_loc_a;
  double *tend_loc_cld;

  double *tend_cml_u;
  double *tend_cml_v;
  double *tend_cml_t;
  double *tend_cml_q;
  double *tend_cml_o3;
  double *tend_cml_a;
  double *tend_cml_cld;

  double ptsphy;            //! Physics timestep

  double *plcrit_aer;
  double *picrit_aer;
  double *pre_ice;
  double *pccn;       //! liquid cloud condensation nuclei
  double *pnice;     //! ice number concentration (cf. CCN)
  double *pt;           //! T at start of callpar
  double *pq;           //! Q at start of callpar
  double *pvfa;       //! CC from VDF scheme
  double *pvfl;       //! Liq from VDF scheme
  double *pvfi;       //! Ice from VDF scheme
  double *pdyna;     //! CC from Dynamics
  double *pdynl;     //! Liq from Dynamics
  double *pdyni;     //! Liq from Dynamics
  double *phrsw;     //! Short-wave heating rate
  double *phrlw;     //! Long-wave heating rate
  double *pvervel; //! Vertical velocity
  double *pap;         //! Pressure on full levels
  double *paph;       //! Pressure on half levels
  double *plsm;           //! Land fraction (0-1)
  int *ldcum;
  int    *ktype;          //! Convection type 0,1,2
  double *plu;         //! Conv. condensate
  double *plude;     //! Conv. detrained water
  double *plude_tmp;
  double *psnde;     //! Conv. detrained snow
  double *pmfu;       //! Conv. mass flux up
  double *pmfd;       //! Conv. mass flux down
  double *pa;           //! Original Cloud fraction (t)

  double *pclv;
  double *psupsat;

  double *pcovptot;   //! Precip fraction
  double *prainfrac_toprfz;
  double *pfsqlf;       //! Flux of liquid
  double *pfsqif;       //! Flux of ice
  double *pfcqlng;     //! -ve corr for liq
  double *pfcqnng;     //! -ve corr for ice
  double *pfsqrf;       //! Flux diagnostics
  double *pfsqsf;       //!    for DDH, generic
  double *pfcqrng;     //! rain
  double *pfcqsng;     //! snow
  double *pfsqltur;   //! liquid flux due to VDF
  double *pfsqitur;   //! ice flux due to VDF
  double *pfplsl;       //! liq+rain sedim flux
  double *pfplsn;       //! ice+snow sedim flux
  double *pfhpsl;       //! Enthalpy flux for liq
  double *pfhpsn;       //! Enthalpy flux for ice

  /* Define or query data dimensions from input file */
  int klon, nlev;
  int nblocks = (numcols / nproma) + min(numcols % nproma, 1);

  double zinfo[4][numthreads];
  const double zhpm = 12482329.0;  // IBM P7 HPM flop count for 100 points at L137

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

  tend_loc_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_loc_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_loc_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_loc_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );
  tend_cml_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_cml_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_cml_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_cml_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );
  tend_tmp_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_tmp_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_tmp_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  tend_tmp_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );

  plcrit_aer = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  picrit_aer = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pre_ice    = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pccn       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pnice      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pt         = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pq         = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pvfa       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pvfl       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pvfi       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pdyna      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pdynl      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pdyni      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  phrsw      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  phrlw      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pvervel    = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pap        = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  paph       = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  plsm       = (double*) malloc( sizeof(double) * nblocks*nproma );
  ldcum      = (int*) malloc( sizeof(int) * nblocks*nproma );
  ktype      = (int*) malloc( sizeof(int) * nblocks*nproma );
  plu        = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  plude      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  psnde      = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pmfu       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pmfd       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pa         = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pclv       = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );
  psupsat    = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  pcovptot   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  prainfrac_toprfz = (double*) malloc( sizeof(double) * nblocks*nproma );
  pfsqlf    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfsqif    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfcqnng   = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfcqlng   = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfsqrf    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfsqsf    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfcqrng   = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfcqsng   = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfsqltur  = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfsqitur  = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfplsl    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfplsn    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfhpsl    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  pfhpsn    = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );


  double rg;
  double rd;
  double rcpd;
  double retv;
  double rlvtt;
  double rlstt;
  double rlmlt;
  double rtt;
  double rv;
  double r2es;
  double r3les;
  double r3ies;
  double r4les;
  double r4ies;
  double r5les;
  double r5ies;
  double r5alvcp;
  double r5alscp;
  double ralvdcp;
  double ralsdcp;
  double ralfdcp;
  double rtwat;
  double rtice;
  double rticecu;
  double rtwat_rtice_r;
  double rtwat_rticecu_r;
  double rkoop1;
  double rkoop2;

  // device declarations
  double *d_plcrit_aer;
  double *d_picrit_aer;
  double *d_pre_ice;
  double *d_pccn;
  double *d_pnice;
  double *d_pt;
  double *d_pq;
  double *d_tend_loc_t;
  double *d_tend_loc_q;
  double *d_tend_loc_a;
  double *d_tend_loc_cld;
  double *d_tend_tmp_t;
  double *d_tend_tmp_q;
  double *d_tend_tmp_a;
  double *d_tend_tmp_cld;
  double *d_pvfa;
  double *d_pvfl;
  double *d_pvfi;
  double *d_pdyna;
  double *d_pdynl;
  double *d_pdyni;
  double *d_phrsw;
  double *d_phrlw;
  double *d_pvervel;
  double *d_pap;
  double *d_paph;
  double *d_plsm;
  int *d_ktype;
  double *d_plu;
  double *d_plude;
  double *d_psnde;
  double *d_pmfu;
  double *d_pmfd;
  double *d_pa;
  double *d_pclv;
  double *d_psupsat;
  struct TECLDP *d_yrecldp;
  double *d_pcovptot;
  double *d_prainfrac_toprfz;
  double *d_pfsqlf;
  double *d_pfsqif;
  double *d_pfcqnng;
  double *d_pfcqlng;
  double *d_pfsqrf;
  double *d_pfsqsf;
  double *d_pfcqrng;
  double *d_pfcqsng;
  double *d_pfsqltur;
  double *d_pfsqitur;
  double *d_pfplsl;
  double *d_pfplsn;
  double *d_pfhpsl;
  double *d_pfhpsn;
  // end device declarations

  hipMalloc(&d_plcrit_aer, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_picrit_aer, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pre_ice, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pccn, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pnice, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pt, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pq, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_t, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_q, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_a, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_loc_cld, sizeof(double) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_tend_tmp_t, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_q, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_a, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_tend_tmp_cld, sizeof(double) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_pvfa, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pvfl, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pvfi, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pdyna, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pdynl, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pdyni, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_phrsw, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_phrlw, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pvervel, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pap, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_paph, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_plsm, sizeof(double) * nblocks*nproma);
  hipMalloc(&d_ktype, sizeof(int) * nblocks*nproma);
  hipMalloc(&d_plu, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_plude, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_psnde, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pmfu, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pmfd, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pa, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_pclv, sizeof(double) * nblocks*nlev*nproma*nclv);
  hipMalloc(&d_psupsat, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_yrecldp, sizeof(struct TECLDP));
  hipMalloc(&d_pcovptot, sizeof(double) * nblocks*nlev*nproma);
  hipMalloc(&d_prainfrac_toprfz, sizeof(double) * nblocks*nproma);
  hipMalloc(&d_pfsqlf, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqif, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqnng, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqlng, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqrf, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqsf, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqrng, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfcqsng, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqltur, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfsqitur, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfplsl, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfplsn, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfhpsl, sizeof(double) * nblocks*(nlev+1)*nproma);
  hipMalloc(&d_pfhpsn, sizeof(double) * nblocks*(nlev+1)*nproma);

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
  hipMemcpy(d_plcrit_aer, plcrit_aer, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_picrit_aer, picrit_aer, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pre_ice, pre_ice, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pccn, pccn, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pnice, pnice, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pt, pt, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pq, pq, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_t, tend_loc_t, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_q, tend_loc_q, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_a, tend_loc_a, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_loc_cld, tend_loc_cld, sizeof(double) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_t, tend_tmp_t, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_q, tend_tmp_q, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_a, tend_tmp_a, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_tend_tmp_cld, tend_tmp_cld, sizeof(double) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfa, pvfa, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfl, pvfl, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvfi, pvfi, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdyna, pdyna, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdynl, pdynl, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pdyni, pdyni, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_phrsw, phrsw, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_phrlw, phrlw, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pvervel, pvervel, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pap, pap, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_paph, paph, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plsm, plsm, sizeof(double) * nblocks*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_ktype, ktype, sizeof(int) * nblocks*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plu, plu, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_plude, plude, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_psnde, psnde, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pmfu, pmfu, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pmfd, pmfd, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pa, pa, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_pclv, pclv, sizeof(double) * nblocks*nlev*nproma*nclv, hipMemcpyHostToDevice);
  hipMemcpy(d_psupsat, psupsat, sizeof(double) * nblocks*nlev*nproma, hipMemcpyHostToDevice);
  hipMemcpy(d_yrecldp, yrecldp, sizeof(TECLDP), hipMemcpyHostToDevice);
  // end host to device

  double t1 = omp_get_wtime();

    int b, bsize, icalls=0, igpc=numcols;
    int coreid = mycpu();
    int tid = omp_get_thread_num();
    double start = omp_get_wtime();

    dim3 blockdim(nproma, 1, 1);
    //dim3 griddim(1, 1, ceil(((double)numcols) / ((double)nproma)));
    dim3 griddim(ceil(((double)numcols) / ((double)nproma)), 1, 1);
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
                rkoop1, rkoop2);


    gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );

    double end = omp_get_wtime();

    // device to host
    hipMemcpy(tend_loc_t, d_tend_loc_t, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_q, d_tend_loc_q, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_a, d_tend_loc_a, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(tend_loc_cld, d_tend_loc_cld, sizeof(double) * nblocks*nlev*nproma*nclv, hipMemcpyDeviceToHost);
    hipMemcpy(phrlw, d_phrlw, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(plude, d_plude, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(yrecldp, d_yrecldp, sizeof(TECLDP), hipMemcpyDeviceToHost);
    hipMemcpy(pcovptot, d_pcovptot, sizeof(double) * nblocks*nlev*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(prainfrac_toprfz, d_prainfrac_toprfz, sizeof(double) * nblocks*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqlf, d_pfsqlf, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqif, d_pfsqif, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqnng, d_pfcqnng, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqlng, d_pfcqlng, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqrf, d_pfsqrf, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqsf, d_pfsqsf, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqrng, d_pfcqrng, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfcqsng, d_pfcqsng, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqltur, d_pfsqltur, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfsqitur, d_pfsqitur, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfplsl, d_pfplsl, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfplsn, d_pfplsn, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfhpsl, d_pfhpsl, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    hipMemcpy(pfhpsn, d_pfhpsn, sizeof(double) * nblocks*(nlev+1)*nproma, hipMemcpyDeviceToHost);
    // end device to host

    /* int msec = diff * 1000 / CLOCKS_PER_SEC; */
    zinfo[0][tid] = end - start;
    zinfo[1][tid] = (double) coreid;
    zinfo[2][tid] = (double) icalls;
    zinfo[3][tid] = (double) igpc;

  double t2 = omp_get_wtime();

  printf("     NUMOMP=%d, NGPTOT=%d, NPROMA=%d, NGPBLKS=%d\n", numthreads, numcols, nproma, nblocks);
  printf(" %+10s%+10s%+10s%+10s%+10s %+4s : %+10s%+10s\n",
    "NUMOMP", "NGPTOT", "#GP-cols", "#BLKS", "NPROMA", "tid#", "Time(msec)", "MFlops/s");
  double zfrac, zmflops;
  for (int t = 0; t < numthreads; t++) {
    const double tloc = zinfo[0][t];
    const int coreid = (int) zinfo[1][t];
    const int icalls = (int) zinfo[2][t];
    const int igpc = (int) zinfo[3][t];
    zfrac = (double)igpc / (double)numcols;
    if (tloc > 0.0) {
      zmflops = 1.0e-06 * zfrac * zhpm * ((double)numcols / 100.) / tloc;
    } else {
      zmflops = 0.;
    }
    printf(" %10d%10d%10d%10d%10d %4d : %10d%10d @ core#\n",
	   numthreads, numcols, igpc, icalls, nproma, t, (int)(tloc * 1000.), (int)zmflops);
  }
  double tdiff = t2 - t1;
  zfrac = 1.0;
  if (tdiff > 0.0) {
    zmflops = 1.0e-06 * zfrac * zhpm * ((double)numcols / 100.) / tdiff;
  } else {
    zmflops = 0.0;
  }
  printf(" %10d%10d%10d%10d%10d %4d : %10d%10d TOTAL\n",
	 numthreads, numcols, numcols, nblocks, nproma, -1, (int)(tdiff * 1000.), (int)zmflops);

  cloudsc_validate(klon, nlev, nclv, numcols, nproma,
		   plude, pcovptot, prainfrac_toprfz, pfsqlf, pfsqif,
		   pfcqlng, pfcqnng, pfsqrf, pfsqsf, pfcqrng, pfcqsng,
		   pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn,
		   tend_loc_a, tend_loc_q, tend_loc_t, tend_loc_cld);

  free(plcrit_aer); // ALLOCATE(PLCRIT_AER(KLON,KLEV))
  free(picrit_aer); // ALLOCATE(PICRIT_AER(KLON,KLEV))
  free(pre_ice);    // ALLOCATE(PRE_ICE(KLON,KLEV))
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
  free(paph); // ALLOCATE(PAPH(KLON,KLEV+1))
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
  // end free device
}

