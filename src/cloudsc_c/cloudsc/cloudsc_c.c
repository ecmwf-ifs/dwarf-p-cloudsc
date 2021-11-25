/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "cloudsc_c.h"

#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>


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
	      double * restrict v_pfhpsn)
{
  /* Array casts for pointer arguments */
  double (*pt)[klon] = (double (*)[klon]) v_pt;
  double (*pq)[klon] = (double (*)[klon]) v_pq;
  double (*tendency_cml_t)[klon] = (double (*)[klon]) v_tendency_cml_t;
  double (*tendency_cml_q)[klon] = (double (*)[klon]) v_tendency_cml_q;
  double (*tendency_cml_a)[klon] = (double (*)[klon]) v_tendency_cml_a;
  double (*tendency_cml_cld)[klev][klon] = (double (*)[klev][klon]) v_tendency_cml_cld;
  double (*tendency_tmp_t)[klon] = (double (*)[klon]) v_tendency_tmp_t;
  double (*tendency_tmp_q)[klon] = (double (*)[klon]) v_tendency_tmp_q;
  double (*tendency_tmp_a)[klon] = (double (*)[klon]) v_tendency_tmp_a;
  double (*tendency_tmp_cld)[klev][klon] = (double (*)[klev][klon]) v_tendency_tmp_cld;
  double (*tendency_loc_t)[klon] = (double (*)[klon]) v_tendency_loc_t;
  double (*tendency_loc_q)[klon] = (double (*)[klon]) v_tendency_loc_q;
  double (*tendency_loc_a)[klon] = (double (*)[klon]) v_tendency_loc_a;
  double (*tendency_loc_cld)[klev][klon] = (double (*)[klev][klon]) v_tendency_loc_cld;
  double (*pvfa)[klon] = (double (*)[klon]) v_pvfa;
  double (*pvfl)[klon] = (double (*)[klon]) v_pvfl;
  double (*pvfi)[klon] = (double (*)[klon]) v_pvfi;
  double (*pdyna)[klon] = (double (*)[klon]) v_pdyna;
  double (*pdynl)[klon] = (double (*)[klon]) v_pdynl;
  double (*pdyni)[klon] = (double (*)[klon]) v_pdyni;
  double (*phrsw)[klon] = (double (*)[klon]) v_phrsw;
  double (*phrlw)[klon] = (double (*)[klon]) v_phrlw;
  double (*pvervel)[klon] = (double (*)[klon]) v_pvervel;
  double (*pap)[klon] = (double (*)[klon]) v_pap;
  double (*paph)[klon] = (double (*)[klon]) v_paph;
  double (*plsm) = (double (*)) v_plsm;
  int (*ktype) = (int (*)) v_ktype;
  double (*plu)[klon] = (double (*)[klon]) v_plu;
  double (*plude)[klon] = (double (*)[klon]) v_plude;
  double (*psnde)[klon] = (double (*)[klon]) v_psnde;
  double (*pmfu)[klon] = (double (*)[klon]) v_pmfu;
  double (*pmfd)[klon] = (double (*)[klon]) v_pmfd;
  double (*pa)[klon] = (double (*)[klon]) v_pa;
  double (*pclv)[klev][klon] = (double (*)[klev][klon]) v_pclv;
  double (*psupsat)[klon] = (double (*)[klon]) v_psupsat;
  double (*plcrit_aer)[klon] = (double (*)[klon]) v_plcrit_aer;
  double (*picrit_aer)[klon] = (double (*)[klon]) v_picrit_aer;
  double (*pre_ice)[klon] = (double (*)[klon]) v_pre_ice;
  double (*pccn)[klon] = (double (*)[klon]) v_pccn;
  double (*pnice)[klon] = (double (*)[klon]) v_pnice;
  double (*pcovptot)[klon] = (double (*)[klon]) v_pcovptot;
  double (*prainfrac_toprfz) = (double (*)) v_prainfrac_toprfz;
  double (*pfsqlf)[klon] = (double (*)[klon]) v_pfsqlf;
  double (*pfsqif)[klon] = (double (*)[klon]) v_pfsqif;
  double (*pfcqnng)[klon] = (double (*)[klon]) v_pfcqnng;
  double (*pfcqlng)[klon] = (double (*)[klon]) v_pfcqlng;
  double (*pfsqrf)[klon] = (double (*)[klon]) v_pfsqrf;
  double (*pfsqsf)[klon] = (double (*)[klon]) v_pfsqsf;
  double (*pfcqrng)[klon] = (double (*)[klon]) v_pfcqrng;
  double (*pfcqsng)[klon] = (double (*)[klon]) v_pfcqsng;
  double (*pfsqltur)[klon] = (double (*)[klon]) v_pfsqltur;
  double (*pfsqitur)[klon] = (double (*)[klon]) v_pfsqitur;
  double (*pfplsl)[klon] = (double (*)[klon]) v_pfplsl;
  double (*pfplsn)[klon] = (double (*)[klon]) v_pfplsn;
  double (*pfhpsl)[klon] = (double (*)[klon]) v_pfhpsl;
  double (*pfhpsn)[klon] = (double (*)[klon]) v_pfhpsn;
  double zlcond1[klon];
  double zlcond2[klon];
  double zlevap;
  double zleros;
  double zlevapl[klon];
  double zlevapi[klon];
  double zrainaut[klon];
  double zsnowaut[klon];
  double zliqcld[klon];
  double zicecld[klon];
  double zfokoop[klon];
  double zfoealfa[klev+1][klon];
  double zicenuclei[klon];
  double zlicld[klon];
  double zacond;
  double zaeros;
  double zlfinalsum[klon];
  double zdqs[klon];
  double ztold[klon];
  double zqold[klon];
  double zdtgdp[klon];
  double zrdtgdp[klon];
  double ztrpaus[klon];
  double zcovpclr[klon];
  double zpreclr;
  double zcovptot[klon];
  double zcovpmax[klon];
  double zqpretot[klon];
  double zdpevap;
  double zdtforc;
  double zdtdiab;
  double ztp1[klev][klon];
  double zldefr[klon];
  double zldifdt[klon];
  double zdtgdpf[klon];
  double zlcust[5][klon];
  double zacust[klon];
  double zmf[klon];
  double zrho[klon];
  double ztmp1[klon];
  double ztmp2[klon];
  double ztmp3[klon];
  double ztmp4[klon];
  double ztmp5[klon];
  double ztmp6[klon];
  double ztmp7[klon];
  double zalfawm[klon];
  double zsolab[klon];
  double zsolac[klon];
  double zanew;
  double zanewm1[klon];
  double zgdp[klon];
  double zda[klon];
  double zli[klev][klon];
  double za[klev][klon];
  double zaorig[klev][klon];
  int llflag[klon];
  int llo1;
  int icall;
  int ik;
  int jk;
  int jl;
  int jm;
  int jn;
  int jo;
  int jlen;
  int is;
  double zdp[klon];
  double zpaphd[klon];
  double zalfa;
  double zalfaw;
  double zbeta;
  double zbeta1;
  double zcfpr;
  double zcor;
  double zcdmax;
  double zmin[klon];
  double zlcondlim;
  double zdenom;
  double zdpmxdt;
  double zdpr;
  double zdtdp;
  double ze;
  double zepsec;
  double zfac;
  double zfaci;
  double zfacw;
  double zgdcp;
  double zinew;
  double zlcrit;
  double zmfdn;
  double zprecip;
  double zqe;
  double zqsat;
  double zqtmst;
  double zrdcp;
  double zrhc;
  double zsig;
  double zsigk;
  double zwtot;
  double zzco;
  double zzdl;
  double zzrh;
  double zzzdt;
  double zqadj;
  double zqnew;
  double ztnew;
  double zrg_r;
  double zgdph_r;
  double zcons1;
  double zcond;
  double zcons1a;
  double zlfinal;
  double zmelt;
  double zevap;
  double zfrz;
  double zvpliq;
  double zvpice;
  double zadd;
  double zbdd;
  double zcvds;
  double zice0;
  double zdepos;
  double zsupsat[klon];
  double zfall;
  double zre_ice;
  double zrldcp;
  double zqp1env;
  int iphase[5];
  int imelt[5];
  int llfall[5];
  int llindex1[5][klon];
  int llindex3[5][5][klon];
  double zmax;
  double zrat;
  int iorder[5][klon];
  double zliqfrac[klev][klon];
  double zicefrac[klev][klon];
  double zqx[5][klev][klon];
  double zqx0[5][klev][klon];
  double zqxn[5][klon];
  double zqxfg[5][klon];
  double zqxnm1[5][klon];
  double zfluxq[5][klon];
  double zpfplsx[5][klev+1][klon];
  double zlneg[5][klev][klon];
  double zmeltmax[klon];
  double zfrzmax[klon];
  double zicetot[klon];
  double zqxn2d[5][klev][klon];
  double zqsmix[klev][klon];
  double zqsliq[klev][klon];
  double zqsice[klev][klon];
  double zfoeewmt[klev][klon];
  double zfoeew[klev][klon];
  double zfoeeliqt[klev][klon];
  double zdqsliqdt[klon];
  double zdqsicedt[klon];
  double zdqsmixdt[klon];
  double zcorqsliq[klon];
  double zcorqsice[klon];
  double zcorqsmix[klon];
  double zevaplimliq[klon];
  double zevaplimice[klon];
  double zevaplimmix[klon];
  double zsolqa[5][5][klon];
  double zsolqb[5][5][klon];
  double zqlhs[5][5][klon];
  double zvqx[5];
  double zexplicit;
  double zratio[5][klon];
  double zsinksum[5][klon];
  double zfallsink[5][klon];
  double zfallsrce[5][klon];
  double zconvsrce[5][klon];
  double zconvsink[5][klon];
  double zpsupsatsrce[5][klon];
  double ztw1 = 1329.31000000000;
  double ztw2 = 0.00746150000000000;
  double ztw3 = 85000.0000000000;
  double ztw4 = 40.6370000000000;
  double ztw5 = 275.000000000000;
  double zsubsat;
  double ztdmtw0;
  double ztcg;
  double zfacx1i;
  double zfacx1s;
  double zaplusb;
  double zcorrfac;
  double zcorrfac2;
  double zpr02;
  double zterm1;
  double zterm2;
  double zcldtopdist[klon];
  double zinfactor;
  int iwarmrain;
  int ievaprain;
  int ievapsnow;
  int idepice;
  double zrainacc[klon];
  double zraincld[klon];
  double zsnowrime[klon];
  double zsnowcld[klon];
  double zesatliq;
  double zfallcorr;
  double zlambda;
  double zevap_denom;
  double zcorr2;
  double zka;
  double zconst;
  double ztemp;
  double zsumq0[klev][klon];
  double zsumq1[klev][klon];
  double zerrorq[klev][klon];
  double zsumh0[klev][klon];
  double zsumh1[klev][klon];
  double zerrorh[klev][klon];
  double zrain;
  double z_tmp1[kfdia-kidia+1];
  double z_tmp2[kfdia-kidia+1];
  double z_tmp3[kfdia-kidia+1];
  double z_tmp4[kfdia-kidia+1];
  double z_tmp6[kfdia-kidia+1];
  double z_tmp7[kfdia-kidia+1];
  double z_tmpk[klev][kfdia-kidia+1];
  double zhook_handle;
  double ztmpl;
  double ztmpi;
  double ztmpa;
  double zmm;
  double zrr;
  double zrg[klon];
  double zzsum;
  double zzratio;
  double zepsilon;
  double zcond1;
  double zqp;
  double psum_solqa[klon];
  int i_klev;
  int i_kfldx;
  int i_5;
  int i_klon;

  // REAL(KIND=JPRB) :: FOEDELTA
  // REAL(KIND=JPRB) :: PTARE
  // FOEDELTA (PTARE) = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-RTT))
  // REAL(KIND=JPRB) :: FOEALFA
  // FOEALFA (PTARE) = MIN(1.0_JPRB,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
  //  &*RTWAT_RTICE_R)**2) 
  // REAL(KIND=JPRB) :: FOEEWM,FOEDEM,FOELDCPM
  // FOEEWM ( PTARE ) = R2ES *&
  //      &(FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
  //   &(1.0_JPRB-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
  // FOEDEM ( PTARE ) = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
  //              &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
  // FOELDCPM ( PTARE ) = FOEALFA(PTARE)*RALVDCP+&
  //             &(1.0_JPRB-FOEALFA(PTARE))*RALSDCP
  // REAL(KIND=JPRB) :: FOEELIQ, FOEEICE 
  // FOEELIQ( PTARE ) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  // FOEEICE( PTARE ) = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))
  // #include "fccld.func.h"
  // REAL(KIND=JPRB) :: FOKOOP 
  // FOKOOP (PTARE) = MIN(RKOOP1-RKOOP2*PTARE,FOEELIQ(PTARE)/FOEEICE(PTARE))
  //===============================================================================
  //IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)
  //===============================================================================
  //  0.0     Beginning of timestep book-keeping
  //----------------------------------------------------------------------
  //######################################################################
  //             0.  *** SET UP CONSTANTS ***
  //######################################################################
  zepsilon = 100.0*DBL_EPSILON;    // --------------------------------------------------------------

  // ---------------------------------------------------------------------
  // Set version of warm-rain autoconversion/accretion
  // IWARMRAIN = 1 // Sundquist
  // IWARMRAIN = 2 // Khairoutdinov and Kogan (2000)
  // ---------------------------------------------------------------------
  iwarmrain = 2;    // ---------------------------------------------------------------------
  // Set version of rain evaporation
  // IEVAPRAIN = 1 // Sundquist
  // IEVAPRAIN = 2 // Abel and Boutle (2013)
  // ---------------------------------------------------------------------
  ievaprain = 2;    // ---------------------------------------------------------------------
  // Set version of snow evaporation
  // IEVAPSNOW = 1 // Sundquist
  // IEVAPSNOW = 2 // New
  // ---------------------------------------------------------------------
  ievapsnow = 1;    // ---------------------------------------------------------------------
  // Set version of ice deposition
  // IDEPICE = 1 // Rotstayn (2001)
  // IDEPICE = 2 // New
  // ---------------------------------------------------------------------
  idepice = 1;    // ---------------------
  // Some simple constants
  // ---------------------
  zqtmst = 1.0/ptsphy;
  zgdcp = rg/rcpd;
  zrdcp = rd/rcpd;
  zcons1a = rcpd/(((rlmlt*rg)*yrecldp->rtaumel));
  zepsec = 1.0e-14;
  zrg_r = 1.0/rg;
  zrldcp = 1.0/(ralsdcp - ralvdcp);    // Note: Defined in module/yoecldp.F90
  // NCLDQL=1    // liquid cloud water
  // NCLDQI=2    // ice cloud water
  // NCLDQR=3    // rain water
  // NCLDQS=4    // snow
  // NCLDQV=5    // vapour
  // -----------------------------------------------
  // Define species phase, 0=vapour, 1=liquid, 2=ice
  // -----------------------------------------------
  iphase[5-1] = 0;
  iphase[1-1] = 1;
  iphase[3-1] = 1;
  iphase[2-1] = 2;
  iphase[4-1] = 2;    // ---------------------------------------------------
  // Set up melting/freezing index, 
  // if an ice category melts/freezes, where does it go?
  // ---------------------------------------------------
  imelt[5-1] = -99;
  imelt[1-1] = 2;
  imelt[3-1] = 4;
  imelt[2-1] = 3;
  imelt[4-1] = 3;    // -----------------------------------------------
  // INITIALIZATION OF OUTPUT TENDENCIES
  // -----------------------------------------------
  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      tendency_loc_t[jk-1][jl-1] = 0.0;
      tendency_loc_q[jk-1][jl-1] = 0.0;
      tendency_loc_a[jk-1][jl-1] = 0.0;
    }

  }

  for (jm=1; jm<=4; jm+=1) {
    for (jk=1; jk<=klev; jk+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
        tendency_loc_cld[jm-1][jk-1][jl-1] = 0.0;
      }

    }

  }

  // -------------------------
  // set up fall speeds in m/s
  // -------------------------
  zvqx[5-1] = 0.0;
  zvqx[1-1] = 0.0;
  zvqx[2-1] = yrecldp->rvice;
  zvqx[3-1] = yrecldp->rvrain;
  zvqx[4-1] = yrecldp->rvsnow;
  for (i_5=1; i_5<=5; i_5++) {
    llfall[i_5-1] = false;
  }

  for (jm=1; jm<=5; jm+=1) {
    if (zvqx[jm-1] > 0.0)
    {
      llfall[jm-1] = true;
    }

  }

  // Set LLFALL to false for ice (but ice still sediments//)
  // Need to rationalise this at some point
  llfall[2-1] = false;    //######################################################################
  //             1.  *** INITIAL VALUES FOR VARIABLES ***
  //######################################################################
  // ----------------------
  // non CLV initialization 
  // ----------------------
  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      ztp1[jk-1][jl-1] = pt[jk-1][jl-1] + ptsphy*tendency_tmp_t[jk-1][jl-1];
      zqx[5-1][jk-1][jl-1] = pq[jk-1][jl-1] + ptsphy*tendency_tmp_q[jk-1][jl-1];
      zqx0[5-1][jk-1][jl-1] = pq[jk-1][jl-1] + ptsphy*tendency_tmp_q[jk-1][jl-1];
      za[jk-1][jl-1] = pa[jk-1][jl-1] + ptsphy*tendency_tmp_a[jk-1][jl-1];
      zaorig[jk-1][jl-1] = pa[jk-1][jl-1] + ptsphy*tendency_tmp_a[jk-1][jl-1];
    }

  }

  // -------------------------------------
  // initialization for CLV family
  // -------------------------------------
  for (jm=1; jm<=4; jm+=1) {
    for (jk=1; jk<=klev; jk+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
        zqx[jm-1][jk-1][jl-1] = pclv[jm-1][jk-1][jl-1] + ptsphy*tendency_tmp_cld[jm-1][jk-1][jl-1];
        zqx0[jm-1][jk-1][jl-1] = pclv[jm-1][jk-1][jl-1] + ptsphy*tendency_tmp_cld[jm-1][jk-1][jl-1];
      }

    }

  }

  //-------------
  // zero arrays
  //-------------
  for (jm=1; jm<=5; jm+=1) {
    for (jk=1; jk<=klev + 1; jk+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zpfplsx[jm-1][jk-1][jl-1] = 0.0;
      }

    }

  }

  for (jm=1; jm<=5; jm+=1) {
    for (jk=1; jk<=klev; jk+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zqxn2d[jm-1][jk-1][jl-1] = 0.0;
	zlneg[jm-1][jk-1][jl-1] = 0.0;
      }

    }

  }

  for (jl=kidia; jl<=kfdia; jl+=1) {
    prainfrac_toprfz[jl-1] = 0.0;
  }


  // ----------------------------------------------------
  // Tidy up very small cloud cover or total cloud water
  // ----------------------------------------------------
  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zqx[1-1][jk-1][jl-1] + zqx[2-1][jk-1][jl-1] < yrecldp->rlmin || za[jk-1][jl-1] < yrecldp->ramin)
      {
	// Evaporate small cloud liquid water amounts
	zlneg[1-1][jk-1][jl-1] = zlneg[1-1][jk-1][jl-1] + zqx[1-1][jk-1][jl-1];
	zqadj = zqx[1-1][jk-1][jl-1]*zqtmst;
	tendency_loc_q[jk-1][jl-1] = tendency_loc_q[jk-1][jl-1] + zqadj;
	tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] - ralvdcp*zqadj;
	zqx[5-1][jk-1][jl-1] = zqx[5-1][jk-1][jl-1] + zqx[1-1][jk-1][jl-1];
	zqx[1-1][jk-1][jl-1] = 0.0;            // Evaporate small cloud ice water amounts
	zlneg[2-1][jk-1][jl-1] = zlneg[2-1][jk-1][jl-1] + zqx[2-1][jk-1][jl-1];
	zqadj = zqx[2-1][jk-1][jl-1]*zqtmst;
	tendency_loc_q[jk-1][jl-1] = tendency_loc_q[jk-1][jl-1] + zqadj;
	tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] - ralsdcp*zqadj;
	zqx[5-1][jk-1][jl-1] = zqx[5-1][jk-1][jl-1] + zqx[2-1][jk-1][jl-1];
	zqx[2-1][jk-1][jl-1] = 0.0;            // Set cloud cover to zero
	za[jk-1][jl-1] = 0.0;
      }

    }

  }

  // ---------------------------------
  // Tidy up small CLV variables
  // ---------------------------------
  //DIR$ IVDEP
  for (jm=1; jm<=4; jm+=1) {
    //DIR$ IVDEP
    for (jk=1; jk<=klev; jk+=1) {
      //DIR$ IVDEP
      for (jl=kidia; jl<=kfdia; jl+=1) {
	if (zqx[jm-1][jk-1][jl-1] < yrecldp->rlmin)
	{
	  zlneg[jm-1][jk-1][jl-1] = zlneg[jm-1][jk-1][jl-1] + zqx[jm-1][jk-1][jl-1];
	  zqadj = zqx[jm-1][jk-1][jl-1]*zqtmst;
	  tendency_loc_q[jk-1][jl-1] = tendency_loc_q[jk-1][jl-1] + zqadj;
	  if (iphase[jm-1] == 1)
	  {
	    tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] - ralvdcp*zqadj;
	  }

	  if (iphase[jm-1] == 2)
	  {
	    tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] - ralsdcp*zqadj;
	  }

	  zqx[5-1][jk-1][jl-1] = zqx[5-1][jk-1][jl-1] + zqx[jm-1][jk-1][jl-1];
	  zqx[jm-1][jk-1][jl-1] = 0.0;
	}

      }

    }

  }

  // ------------------------------
  // Define saturation values
  // ------------------------------
  // IF(N_VMASS > 1) THEN
  //   JLEN=KFDIA-KIDIA+1
  // ENDIF
  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      //----------------------------------------
      // old *diagnostic* mixed phase saturation
      //---------------------------------------- 
      zfoealfa[jk-1][jl-1] = (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)));
      zfoeewmt[jk-1][jl-1] = fmin((double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies))))/pap[jk-1][jl-1], 0.5);
      zqsmix[jk-1][jl-1] = zfoeewmt[jk-1][jl-1];
      zqsmix[jk-1][jl-1] = zqsmix[jk-1][jl-1]/(1.0 - retv*zqsmix[jk-1][jl-1]);        //---------------------------------------------
      // ice saturation T<273K
      // liquid water saturation for T>273K 
      //---------------------------------------------
      zalfa = (double)(fmax(0.0, copysign(1.0, ztp1[jk-1][jl-1] - rtt)));
      zfoeew[jk-1][jl-1] = fmin((zalfa*(double)(r2es*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les))) + (1.0 - zalfa)*(double)(r2es*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies))))/pap[jk-1][jl-1], 0.5);
      zfoeew[jk-1][jl-1] = fmin(0.5, zfoeew[jk-1][jl-1]);
      zqsice[jk-1][jl-1] = zfoeew[jk-1][jl-1]/(1.0 - retv*zfoeew[jk-1][jl-1]);        //----------------------------------
      // liquid water saturation
      //---------------------------------- 
      zfoeeliqt[jk-1][jl-1] = fmin((double)(r2es*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)))/pap[jk-1][jl-1], 0.5);
      zqsliq[jk-1][jl-1] = zfoeeliqt[jk-1][jl-1];
      zqsliq[jk-1][jl-1] = zqsliq[jk-1][jl-1]/(1.0 - retv*zqsliq[jk-1][jl-1]);        //      //----------------------------------
      //      // ice water saturation
      //      //---------------------------------- 
      //      ZFOEEICET(JL,JK)=MIN(FOEEICE(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      //      ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
      //      ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-RETV*ZQSICE(JL,JK))
    }

    // Calculation of relative humidity for diagnostics
    // Not used at present
    //  DO JL=KIDIA,KFDIA
    //    ZRHM(JL,JK)=ZQX(JL,JK,NCLDQV)/ZQSMIX(JL,JK)
    //    ZRHL(JL,JK)=ZQX(JL,JK,NCLDQV)/ZQSLIQ(JL,JK)
    //    ZRHI(JL,JK)=ZQX(JL,JK,NCLDQV)/ZQSICE(JL,JK)
    //  ENDDO
  }

  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      //------------------------------------------
      // Ensure cloud fraction is between 0 and 1
      //------------------------------------------
      za[jk-1][jl-1] = fmax(0.0, fmin(1.0, za[jk-1][jl-1]));        //-------------------------------------------------------------------
      // Calculate liq/ice fractions (no longer a diagnostic relationship)
      //-------------------------------------------------------------------
      zli[jk-1][jl-1] = zqx[1-1][jk-1][jl-1] + zqx[2-1][jk-1][jl-1];
      if (zli[jk-1][jl-1] > yrecldp->rlmin)
      {
        zliqfrac[jk-1][jl-1] = zqx[1-1][jk-1][jl-1]/zli[jk-1][jl-1];
        zicefrac[jk-1][jl-1] = 1.0 - zliqfrac[jk-1][jl-1];
      } else {
        zliqfrac[jk-1][jl-1] = 0.0;
        zicefrac[jk-1][jl-1] = 0.0;
      }

    }

  }

  //----------------------------------------------------------------------
  // This is test code. Instead of resetting cloud water or cloud cover
  // to zero if one of the two is zero, here we slave the cloud cover
  // to the cloud water variable. I.e. if cloud cover is zero it is 
  // set to an appropriate non-zero value.
  // It uses a Beta curve to get the variance and then derive cloud cover.
  // It is quite slow since it involves iteration, and should be left 
  // until a fully prognostic variance equation for total water is 
  // implemented
  //-----------------------------------------------------------------------
  //IF (.FALSE.) THEN
  //  CALL CLOUDVAR &
  //---input
  // & ( KIDIA, KFDIA, KLON  , KLEV  , 1    , KLEV, &
  // &   ZTP1, ZQP1, ZQSMIX, ZLI, PAP, &
  //---output
  // &   ZVAR, ZQTMIN, ZQTMAX ) //last two are dummy args  
  //  CALL COVER &
  //---input
  // & ( KIDIA, KFDIA , KLON, KLEV, 1, KLEV, &
  // &   ZA, ZTP1,  ZQP1, ZQSMIX, ZLI, PAP, ZVAR, &
  //---output
  // &   ZQTMAX, ZABETA )  
  //  DO JK=1,KLEV
  //    DO JL=KIDIA,KFDIA
  //      IF (ZLI(JL,JK)/MAX(ZA(JL,JK),ZEPSEC)>RCLDMAX) THEN
  //        ZA(JL,JK)=ZABETA(JL,JK) // not part of tendency       
  //      ENDIF
  //    ENDDO
  //  ENDDO
  //ENDIF
  //--------------------------------------
  // NPM
  // Initialize liq water temperature T_L 
  // Not used at present
  //--------------------------------------
  //ZTL(:,:)=ZTP1(:,:)
  //DO JM=1,NCLV
  //  DO JK=1,KLEV
  //    DO JL=KIDIA,KFDIA
  //      IF (IPHASE(JM)==1) ZTL(JL,JK)=ZTL(JL,JK)-RALVDCP*ZQX(JL,JK,JM)
  //      IF (IPHASE(JM)==2) ZTL(JL,JK)=ZTL(JL,JK)-RALSDCP*ZQX(JL,JK,JM)
  //    ENDDO
  //  ENDDO
  //ENDDO
  //######################################################################
  //        2.       *** CONSTANTS AND PARAMETERS ***
  //######################################################################
  //  Calculate L in updrafts of bl-clouds
  //  Specify QS, P/PS for tropopause (for c2)
  //  And initialize variables
  //------------------------------------------
  //---------------------------------
  // Find tropopause level (ZTRPAUS)
  //---------------------------------
  for (jl=kidia; jl<=kfdia; jl+=1) {
    ztrpaus[jl-1] = 0.1;
    zpaphd[jl-1] = 1.0/paph[klev+1-1][jl-1];
  }

  for (jk=1; jk<=klev - 1; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zsig = pap[jk-1][jl-1]*zpaphd[jl-1];
      if (zsig > 0.1 && ztp1[jk-1][jl-1] > ztp1[jk+1-1][jl-1] && zsig < 0.4)
      {
        ztrpaus[jl-1] = zsig;
      }

    }

  }

  //-----------------------------
  // Reset single level variables
  //-----------------------------
  for (jl=kidia; jl<=kfdia; jl+=1) {
    zanewm1[jl-1] = 0.0;
    zda[jl-1] = 0.0;
    zcovpclr[jl-1] = 0.0;
    zcovpmax[jl-1] = 0.0;
    zcovptot[jl-1] = 0.0;
    zcldtopdist[jl-1] = 0.0;
  }

  //######################################################################
  //           3.       *** PHYSICS ***
  //######################################################################
  //----------------------------------------------------------------------
  //                       START OF VERTICAL LOOP
  //----------------------------------------------------------------------
  for (jk=yrecldp->ncldtop; jk<=klev; jk+=1) {
    //----------------------------------------------------------------------
    // 3.0 INITIALIZE VARIABLES
    //----------------------------------------------------------------------
    //---------------------------------
    // First guess microphysics
    //---------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
        zqxfg[jm-1][jl-1] = zqx[jm-1][jk-1][jl-1];
      }

    }

    //---------------------------------
    // Set KLON arrays to zero
    //---------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zlicld[jl-1] = 0.0;
      zrainaut[jl-1] = 0.0;
      zrainacc[jl-1] = 0.0;
      zsnowaut[jl-1] = 0.0;
      zldefr[jl-1] = 0.0;
      zacust[jl-1] = 0.0;
      zqpretot[jl-1] = 0.0;
      zlfinalsum[jl-1] = 0.0;        // Required for first guess call
      zlcond1[jl-1] = 0.0;
      zlcond2[jl-1] = 0.0;
      zsupsat[jl-1] = 0.0;
      zlevapl[jl-1] = 0.0;
      zlevapi[jl-1] = 0.0;        //-------------------------------------                
      // solvers for cloud fraction                          
      //-------------------------------------                
      zsolab[jl-1] = 0.0;
      zsolac[jl-1] = 0.0;
      zicetot[jl-1] = 0.0;
    }

    //------------------------------------------           
    // reset matrix so missing pathways are set            
    //------------------------------------------           
    for (jm=1; jm<=5; jm+=1) {
      for (jn=1; jn<=5; jn+=1) {
        for (jl=kidia; jl<=kfdia; jl+=1) {
          zsolqb[jm-1][jn-1][jl-1] = 0.0;
          zsolqa[jm-1][jn-1][jl-1] = 0.0;
        }

      }

    }

    //----------------------------------                   
    // reset new microphysics variables                    
    //----------------------------------                   
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
        zfallsrce[jm-1][jl-1] = 0.0;
        zfallsink[jm-1][jl-1] = 0.0;
        zconvsrce[jm-1][jl-1] = 0.0;
        zconvsink[jm-1][jl-1] = 0.0;
        zpsupsatsrce[jm-1][jl-1] = 0.0;
        zratio[jm-1][jl-1] = 0.0;
      }

    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      //-------------------------
      // derived variables needed
      //-------------------------
      zdp[jl-1] = paph[jk+1-1][jl-1] - paph[jk-1][jl-1];
      zgdp[jl-1] = rg/zdp[jl-1];
      zrho[jl-1] = pap[jk-1][jl-1]/((rd*ztp1[jk-1][jl-1]));
      zdtgdp[jl-1] = ptsphy*zgdp[jl-1];
      zrdtgdp[jl-1] = zdp[jl-1]*(1.0/((ptsphy*rg)));
      if (jk > 1)
      {
        zdtgdpf[jl-1] = (ptsphy*rg)/(pap[jk-1][jl-1] - pap[jk-1-1][jl-1]);
      }

      //------------------------------------
      // Calculate dqs/dT correction factor
      //------------------------------------
      // Reminder: RETV=RV/RD-1
      // liquid
      zfacw = r5les/pow(ztp1[jk-1][jl-1] - r4les, 2);
      zcor = 1.0/(1.0 - retv*zfoeeliqt[jk-1][jl-1]);
      zdqsliqdt[jl-1] = (zfacw*zcor)*zqsliq[jk-1][jl-1];
      zcorqsliq[jl-1] = 1.0 + ralvdcp*zdqsliqdt[jl-1];        // ice
      zfaci = r5ies/pow(ztp1[jk-1][jl-1] - r4ies, 2);
      zcor = 1.0/(1.0 - retv*zfoeew[jk-1][jl-1]);
      zdqsicedt[jl-1] = (zfaci*zcor)*zqsice[jk-1][jl-1];
      zcorqsice[jl-1] = 1.0 + ralsdcp*zdqsicedt[jl-1];        // diagnostic mixed
      zalfaw = zfoealfa[jk-1][jl-1];
      zalfawm[jl-1] = zalfaw;
      zfac = zalfaw*zfacw + (1.0 - zalfaw)*zfaci;
      zcor = 1.0/(1.0 - retv*zfoeewmt[jk-1][jl-1]);
      zdqsmixdt[jl-1] = (zfac*zcor)*zqsmix[jk-1][jl-1];
      zcorqsmix[jl-1] = 1.0 + (double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)*zdqsmixdt[jl-1];        // evaporation/sublimation limits
      zevaplimmix[jl-1] = fmax((zqsmix[jk-1][jl-1] - zqx[5-1][jk-1][jl-1])/zcorqsmix[jl-1], 0.0);
      zevaplimliq[jl-1] = fmax((zqsliq[jk-1][jl-1] - zqx[5-1][jk-1][jl-1])/zcorqsliq[jl-1], 0.0);
      zevaplimice[jl-1] = fmax((zqsice[jk-1][jl-1] - zqx[5-1][jk-1][jl-1])/zcorqsice[jl-1], 0.0);        //--------------------------------
      // in-cloud consensate amount
      //--------------------------------
      ztmpa = 1.0*1.0/fmax(za[jk-1][jl-1], zepsec);
      zliqcld[jl-1] = zqx[1-1][jk-1][jl-1]*ztmpa;
      zicecld[jl-1] = zqx[2-1][jk-1][jl-1]*ztmpa;
      zlicld[jl-1] = zliqcld[jl-1] + zicecld[jl-1];
    }

    //------------------------------------------------
    // Evaporate very small amounts of liquid and ice
    //------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zqx[1-1][jk-1][jl-1] < yrecldp->rlmin)
      {
        zsolqa[1-1][5-1][jl-1] = zqx[1-1][jk-1][jl-1];
        zsolqa[5-1][1-1][jl-1] = -zqx[1-1][jk-1][jl-1];
      }

      if (zqx[2-1][jk-1][jl-1] < yrecldp->rlmin)
      {
        zsolqa[2-1][5-1][jl-1] = zqx[2-1][jk-1][jl-1];
        zsolqa[5-1][2-1][jl-1] = -zqx[2-1][jk-1][jl-1];
      }

    }

    //---------------------------------------------------------------------
    //  3.1  ICE SUPERSATURATION ADJUSTMENT
    //---------------------------------------------------------------------
    // Note that the supersaturation adjustment is made with respect to 
    // liquid saturation:  when T>0C 
    // ice saturation:     when T<0C
    //                     with an adjustment made to allow for ice 
    //                     supersaturation in the clear sky
    // Note also that the KOOP factor automatically clips the supersaturation
    // to a maximum set by the liquid water saturation mixing ratio
    // important for temperatures near to but below 0C
    //----------------------------------------------------------------------- 
    //DIR$ NOFUSION
    for (jl=kidia; jl<=kfdia; jl+=1) {
      //-----------------------------------
      // 3.1.1 Supersaturation limit (from Koop)
      //-----------------------------------
      // Needs to be set for all temperatures
      zfokoop[jl-1] = (double)(fmin(rkoop1 - rkoop2*ztp1[jk-1][jl-1], (double)(r2es*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)))*1.0/(double)(r2es*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies)))));
    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (yrecldp->nssopt == 0 || ztp1[jk-1][jl-1] >= rtt)
      {
        zfac = 1.0;
        zfaci = 1.0;
      } else {
        zfac = za[jk-1][jl-1] + zfokoop[jl-1]*(1.0 - za[jk-1][jl-1]);
        zfaci = ptsphy/yrecldp->rkooptau;
      }

      //-------------------------------------------------------------------
      // 3.1.2 Calculate supersaturation wrt Koop including dqs/dT 
      //       correction factor
      // [#Note: QSICE or QSLIQ]
      //-------------------------------------------------------------------
      // Calculate supersaturation to add to cloud
      if (za[jk-1][jl-1] > 1.0 - yrecldp->ramin)
      {
        zsupsat[jl-1] = fmax((zqx[5-1][jk-1][jl-1] - zfac*zqsice[jk-1][jl-1])/zcorqsice[jl-1], 0.0);
      } else {
        // Calculate environmental humidity supersaturation
        zqp1env = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsice[jk-1][jl-1])*1.0/fmax(1.0 - za[jk-1][jl-1], zepsilon);          //& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
        zsupsat[jl-1] = fmax(((1.0 - za[jk-1][jl-1])*(zqp1env - zfac*zqsice[jk-1][jl-1]))/zcorqsice[jl-1], 0.0);
      }

      //-------------------------------------------------------------------
      // Here the supersaturation is turned into liquid water
      // However, if the temperature is below the threshold for homogeneous
      // freezing then the supersaturation is turned instantly to ice.
      //--------------------------------------------------------------------
      if (zsupsat[jl-1] > zepsec)
      {
        if (ztp1[jk-1][jl-1] > yrecldp->rthomo)
        {
          // Turn supersaturation into liquid water        
          zsolqa[5-1][1-1][jl-1] = zsolqa[5-1][1-1][jl-1] + zsupsat[jl-1];
          zsolqa[1-1][5-1][jl-1] = zsolqa[1-1][5-1][jl-1] - zsupsat[jl-1];            // Include liquid in first guess
          zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] + zsupsat[jl-1];
        } else {
          // Turn supersaturation into ice water        
          zsolqa[5-1][2-1][jl-1] = zsolqa[5-1][2-1][jl-1] + zsupsat[jl-1];
          zsolqa[2-1][5-1][jl-1] = zsolqa[2-1][5-1][jl-1] - zsupsat[jl-1];            // Add ice to first guess for deposition term 
          zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + zsupsat[jl-1];
        }

        // Increase cloud amount using RKOOPTAU timescale
        zsolac[jl-1] = (1.0 - za[jk-1][jl-1])*zfaci;          // Store cloud budget diagnostics if required
      }

      //-------------------------------------------------------
      // 3.1.3 Include supersaturation from previous timestep
      // (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
      //-------------------------------------------------------
      if (psupsat[jk-1][jl-1] > zepsec)
      {
	if (ztp1[jk-1][jl-1] > yrecldp->rthomo)
	{
	  // Turn supersaturation into liquid water
	  zsolqa[1-1][1-1][jl-1] = zsolqa[1-1][1-1][jl-1] + psupsat[jk-1][jl-1];
	  zpsupsatsrce[1-1][jl-1] = psupsat[jk-1][jl-1];                // Add liquid to first guess for deposition term 
	  zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] + psupsat[jk-1][jl-1];                // Store cloud budget diagnostics if required
	} else {
	  // Turn supersaturation into ice water
	  zsolqa[2-1][2-1][jl-1] = zsolqa[2-1][2-1][jl-1] + psupsat[jk-1][jl-1];
	  zpsupsatsrce[2-1][jl-1] = psupsat[jk-1][jl-1];                // Add ice to first guess for deposition term 
	  zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + psupsat[jk-1][jl-1];                // Store cloud budget diagnostics if required
	}

	// Increase cloud amount using RKOOPTAU timescale
	zsolac[jl-1] = (1.0 - za[jk-1][jl-1])*zfaci;              // Store cloud budget diagnostics if required
      }

    }

    //---------------------------------------------------------------------
    //  3.2  DETRAINMENT FROM CONVECTION
    //---------------------------------------------------------------------
    // * Diagnostic T-ice/liq split retained for convection
    //    Note: This link is now flexible and a future convection 
    //    scheme can detrain explicit seperate budgets of:
    //    cloud water, ice, rain and snow
    // * There is no (1-ZA) multiplier term on the cloud detrainment 
    //    term, since is now written in mass-flux terms  
    // [#Note: Should use ZFOEALFACU used in convection rather than ZFOEALFA]
    //---------------------------------------------------------------------
    if (jk >= yrecldp->ncldtop && jk < klev)
    {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	plude[jk-1][jl-1] = plude[jk-1][jl-1]*zdtgdp[jl-1];
	if (plu[jk+1-1][jl-1] > zepsec && plude[jk-1][jl-1] > yrecldp->rlmin)
	{
	  zsolac[jl-1] = zsolac[jl-1] + plude[jk-1][jl-1]/plu[jk+1-1][jl-1];              // *diagnostic temperature split*
	  zalfaw = zfoealfa[jk-1][jl-1];
	  zconvsrce[1-1][jl-1] = zalfaw*plude[jk-1][jl-1];
	  zconvsrce[2-1][jl-1] = (1.0 - zalfaw)*plude[jk-1][jl-1];
	  zsolqa[1-1][1-1][jl-1] = zsolqa[1-1][1-1][jl-1] + zconvsrce[1-1][jl-1];
	  zsolqa[2-1][2-1][jl-1] = zsolqa[2-1][2-1][jl-1] + zconvsrce[2-1][jl-1];              // Store cloud budget diagnostics if required
	} else {
	  plude[jk-1][jl-1] = 0.0;
	}

	// *convective snow detrainment source
	zsolqa[4-1][4-1][jl-1] = zsolqa[4-1][4-1][jl-1] + psnde[jk-1][jl-1]*zdtgdp[jl-1];
      }

    }

    //---------------------------------------------------------------------
    //  3.3  SUBSIDENCE COMPENSATING CONVECTIVE UPDRAUGHTS
    //---------------------------------------------------------------------
    // Three terms:
    // * Convective subsidence source of cloud from layer above
    // * Evaporation of cloud within the layer
    // * Subsidence sink of cloud to the layer below (Implicit solution)
    //---------------------------------------------------------------------
    //-----------------------------------------------
    // Subsidence source from layer above
    //               and 
    // Evaporation of cloud within the layer
    //-----------------------------------------------
    if (jk > yrecldp->ncldtop)
    {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zmf[jl-1] = fmax(0.0, (pmfu[jk-1][jl-1] + pmfd[jk-1][jl-1])*zdtgdp[jl-1]);
	zacust[jl-1] = zmf[jl-1]*zanewm1[jl-1];
      }

      for (jm=1; jm<=5; jm+=1) {
	if (!llfall[jm-1] && iphase[jm-1] > 0)
	{
	  for (jl=kidia; jl<=kfdia; jl+=1) {
	    zlcust[jm-1][jl-1] = zmf[jl-1]*zqxnm1[jm-1][jl-1];                // record total flux for enthalpy budget:
	    zconvsrce[jm-1][jl-1] = zconvsrce[jm-1][jl-1] + zlcust[jm-1][jl-1];
	  }

	}

      }

      // Now have to work out how much liquid evaporates at arrival point 
      // since there is no prognostic memory for in-cloud humidity, i.e. 
      // we always assume cloud is saturated. 
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zdtdp = ((zrdcp*0.5)*(ztp1[jk-1-1][jl-1] + ztp1[jk-1][jl-1]))/paph[jk-1][jl-1];
	zdtforc = zdtdp*(pap[jk-1][jl-1] - pap[jk-1-1][jl-1]);            //[#Note: Diagnostic mixed phase should be replaced below]
	zdqs[jl-1] = (zanewm1[jl-1]*zdtforc)*zdqsmixdt[jl-1];
      }

      for (jm=1; jm<=5; jm+=1) {
	if (!llfall[jm-1] && iphase[jm-1] > 0)
	{
	  for (jl=kidia; jl<=kfdia; jl+=1) {
	    zlfinal = fmax(0.0, zlcust[jm-1][jl-1] - zdqs[jl-1]);                // no supersaturation allowed incloud ---V
	    zevap = fmin(zlcust[jm-1][jl-1] - zlfinal, zevaplimmix[jl-1]);                //          ZEVAP=0.0_JPRB
	    zlfinal = zlcust[jm-1][jl-1] - zevap;
	    zlfinalsum[jl-1] = zlfinalsum[jl-1] + zlfinal;
	    zsolqa[jm-1][jm-1][jl-1] = zsolqa[jm-1][jm-1][jl-1] + zlcust[jm-1][jl-1];
	    zsolqa[jm-1][5-1][jl-1] = zsolqa[jm-1][5-1][jl-1] + zevap;
	    zsolqa[5-1][jm-1][jl-1] = zsolqa[5-1][jm-1][jl-1] - zevap;                // Store cloud liquid diagnostic if required
	  }

	}

      }

      //  Reset the cloud contribution if no cloud water survives to this level:
      for (jl=kidia; jl<=kfdia; jl+=1) {
	if (zlfinalsum[jl-1] < zepsec)
	{
	  zacust[jl-1] = 0.0;
	}

	zsolac[jl-1] = zsolac[jl-1] + zacust[jl-1];            // Store cloud fraction diagnostic if required
	//ZBUDCC(JL,5) isn't included as only reduced if cloud->zero
      }

    }

    //---------------------------------------------------------------------
    // Subsidence sink of cloud to the layer below 
    // (Implicit - re. CFL limit on convective mass flux)
    //---------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (jk < klev)
      {
	zmfdn = fmax(0.0, (pmfu[jk+1-1][jl-1] + pmfd[jk+1-1][jl-1])*zdtgdp[jl-1]);
	zsolab[jl-1] = zsolab[jl-1] + zmfdn;
	zsolqb[1-1][1-1][jl-1] = zsolqb[1-1][1-1][jl-1] + zmfdn;
	zsolqb[2-1][2-1][jl-1] = zsolqb[2-1][2-1][jl-1] + zmfdn;            // Record sink for cloud budget and enthalpy budget diagnostics
	zconvsink[1-1][jl-1] = zmfdn;
	zconvsink[2-1][jl-1] = zmfdn;
      }

    }

    //----------------------------------------------------------------------
    // 3.4  EROSION OF CLOUDS BY TURBULENT MIXING
    //----------------------------------------------------------------------
    // NOTE: In default tiedtke scheme this process decreases the cloud 
    //       area but leaves the specific cloud water content
    //       within clouds unchanged
    //----------------------------------------------------------------------
    // ------------------------------
    // Define turbulent erosion rate
    // ------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zldifdt[jl-1] = yrecldp->rcldiff*ptsphy;          //Increase by factor of 5 for convective points
      if (ktype[jl-1] > 0 && plude[jk-1][jl-1] > zepsec)
      {
	zldifdt[jl-1] = yrecldp->rcldiff_convi*zldifdt[jl-1];
      }

    }

    // At the moment, works on mixed RH profile and partitioned ice/liq fraction
    // so that it is similar to previous scheme
    // Should apply RHw for liquid cloud and RHi for ice cloud separately 
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zli[jk-1][jl-1] > zepsec)
      {
	// Calculate environmental humidity
	//      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
	//    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
	//      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
	ze = zldifdt[jl-1]*fmax(zqsmix[jk-1][jl-1] - zqx[5-1][jk-1][jl-1], 0.0);
	zleros = za[jk-1][jl-1]*ze;
	zleros = fmin(zleros, zevaplimmix[jl-1]);
	zleros = fmin(zleros, zli[jk-1][jl-1]);
	zaeros = zleros/zlicld[jl-1];            // Erosion is -ve LINEAR in L,A
	zsolac[jl-1] = zsolac[jl-1] - zaeros;
	zsolqa[1-1][5-1][jl-1] = zsolqa[1-1][5-1][jl-1] + zliqfrac[jk-1][jl-1]*zleros;
	zsolqa[5-1][1-1][jl-1] = zsolqa[5-1][1-1][jl-1] - zliqfrac[jk-1][jl-1]*zleros;
	zsolqa[2-1][5-1][jl-1] = zsolqa[2-1][5-1][jl-1] + zicefrac[jk-1][jl-1]*zleros;
	zsolqa[5-1][2-1][jl-1] = zsolqa[5-1][2-1][jl-1] - zicefrac[jk-1][jl-1]*zleros;            // Store cloud budget diagnostics if required
      }

    }

    //----------------------------------------------------------------------
    // 3.4  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
    //----------------------------------------------------------------------
    //  calculate dqs/dt
    //  Note: For the separate prognostic Qi and Ql, one would ideally use
    //  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
    //  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
    //  These would then instantaneous freeze if T<-38C or lead to ice growth 
    //  by deposition in warmer mixed phase clouds.  However, since we do 
    //  not have a separate prognostic equation for in-cloud humidity or a 
    //  statistical scheme approach in place, the depositional growth of ice 
    //  in the mixed phase can not be modelled and we resort to supersaturation  
    //  wrt ice instanteously converting to ice over one timestep 
    //  (see Tompkins et al. QJRMS 2007 for details)
    //  Thus for the initial implementation the diagnostic mixed phase is 
    //  retained for the moment, and the level of approximation noted.  
    //----------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zdtdp = (zrdcp*ztp1[jk-1][jl-1])/pap[jk-1][jl-1];
      zdpmxdt = zdp[jl-1]*zqtmst;
      zmfdn = 0.0;
      if (jk < klev)
      {
        zmfdn = pmfu[jk+1-1][jl-1] + pmfd[jk+1-1][jl-1];
      }

      zwtot = pvervel[jk-1][jl-1] + (0.5*rg)*(pmfu[jk-1][jl-1] + pmfd[jk-1][jl-1] + zmfdn);
      zwtot = fmin(zdpmxdt, fmax(-zdpmxdt, zwtot));
      zzzdt = phrsw[jk-1][jl-1] + phrlw[jk-1][jl-1];
      zdtdiab = fmin(zdpmxdt*zdtdp, fmax(-zdpmxdt*zdtdp, zzzdt))*ptsphy + ralfdcp*zldefr[jl-1];        // Note: ZLDEFR should be set to the difference between the mixed phase functions
      // in the convection and cloud scheme, but this is not calculated, so is zero and
      // the functions must be the same
      zdtforc = (zdtdp*zwtot)*ptsphy + zdtdiab;
      zqold[jl-1] = zqsmix[jk-1][jl-1];
      ztold[jl-1] = ztp1[jk-1][jl-1];
      ztp1[jk-1][jl-1] = ztp1[jk-1][jl-1] + zdtforc;
      ztp1[jk-1][jl-1] = fmax(ztp1[jk-1][jl-1], 160.0);
      llflag[jl-1] = true;
    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      zqp = 1.0/pap[jk-1][jl-1];
      zqsat = (double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies))))*zqp;
      zqsat = fmin(0.5, zqsat);
      zcor = 1.0/(1.0 - retv*zqsat);
      zqsat = zqsat*zcor;
      zcond = (zqsmix[jk-1][jl-1] - zqsat)*1.0/(1.0 + (zqsat*zcor)*(double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1][jl-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1][jl-1] - r4ies, 2))));
      ztp1[jk-1][jl-1] = ztp1[jk-1][jl-1] + (double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)*zcond;
      zqsmix[jk-1][jl-1] = zqsmix[jk-1][jl-1] - zcond;
      zqsat = (double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies))))*zqp;
      zqsat = fmin(0.5, zqsat);
      zcor = 1.0/(1.0 - retv*zqsat);
      zqsat = zqsat*zcor;
      zcond1 = (zqsmix[jk-1][jl-1] - zqsat)*1.0/(1.0 + (zqsat*zcor)*(double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1][jl-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1][jl-1] - r4ies, 2))));
      ztp1[jk-1][jl-1] = ztp1[jk-1][jl-1] + (double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)*zcond1;
      zqsmix[jk-1][jl-1] = zqsmix[jk-1][jl-1] - zcond1;
    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      zdqs[jl-1] = zqsmix[jk-1][jl-1] - zqold[jl-1];
      zqsmix[jk-1][jl-1] = zqold[jl-1];
      ztp1[jk-1][jl-1] = ztold[jl-1];
    }

    //----------------------------------------------------------------------
    // 3.4a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
    // ----------------------------------------------------------------------
    // Erosion term is LINEAR in L
    // Changed to be uniform distribution in cloud region
    for (jl=kidia; jl<=kfdia; jl+=1) {
      // Previous function based on DELTA DISTRIBUTION in cloud:
      if (zdqs[jl-1] > 0.0)
      {
        //    If subsidence evaporation term is turned off, then need to use updated
        //    liquid and cloud here?
        //    ZLEVAP = MAX(ZA(JL,JK)+ZACUST(JL),1.0_JPRB)*MIN(ZDQS(JL),ZLICLD(JL)+ZLFINALSUM(JL))
        zlevap = za[jk-1][jl-1]*fmin(zdqs[jl-1], zlicld[jl-1]);
        zlevap = fmin(zlevap, zevaplimmix[jl-1]);
        zlevap = fmin(zlevap, fmax(zqsmix[jk-1][jl-1] - zqx[5-1][jk-1][jl-1], 0.0));          // For first guess call
        zlevapl[jl-1] = zliqfrac[jk-1][jl-1]*zlevap;
        zlevapi[jl-1] = zicefrac[jk-1][jl-1]*zlevap;
        zsolqa[1-1][5-1][jl-1] = zsolqa[1-1][5-1][jl-1] + zliqfrac[jk-1][jl-1]*zlevap;
        zsolqa[5-1][1-1][jl-1] = zsolqa[5-1][1-1][jl-1] - zliqfrac[jk-1][jl-1]*zlevap;
        zsolqa[2-1][5-1][jl-1] = zsolqa[2-1][5-1][jl-1] + zicefrac[jk-1][jl-1]*zlevap;
        zsolqa[5-1][2-1][jl-1] = zsolqa[5-1][2-1][jl-1] - zicefrac[jk-1][jl-1]*zlevap;          // Store cloud budget diagnostics if required
      }

    }

    //----------------------------------------------------------------------
    // 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
    //----------------------------------------------------------------------
    // (1) Increase of cloud water in existing clouds
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zdqs[jl-1] <= -yrecldp->rlmin && za[jk-1][jl-1] > zepsec)
      {
        zlcond1[jl-1] = fmax(-zdqs[jl-1], 0.0);          //old limiter (significantly improves upper tropospheric humidity rms)
        if (za[jk-1][jl-1] > 0.99)
        {
          zcor = 1.0/(1.0 - retv*zqsmix[jk-1][jl-1]);
          zcdmax = (zqx[5-1][jk-1][jl-1] - zqsmix[jk-1][jl-1])*1.0/(1.0 + (zcor*zqsmix[jk-1][jl-1])*(double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1][jl-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1][jl-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1][jl-1] - r4ies, 2))));
        } else {
          zcdmax = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsmix[jk-1][jl-1])/za[jk-1][jl-1];
        }

        zlcond1[jl-1] = fmax(fmin(zlcond1[jl-1], zcdmax), 0.0);          // end old limiter
        zlcond1[jl-1] = za[jk-1][jl-1]*zlcond1[jl-1];
        if (zlcond1[jl-1] < yrecldp->rlmin)
        {
          zlcond1[jl-1] = 0.0;
        }

        //-------------------------------------------------------------------------
        // All increase goes into liquid unless so cold cloud homogeneously freezes
        // Include new liquid formation in first guess value, otherwise liquid 
        // remains at cold temperatures until next timestep.
        //-------------------------------------------------------------------------
        if (ztp1[jk-1][jl-1] > yrecldp->rthomo)
        {
          zsolqa[5-1][1-1][jl-1] = zsolqa[5-1][1-1][jl-1] + zlcond1[jl-1];
          zsolqa[1-1][5-1][jl-1] = zsolqa[1-1][5-1][jl-1] - zlcond1[jl-1];
          zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] + zlcond1[jl-1];            // Store cloud liquid diagnostic if required
        } else {
          zsolqa[5-1][2-1][jl-1] = zsolqa[5-1][2-1][jl-1] + zlcond1[jl-1];
          zsolqa[2-1][5-1][jl-1] = zsolqa[2-1][5-1][jl-1] - zlcond1[jl-1];
          zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + zlcond1[jl-1];            // Store cloud ice diagnostic if required
        }

      }

    }

    // (2) Generation of new clouds (da/dt>0)
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zdqs[jl-1] <= -yrecldp->rlmin && za[jk-1][jl-1] < 1.0 - zepsec)
      {
        //---------------------------
        // Critical relative humidity
        //---------------------------
        zrhc = yrecldp->ramid;
        zsigk = pap[jk-1][jl-1]/paph[klev+1-1][jl-1];          // Increase RHcrit to 1.0 towards the surface (eta>0.8)
        if (zsigk > 0.8)
        {
          zrhc = yrecldp->ramid + (1.0 - yrecldp->ramid)*pow((zsigk - 0.8)/0.2, 2);
        }

        // Commented out for CY37R1 to reduce humidity in high trop and strat
        //      // Increase RHcrit to 1.0 towards the tropopause (trop-0.2) and above
        //      ZBOTT=ZTRPAUS(JL)+0.2_JPRB
        //      IF(ZSIGK < ZBOTT) THEN
        //        ZRHC=RAMID+(1.0_JPRB-RAMID)*MIN(((ZBOTT-ZSIGK)/0.2_JPRB)**2,1.0_JPRB)
        //      ENDIF
        //---------------------------
        // Supersaturation options
        //---------------------------      
        if (yrecldp->nssopt == 0)
        {
          // No scheme
          zqe = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsice[jk-1][jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
          zqe = fmax(0.0, zqe);
        } else {
          if (yrecldp->nssopt == 1)
          {
            // Tompkins 
            zqe = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsice[jk-1][jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
            zqe = fmax(0.0, zqe);
          } else {
            if (yrecldp->nssopt == 2)
            {
              // Lohmann and Karcher
              zqe = zqx[5-1][jk-1][jl-1];
            } else {
              if (yrecldp->nssopt == 3)
              {
                // Gierens
                zqe = zqx[5-1][jk-1][jl-1] + zli[jk-1][jl-1];
              }

            }

          }

        }

        if (yrecldp->nssopt == 0 || ztp1[jk-1][jl-1] >= rtt)
        {
          // No ice supersaturation allowed
          zfac = 1.0;
        } else {
          // Ice supersaturation
          zfac = zfokoop[jl-1];
        }

        if (zqe >= zqsice[jk-1][jl-1]*zfac*zrhc && zqe < zqsice[jk-1][jl-1]*zfac)
        {
          // note: not **2 on 1-a term if ZQE is used. 
          // Added correction term ZFAC to numerator 15/03/2010
          zacond = -((1.0 - za[jk-1][jl-1])*zfac)*zdqs[jl-1]*1.0/fmax(2.0*(zfac*zqsice[jk-1][jl-1] - zqe), zepsec);
          zacond = fmin(zacond, 1.0 - za[jk-1][jl-1]);            // Linear term:
          // Added correction term ZFAC 15/03/2010
          zlcond2[jl-1] = -(zfac*zdqs[jl-1])*0.5*zacond;            // new limiter formulation
          zzdl = (2.0*(zfac*zqsice[jk-1][jl-1] - zqe))*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);            // Added correction term ZFAC 15/03/2010
          if (zdqs[jl-1]*zfac < -zzdl)
          {
            // ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
            zlcondlim = ((za[jk-1][jl-1] - 1.0)*zfac)*zdqs[jl-1] - zfac*zqsice[jk-1][jl-1] + zqx[5-1][jk-1][jl-1];
            zlcond2[jl-1] = fmin(zlcond2[jl-1], zlcondlim);
          }

          zlcond2[jl-1] = fmax(zlcond2[jl-1], 0.0);
          if (1.0 - za[jk-1][jl-1] < zepsec || zlcond2[jl-1] < yrecldp->rlmin)
          {
            zlcond2[jl-1] = 0.0;
            zacond = 0.0;
          }

          if (zlcond2[jl-1] == 0.0)
          {
            zacond = 0.0;
          }

          // Large-scale generation is LINEAR in A and LINEAR in L
          zsolac[jl-1] = zsolac[jl-1] + zacond;            // Store cloud fraction diagnostic if required
          //------------------------------------------------------------------------
          // All increase goes into liquid unless so cold cloud homogeneously freezes
          // Include new liquid formation in first guess value, otherwise liquid 
          // remains at cold temperatures until next timestep.
          //------------------------------------------------------------------------
          if (ztp1[jk-1][jl-1] > yrecldp->rthomo)
          {
            zsolqa[5-1][1-1][jl-1] = zsolqa[5-1][1-1][jl-1] + zlcond2[jl-1];
            zsolqa[1-1][5-1][jl-1] = zsolqa[1-1][5-1][jl-1] - zlcond2[jl-1];
            zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] + zlcond2[jl-1];              // Store cloud liquid diagnostic if required
          } else {
            zsolqa[5-1][2-1][jl-1] = zsolqa[5-1][2-1][jl-1] + zlcond2[jl-1];
            zsolqa[2-1][5-1][jl-1] = zsolqa[2-1][5-1][jl-1] - zlcond2[jl-1];
            zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + zlcond2[jl-1];              // Store cloud ice diagnostic if required
          }

        }

      }

    }

    //----------------------------------------------------------------------
    // 3.7 Growth of ice by vapour deposition 
    //----------------------------------------------------------------------
    // Following Rotstayn et al. 2001:
    // does not use the ice nuclei number from cloudaer.F90
    // but rather a simple Meyers et al. 1992 form based on the 
    // supersaturation and assuming clouds are saturated with 
    // respect to liquid water (well mixed), (or Koop adjustment)
    // Growth considered as sink of liquid water if present so 
    // Bergeron-Findeisen adjustment in autoconversion term no longer needed
    //----------------------------------------------------------------------
    //--------------------------------------------------------
    //-
    //- Ice deposition following Rotstayn et al. (2001)
    //-  (monodisperse ice particle size distribution)
    //-
    //--------------------------------------------------------
    if (idepice == 1)
    {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	//--------------------------------------------------------------
	// Calculate distance from cloud top 
	// defined by cloudy layer below a layer with cloud frac <0.01
	// ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
	//--------------------------------------------------------------
	if (za[jk-1][jl-1] >= yrecldp->rcldtopcf && za[jk-1-1][jl-1] < yrecldp->rcldtopcf)
	{
	  zcldtopdist[jl-1] = 0.0;
	} else {
	  zcldtopdist[jl-1] = zcldtopdist[jl-1] + zdp[jl-1]/((zrho[jl-1]*rg));
	}

	//--------------------------------------------------------------
	// only treat depositional growth if liquid present. due to fact 
	// that can not model ice growth from vapour without additional 
	// in-cloud water vapour variable
	//--------------------------------------------------------------
	if (zqxfg[1-1][jl-1] > yrecldp->rlmin && ztp1[jk-1][jl-1] < rtt)
	{
	  zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies)))*rv)/rd;
	  zvpliq = zvpice*zfokoop[jl-1];
	  zicenuclei[jl-1] = 1000.0*exp((12.96*(zvpliq - zvpice))/zvpliq - 0.639);              //------------------------------------------------
	  //   2.4e-2 is conductivity of air
	  //   8.8 = 700**1/3 = density of ice to the third
	  //------------------------------------------------
	  zadd = (rlstt*(rlstt/((rv*ztp1[jk-1][jl-1])) - 1.0))/((0.024*ztp1[jk-1][jl-1]));
	  zbdd = ((rv*ztp1[jk-1][jl-1])*pap[jk-1][jl-1])/((2.21*zvpice));
	  zcvds = ((7.8*pow(zicenuclei[jl-1]/zrho[jl-1], 0.666))*(zvpliq - zvpice))/(((8.87*(zadd + zbdd))*zvpice));              //-----------------------------------------------------
	  // RICEINIT=1.E-12_JPRB is initial mass of ice particle
	  //-----------------------------------------------------
	  zice0 = fmax(zicecld[jl-1], (zicenuclei[jl-1]*yrecldp->riceinit)/zrho[jl-1]);              //------------------
	  // new value of ice:
	  //------------------
	  zinew = pow((0.666*zcvds)*ptsphy + pow(zice0, 0.666), 1.5);              //---------------------------
	  // grid-mean deposition rate:
	  //--------------------------- 
	  zdepos = fmax(za[jk-1][jl-1]*(zinew - zice0), 0.0);              //--------------------------------------------------------------------
	  // Limit deposition to liquid water amount
	  // If liquid is all frozen, ice would use up reservoir of water 
	  // vapour in excess of ice saturation mixing ratio - However this 
	  // can not be represented without a in-cloud humidity variable. Using 
	  // the grid-mean humidity would imply a large artificial horizontal 
	  // flux from the clear sky to the cloudy area. We thus rely on the 
	  // supersaturation check to clean up any remaining supersaturation
	  //--------------------------------------------------------------------
	  zdepos = fmin(zdepos, zqxfg[1-1][jl-1]);              //--------------------------------------------------------------------
	  // At top of cloud, reduce deposition rate near cloud top to account for
	  // small scale turbulent processes, limited ice nucleation and ice fallout 
	  //--------------------------------------------------------------------
	  //      ZDEPOS = ZDEPOS*MIN(RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH,1.0_JPRB)
	  // Change to include dependence on ice nuclei concentration
	  // to increase deposition rate with decreasing temperatures 
	  zinfactor = fmin(zicenuclei[jl-1]/15000.0, 1.0);
	  zdepos = zdepos*fmin(zinfactor + (1.0 - zinfactor)*(yrecldp->rdepliqrefrate + zcldtopdist[jl-1]/yrecldp->rdepliqrefdepth), 1.0);              //--------------
	  // add to matrix 
	  //--------------
	  zsolqa[1-1][2-1][jl-1] = zsolqa[1-1][2-1][jl-1] + zdepos;
	  zsolqa[2-1][1-1][jl-1] = zsolqa[2-1][1-1][jl-1] - zdepos;
	  zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + zdepos;
	  zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] - zdepos;              // Store cloud budget diagnostics if required
	}

      }

      //--------------------------------------------------------
      //-
      //- Ice deposition assuming ice PSD
      //-
      //--------------------------------------------------------
    } else {
      if (idepice == 2)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  //--------------------------------------------------------------
	  // Calculate distance from cloud top 
	  // defined by cloudy layer below a layer with cloud frac <0.01
	  // ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
	  //--------------------------------------------------------------
	  if (za[jk-1][jl-1] >= yrecldp->rcldtopcf && za[jk-1-1][jl-1] < yrecldp->rcldtopcf)
	  {
	    zcldtopdist[jl-1] = 0.0;
	  } else {
	    zcldtopdist[jl-1] = zcldtopdist[jl-1] + zdp[jl-1]/((zrho[jl-1]*rg));
	  }

	  //--------------------------------------------------------------
	  // only treat depositional growth if liquid present. due to fact 
	  // that can not model ice growth from vapour without additional 
	  // in-cloud water vapour variable
	  //--------------------------------------------------------------
	  if (zqxfg[1-1][jl-1] > yrecldp->rlmin && ztp1[jk-1][jl-1] < rtt)
	  {
	    zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies)))*rv)/rd;
	    zvpliq = zvpice*zfokoop[jl-1];
	    zicenuclei[jl-1] = 1000.0*exp((12.96*(zvpliq - zvpice))/zvpliq - 0.639);                //-----------------------------------------------------
	    // RICEINIT=1.E-12_JPRB is initial mass of ice particle
	    //-----------------------------------------------------
	    zice0 = fmax(zicecld[jl-1], (zicenuclei[jl-1]*yrecldp->riceinit)/zrho[jl-1]);                // Particle size distribution
	    ztcg = 1.0;
	    zfacx1i = 1.0;
	    zaplusb = yrecldp->rcl_apb1*zvpice - yrecldp->rcl_apb2*zvpice*ztp1[jk-1][jl-1] + (pap[jk-1][jl-1]*yrecldp->rcl_apb3)*pow(ztp1[jk-1][jl-1], 3.0);
	    zcorrfac = sqrt(1.0/zrho[jl-1]);
	    zcorrfac2 = pow(ztp1[jk-1][jl-1]/273.0, 1.5)*(393.0/(ztp1[jk-1][jl-1] + 120.0));
	    zpr02 = ((zrho[jl-1]*zice0)*yrecldp->rcl_const1i)/((ztcg*zfacx1i));
	    zterm1 = (((((((zvpliq - zvpice)*pow(ztp1[jk-1][jl-1], 2.0))*zvpice)*zcorrfac2)*ztcg)*yrecldp->rcl_const2i)*zfacx1i)/(((zrho[jl-1]*zaplusb)*zvpice));
	    zterm2 = (0.65*yrecldp->rcl_const6i)*pow(zpr02, yrecldp->rcl_const4i) + (((yrecldp->rcl_const3i*sqrt(zcorrfac))*sqrt(zrho[jl-1]))*pow(zpr02, yrecldp->rcl_const5i))/sqrt(zcorrfac2);
	    zdepos = fmax(((za[jk-1][jl-1]*zterm1)*zterm2)*ptsphy, 0.0);                //--------------------------------------------------------------------
	    // Limit deposition to liquid water amount
	    // If liquid is all frozen, ice would use up reservoir of water 
	    // vapour in excess of ice saturation mixing ratio - However this 
	    // can not be represented without a in-cloud humidity variable. Using 
	    // the grid-mean humidity would imply a large artificial horizontal 
	    // flux from the clear sky to the cloudy area. We thus rely on the 
	    // supersaturation check to clean up any remaining supersaturation
	    //--------------------------------------------------------------------
	    zdepos = fmin(zdepos, zqxfg[1-1][jl-1]);                //--------------------------------------------------------------------
	    // At top of cloud, reduce deposition rate near cloud top to account for
	    // small scale turbulent processes, limited ice nucleation and ice fallout 
	    //--------------------------------------------------------------------
	    // Change to include dependence on ice nuclei concentration
	    // to increase deposition rate with decreasing temperatures 
	    zinfactor = fmin(zicenuclei[jl-1]/15000.0, 1.0);
	    zdepos = zdepos*fmin(zinfactor + (1.0 - zinfactor)*(yrecldp->rdepliqrefrate + zcldtopdist[jl-1]/yrecldp->rdepliqrefdepth), 1.0);                //--------------
	    // add to matrix 
	    //--------------
	    zsolqa[1-1][2-1][jl-1] = zsolqa[1-1][2-1][jl-1] + zdepos;
	    zsolqa[2-1][1-1][jl-1] = zsolqa[2-1][1-1][jl-1] - zdepos;
	    zqxfg[2-1][jl-1] = zqxfg[2-1][jl-1] + zdepos;
	    zqxfg[1-1][jl-1] = zqxfg[1-1][jl-1] - zdepos;                // Store cloud budget diagnostics if required
	  }

	}

      }

    }

    //######################################################################
    //              4  *** PRECIPITATION PROCESSES ***
    //######################################################################
    //----------------------------------
    // revise in-cloud consensate amount
    //----------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      ztmpa = 1.0*1.0/fmax(za[jk-1][jl-1], zepsec);
      zliqcld[jl-1] = zqxfg[1-1][jl-1]*ztmpa;
      zicecld[jl-1] = zqxfg[2-1][jl-1]*ztmpa;
      zlicld[jl-1] = zliqcld[jl-1] + zicecld[jl-1];
    }

    //----------------------------------------------------------------------
    // 4.2 SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
    //     now that rain, snow, graupel species are prognostic
    //     the precipitation flux can be defined directly level by level
    //     There is no vertical memory required from the flux variable
    //----------------------------------------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      if (llfall[jm-1] || jm == 2)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  //------------------------
	  // source from layer above 
	  //------------------------
	  if (jk > yrecldp->ncldtop)
	  {
	    zfallsrce[jm-1][jl-1] = zpfplsx[jm-1][jk-1][jl-1]*zdtgdp[jl-1];
	    zsolqa[jm-1][jm-1][jl-1] = zsolqa[jm-1][jm-1][jl-1] + zfallsrce[jm-1][jl-1];
	    zqxfg[jm-1][jl-1] = zqxfg[jm-1][jl-1] + zfallsrce[jm-1][jl-1];                // use first guess precip----------V
	    zqpretot[jl-1] = zqpretot[jl-1] + zqxfg[jm-1][jl-1];
	  }

	  //-------------------------------------------------
	  // sink to next layer, constant fall speed
	  //-------------------------------------------------
	  // if aerosol effect then override 
	  //  note that for T>233K this is the same as above.
	  if (yrecldp->laericesed && jm == 2)
	  {
	    zre_ice = pre_ice[jk-1][jl-1];                // The exponent value is from 
	    // Morrison et al. JAS 2005 Appendix
	    zvqx[2-1] = 0.002*pow(zre_ice, 1.0);
	  }

	  zfall = zvqx[jm-1]*zrho[jl-1];              //-------------------------------------------------
	  // modified by Heymsfield and Iaquinta JAS 2000
	  //-------------------------------------------------
	  zfallsink[jm-1][jl-1] = zdtgdp[jl-1]*zfall;              // Cloud budget diagnostic stored at end as implicit
	}

      }

    }

    //---------------------------------------------------------------
    // Precip cover overlap using MAX-RAN Overlap
    // Since precipitation is now prognostic we must 
    //   1) apply an arbitrary minimum coverage (0.3) if precip>0
    //   2) abandon the 2-flux clr/cld treatment
    //   3) Thus, since we have no memory of the clear sky precip
    //      fraction, we mimic the previous method by reducing 
    //      ZCOVPTOT(JL), which has the memory, proportionally with 
    //      the precip evaporation rate, taking cloud fraction 
    //      into account
    //   #3 above leads to much smoother vertical profiles of 
    //   precipitation fraction than the Klein-Jakob scheme which 
    //   monotonically increases precip fraction and then resets 
    //   it to zero in a step function once clear-sky precip reaches
    //   zero.
    //---------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zqpretot[jl-1] > zepsec)
      {
        zcovptot[jl-1] = 1.0 - (1.0 - zcovptot[jl-1])*(1.0 - fmax(za[jk-1][jl-1], za[jk-1-1][jl-1]))*1.0/(1.0 - fmin(za[jk-1-1][jl-1], 1.0 - 1.0e-6));
        zcovptot[jl-1] = fmax(zcovptot[jl-1], yrecldp->rcovpmin);
        zcovpclr[jl-1] = fmax(0.0, zcovptot[jl-1] - za[jk-1][jl-1]);
        zraincld[jl-1] = zqxfg[3-1][jl-1]/zcovptot[jl-1];
        zsnowcld[jl-1] = zqxfg[4-1][jl-1]/zcovptot[jl-1];
        zcovpmax[jl-1] = fmax(zcovptot[jl-1], zcovpmax[jl-1]);
      } else {
        zraincld[jl-1] = 0.0;
        zsnowcld[jl-1] = 0.0;
        zcovptot[jl-1] = 0.0;
        zcovpclr[jl-1] = 0.0;
        zcovpmax[jl-1] = 0.0;
      }

    }

    //----------------------------------------------------------------------
    // 4.3a AUTOCONVERSION TO SNOW
    //----------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (ztp1[jk-1][jl-1] <= rtt)
      {
	//-----------------------------------------------------
	//     Snow Autoconversion rate follow Lin et al. 1983
	//-----------------------------------------------------
	if (zicecld[jl-1] > zepsec)
	{
	  zzco = (ptsphy*yrecldp->rsnowlin1)*exp(yrecldp->rsnowlin2*(ztp1[jk-1][jl-1] - rtt));
	  if (yrecldp->laericeauto)
	  {
	    zlcrit = picrit_aer[jk-1][jl-1];                // 0.3 = N**0.333 with N=0.027 
	    zzco = zzco*pow(yrecldp->rnice/pnice[jk-1][jl-1], 0.333);
	  } else {
	    zlcrit = yrecldp->rlcritsnow;
	  }

	  zsnowaut[jl-1] = zzco*(1.0 - exp(-pow(zicecld[jl-1]/zlcrit, 2)));
	  zsolqb[2-1][4-1][jl-1] = zsolqb[2-1][4-1][jl-1] + zsnowaut[jl-1];
	}

      }

      //----------------------------------------------------------------------
      // 4.3b AUTOCONVERSION WARM CLOUDS
      //   Collection and accretion will require separate treatment
      //   but for now we keep this simple treatment
      //----------------------------------------------------------------------
      if (zliqcld[jl-1] > zepsec)
      {
	//--------------------------------------------------------
	//-
	//- Warm-rain process follow Sundqvist (1989)
	//-
	//--------------------------------------------------------
	if (iwarmrain == 1)
	{
	  zzco = yrecldp->rkconv*ptsphy;
	  if (yrecldp->laerliqautolsp)
	  {
	    zlcrit = plcrit_aer[jk-1][jl-1];                // 0.3 = N**0.333 with N=125 cm-3 
	    zzco = zzco*pow(yrecldp->rccn/pccn[jk-1][jl-1], 0.333);
	  } else {
	    // Modify autoconversion threshold dependent on: 
	    //  land (polluted, high CCN, smaller droplets, higher threshold)
	    //  sea  (clean, low CCN, larger droplets, lower threshold)
	    if (plsm[jl-1] > 0.5)
	    {
	      zlcrit = yrecldp->rclcrit_land;
	    } else {
	      zlcrit = yrecldp->rclcrit_sea;
	    }

	  }

	  //------------------------------------------------------------------
	  // Parameters for cloud collection by rain and snow.
	  // Note that with new prognostic variable it is now possible 
	  // to REPLACE this with an explicit collection parametrization
	  //------------------------------------------------------------------   
	  zprecip = (zpfplsx[4-1][jk-1][jl-1] + zpfplsx[3-1][jk-1][jl-1])*1.0/fmax(zepsec, zcovptot[jl-1]);
	  zcfpr = 1.0 + yrecldp->rprc1*sqrt(fmax(zprecip, 0.0));              //      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
	  //       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))
	  if (yrecldp->laerliqcoll)
	  {
	    // 5.0 = N**0.333 with N=125 cm-3 
	    zcfpr = zcfpr*pow(yrecldp->rccn/pccn[jk-1][jl-1], 0.333);
	  }

	  zzco = zzco*zcfpr;
	  zlcrit = zlcrit*1.0/fmax(zcfpr, zepsec);
	  if (zliqcld[jl-1]/zlcrit < 20.0)
	  {
	    zrainaut[jl-1] = zzco*(1.0 - exp(-pow(zliqcld[jl-1]/zlcrit, 2)));
	  } else {
	    zrainaut[jl-1] = zzco;
	  }

	  // rain freezes instantly
	  if (ztp1[jk-1][jl-1] <= rtt)
	  {
	    zsolqb[1-1][4-1][jl-1] = zsolqb[1-1][4-1][jl-1] + zrainaut[jl-1];
	  } else {
	    zsolqb[1-1][3-1][jl-1] = zsolqb[1-1][3-1][jl-1] + zrainaut[jl-1];
	  }

	  //--------------------------------------------------------
	  //-
	  //- Warm-rain process follow Khairoutdinov and Kogan (2000)
	  //-
	  //--------------------------------------------------------
	} else {
	  if (iwarmrain == 2)
	  {
	    if (plsm[jl-1] > 0.5)
	    {
	      zconst = yrecldp->rcl_kk_cloud_num_land;
	      zlcrit = yrecldp->rclcrit_land;
	    } else {
	      zconst = yrecldp->rcl_kk_cloud_num_sea;
	      zlcrit = yrecldp->rclcrit_sea;
	    }

	    if (zliqcld[jl-1] > zlcrit)
	    {
	      zrainaut[jl-1] = ((((1.5*za[jk-1][jl-1])*ptsphy)*yrecldp->rcl_kkaau)*pow(zliqcld[jl-1], yrecldp->rcl_kkbauq))*pow(zconst, yrecldp->rcl_kkbaun);
	      zrainaut[jl-1] = fmin(zrainaut[jl-1], zqxfg[1-1][jl-1]);
	      if (zrainaut[jl-1] < zepsec)
	      {
		zrainaut[jl-1] = 0.0;
	      }

	      zrainacc[jl-1] = (((2.0*za[jk-1][jl-1])*ptsphy)*yrecldp->rcl_kkaac)*pow(zliqcld[jl-1]*zraincld[jl-1], yrecldp->rcl_kkbac);
	      zrainacc[jl-1] = fmin(zrainacc[jl-1], zqxfg[1-1][jl-1]);
	      if (zrainacc[jl-1] < zepsec)
	      {
		zrainacc[jl-1] = 0.0;
	      }

	    } else {
	      zrainaut[jl-1] = 0.0;
	      zrainacc[jl-1] = 0.0;
	    }

	    // If temperature < 0, then autoconversion produces snow rather than rain
	    // Explicit
	    if (ztp1[jk-1][jl-1] <= rtt)
	    {
	      zsolqa[1-1][4-1][jl-1] = zsolqa[1-1][4-1][jl-1] + zrainaut[jl-1];
	      zsolqa[1-1][4-1][jl-1] = zsolqa[1-1][4-1][jl-1] + zrainacc[jl-1];
	      zsolqa[4-1][1-1][jl-1] = zsolqa[4-1][1-1][jl-1] - zrainaut[jl-1];
	      zsolqa[4-1][1-1][jl-1] = zsolqa[4-1][1-1][jl-1] - zrainacc[jl-1];
	    } else {
	      zsolqa[1-1][3-1][jl-1] = zsolqa[1-1][3-1][jl-1] + zrainaut[jl-1];
	      zsolqa[1-1][3-1][jl-1] = zsolqa[1-1][3-1][jl-1] + zrainacc[jl-1];
	      zsolqa[3-1][1-1][jl-1] = zsolqa[3-1][1-1][jl-1] - zrainaut[jl-1];
	      zsolqa[3-1][1-1][jl-1] = zsolqa[3-1][1-1][jl-1] - zrainacc[jl-1];
	    }

	  }

	}

      }

    }

    //----------------------------------------------------------------------
    // RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
    //      only active if T<0degC and supercooled liquid water is present
    //      AND if not Sundquist autoconversion (as this includes riming)
    //----------------------------------------------------------------------
    if (iwarmrain > 1)
    {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	if (ztp1[jk-1][jl-1] <= rtt && zliqcld[jl-1] > zepsec)
	{
	  // Fallspeed air density correction 
	  // TODO: THIS IS A BUG! Due to a missing ``_JPRB`` in the original,
	  // we need to cast the exponent down to single precision to re-create.
	  zfallcorr = pow(yrecldp->rdensref/zrho[jl-1], (float)0.4);
	  //------------------------------------------------------------------
	  // Riming of snow by cloud water - implicit in lwc
	  //------------------------------------------------------------------
	  if (zcovptot[jl-1] > 0.01 && zsnowcld[jl-1] > zepsec)
	  {
	    // Calculate riming term
	    // Factor of liq water taken out because implicit
	    zsnowrime[jl-1] = ((((0.3*zcovptot[jl-1])*ptsphy)*yrecldp->rcl_const7s)*zfallcorr)*pow((zrho[jl-1]*zsnowcld[jl-1])*yrecldp->rcl_const1s, yrecldp->rcl_const8s);                // Limit snow riming term
	    zsnowrime[jl-1] = fmin(zsnowrime[jl-1], 1.0);
	    zsolqb[1-1][4-1][jl-1] = zsolqb[1-1][4-1][jl-1] + zsnowrime[jl-1];
	  }

	  //------------------------------------------------------------------
	  // Riming of ice by cloud water - implicit in lwc
	  // NOT YET ACTIVE
	  //------------------------------------------------------------------
	  //      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
	  //
	  //        // Calculate riming term
	  //        // Factor of liq water taken out because implicit
	  //        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
	  //     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
	  //
	  //        // Limit ice riming term
	  //        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
	  //
	  //        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
	  //
	  //      ENDIF
	}

      }

    }

    //----------------------------------------------------------------------
    // 4.4a  MELTING OF SNOW and ICE
    //       with new implicit solver this also has to treat snow or ice
    //       precipitating from the level above... i.e. local ice AND flux.
    //       in situ ice and snow: could arise from LS advection or warming
    //       falling ice and snow: arrives by precipitation process
    //----------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zicetot[jl-1] = zqxfg[2-1][jl-1] + zqxfg[4-1][jl-1];
      zmeltmax[jl-1] = 0.0;          // If there are frozen hydrometeors present and dry-bulb temperature > 0degC
      if (zicetot[jl-1] > zepsec && ztp1[jk-1][jl-1] > rtt)
      {
	// Calculate subsaturation
	zsubsat = fmax(zqsice[jk-1][jl-1] - zqx[5-1][jk-1][jl-1], 0.0);            // Calculate difference between dry-bulb (ZTP1) and the temperature 
	// at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
	// Melting only occurs if the wet-bulb temperature >0
	// i.e. warming of ice particle due to melting > cooling 
	// due to evaporation.
	ztdmtw0 = ztp1[jk-1][jl-1] - rtt - zsubsat*(ztw1 + ztw2*(pap[jk-1][jl-1] - ztw3) - ztw4*(ztp1[jk-1][jl-1] - ztw5));            // Not implicit yet... 
	// Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
	zcons1 = fabs((ptsphy*(1.0 + 0.5*ztdmtw0))/yrecldp->rtaumel);
	zmeltmax[jl-1] = fmax((ztdmtw0*zcons1)*zrldcp, 0.0);
      }

    }

    // Loop over frozen hydrometeors (ice, snow)
    for (jm=1; jm<=5; jm+=1) {
      if (iphase[jm-1] == 2)
      {
	jn = imelt[jm-1];
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  if (zicetot[jl-1] > zepsec && zmeltmax[jl-1] > zepsec)
	  {
	    // Apply melting in same proportion as frozen hydrometeor fractions 
	    zalfa = zqxfg[jm-1][jl-1]/zicetot[jl-1];
	    zmelt = fmin(zqxfg[jm-1][jl-1], zalfa*zmeltmax[jl-1]);                // needed in first guess
	    // This implies that zqpretot has to be recalculated below
	    // since is not conserved here if ice falls and liquid doesn't
	    zqxfg[jm-1][jl-1] = zqxfg[jm-1][jl-1] - zmelt;
	    zqxfg[jn-1][jl-1] = zqxfg[jn-1][jl-1] + zmelt;
	    zsolqa[jm-1][jn-1][jl-1] = zsolqa[jm-1][jn-1][jl-1] + zmelt;
	    zsolqa[jn-1][jm-1][jl-1] = zsolqa[jn-1][jm-1][jl-1] - zmelt;
	  }

	}

      }

    }

    //----------------------------------------------------------------------
    // 4.4b  FREEZING of RAIN
    //----------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      // If rain present
      if (zqx[3-1][jk-1][jl-1] > zepsec)
      {
	if (ztp1[jk-1][jl-1] <= rtt && ztp1[jk-1-1][jl-1] > rtt)
	{
	  // Base of melting layer/top of refreezing layer so
	  // store rain/snow fraction for precip type diagnosis
	  // If mostly rain, then supercooled rain slow to freeze
	  // otherwise faster to freeze (snow or ice pellets)
	  zqpretot[jl-1] = fmax(zqx[4-1][jk-1][jl-1] + zqx[3-1][jk-1][jl-1], zepsec);
	  prainfrac_toprfz[jl-1] = zqx[3-1][jk-1][jl-1]/zqpretot[jl-1];
	}

	// If temperature less than zero
	if (ztp1[jk-1][jl-1] < rtt)
	{
	  if (prainfrac_toprfz[jl-1] > 0.8)
	  {
	    // Majority of raindrops completely melted
	    // Refreezing is by slow heterogeneous freezing
	    // Slope of rain particle size distribution
	    zlambda = pow(yrecldp->rcl_fac1/((zrho[jl-1]*zqx[3-1][jk-1][jl-1])), yrecldp->rcl_fac2);                // Calculate freezing rate based on Bigg(1953) and Wisner(1972)
	    ztemp = yrecldp->rcl_fzrab*(ztp1[jk-1][jl-1] - rtt);
	    zfrz = ((ptsphy*(yrecldp->rcl_const5r/zrho[jl-1]))*(exp(ztemp) - 1.0))*pow(zlambda, yrecldp->rcl_const6r);
	    zfrzmax[jl-1] = fmax(zfrz, 0.0);
	  } else {
	    // Majority of raindrops only partially melted 
	    // Refreeze with a shorter timescale (reverse of melting...for now)
	    zcons1 = fabs((ptsphy*(1.0 + 0.5*(rtt - ztp1[jk-1][jl-1])))/yrecldp->rtaumel);
	    zfrzmax[jl-1] = fmax(((rtt - ztp1[jk-1][jl-1])*zcons1)*zrldcp, 0.0);
	  }

	  if (zfrzmax[jl-1] > zepsec)
	  {
	    zfrz = fmin(zqx[3-1][jk-1][jl-1], zfrzmax[jl-1]);
	    zsolqa[3-1][4-1][jl-1] = zsolqa[3-1][4-1][jl-1] + zfrz;
	    zsolqa[4-1][3-1][jl-1] = zsolqa[4-1][3-1][jl-1] - zfrz;
	  }

	}

      }

    }

    //----------------------------------------------------------------------
    // 4.4c  FREEZING of LIQUID 
    //----------------------------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      // not implicit yet... 
      zfrzmax[jl-1] = fmax((yrecldp->rthomo - ztp1[jk-1][jl-1])*zrldcp, 0.0);
    }

    jm = 1;
    jn = imelt[jm-1];
    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zfrzmax[jl-1] > zepsec && zqxfg[jm-1][jl-1] > zepsec)
      {
	zfrz = fmin(zqxfg[jm-1][jl-1], zfrzmax[jl-1]);
	zsolqa[jm-1][jn-1][jl-1] = zsolqa[jm-1][jn-1][jl-1] + zfrz;
	zsolqa[jn-1][jm-1][jl-1] = zsolqa[jn-1][jm-1][jl-1] - zfrz;
      }

    }

    //----------------------------------------------------------------------
    // 4.5   EVAPORATION OF RAIN/SNOW
    //----------------------------------------------------------------------
    //----------------------------------------
    // Rain evaporation scheme from Sundquist
    //----------------------------------------
    if (ievaprain == 1)
    {
      // Rain
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zzrh = yrecldp->rprecrhmax + ((1.0 - yrecldp->rprecrhmax)*zcovpmax[jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
	zzrh = fmin(fmax(zzrh, yrecldp->rprecrhmax), 1.0);
	zqe = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsliq[jk-1][jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);            //---------------------------------------------
	// humidity in moistest ZCOVPCLR part of domain
	//---------------------------------------------
	zqe = fmax(0.0, fmin(zqe, zqsliq[jk-1][jl-1]));
	llo1 = zcovpclr[jl-1] > zepsec && zqxfg[3-1][jl-1] > zepsec && zqe < zzrh*zqsliq[jk-1][jl-1];
	if (llo1)
	{
	  // note: zpreclr is a rain flux
	  zpreclr = (zqxfg[3-1][jl-1]*zcovpclr[jl-1])*1.0/copysign(fmax(fabs(zcovptot[jl-1]*zdtgdp[jl-1]), zepsilon), zcovptot[jl-1]*zdtgdp[jl-1]);              //--------------------------------------
	  // actual microphysics formula in zbeta
	  //--------------------------------------
	  zbeta1 = ((sqrt(pap[jk-1][jl-1]/paph[klev+1-1][jl-1])/yrecldp->rvrfactor)*zpreclr)*1.0/fmax(zcovpclr[jl-1], zepsec);
	  zbeta = ((rg*yrecldp->rpecons)*0.5)*pow(zbeta1, 0.5777);
	  zdenom = 1.0 + (zbeta*ptsphy)*zcorqsliq[jl-1];
	  zdpr = ((((zcovpclr[jl-1]*zbeta)*(zqsliq[jk-1][jl-1] - zqe))/zdenom)*zdp[jl-1])*zrg_r;
	  zdpevap = zdpr*zdtgdp[jl-1];              //---------------------------------------------------------
	  // add evaporation term to explicit sink.
	  // this has to be explicit since if treated in the implicit
	  // term evaporation can not reduce rain to zero and model
	  // produces small amounts of rainfall everywhere. 
	  //---------------------------------------------------------
	  // Evaporate rain
	  zevap = fmin(zdpevap, zqxfg[3-1][jl-1]);
	  zsolqa[3-1][5-1][jl-1] = zsolqa[3-1][5-1][jl-1] + zevap;
	  zsolqa[5-1][3-1][jl-1] = zsolqa[5-1][3-1][jl-1] - zevap;              //IF (LLCLDBUDL) ZBUDL(JL,19)=-ZEVAP*ZQTMST
	  //-------------------------------------------------------------
	  // Reduce the total precip coverage proportional to evaporation
	  // to mimic the previous scheme which had a diagnostic
	  // 2-flux treatment, abandoned due to the new prognostic precip
	  //-------------------------------------------------------------
	  zcovptot[jl-1] = fmax(yrecldp->rcovpmin, zcovptot[jl-1] - fmax(0.0, ((zcovptot[jl-1] - za[jk-1][jl-1])*zevap)/zqxfg[3-1][jl-1]));              // Update fg field
	  zqxfg[3-1][jl-1] = zqxfg[3-1][jl-1] - zevap;
	}

      }

      //---------------------------------------------------------
      // Rain evaporation scheme based on Abel and Boutle (2013)
      //---------------------------------------------------------
    } else {
      if (ievaprain == 2)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  //-----------------------------------------------------------------------
	  // Calculate relative humidity limit for rain evaporation 
	  // to avoid cloud formation and saturation of the grid box
	  //-----------------------------------------------------------------------
	  // Limit RH for rain evaporation dependent on precipitation fraction 
	  zzrh = yrecldp->rprecrhmax + ((1.0 - yrecldp->rprecrhmax)*zcovpmax[jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
	  zzrh = fmin(fmax(zzrh, yrecldp->rprecrhmax), 1.0);              // Critical relative humidity
	  //ZRHC=RAMID
	  //ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
	  // Increase RHcrit to 1.0 towards the surface (eta>0.8)
	  //IF(ZSIGK > 0.8_JPRB) THEN
	  //  ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
	  //ENDIF
	  //ZZRH = MIN(ZRHC,ZZRH)
	  // Further limit RH for rain evaporation to 80% (RHcrit in free troposphere)
	  zzrh = fmin(0.8, zzrh);
	  zqe = fmax(0.0, fmin(zqx[5-1][jk-1][jl-1], zqsliq[jk-1][jl-1]));
	  llo1 = zcovpclr[jl-1] > zepsec && zqxfg[3-1][jl-1] > zepsec && zqe < zzrh*zqsliq[jk-1][jl-1];
	  if (llo1)
	  {
	    //-------------------------------------------
	    // Abel and Boutle (2012) evaporation
	    //-------------------------------------------
	    // Calculate local precipitation (kg/kg)
	    zpreclr = zqxfg[3-1][jl-1]/zcovptot[jl-1];                // Fallspeed air density correction 
	    zfallcorr = pow(yrecldp->rdensref/zrho[jl-1], 0.4);                // Saturation vapour pressure with respect to liquid phase
	    zesatliq = (rv/rd)*(double)(r2es*exp((r3les*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4les)));                // Slope of particle size distribution
	    zlambda = pow(yrecldp->rcl_fac1/((zrho[jl-1]*zpreclr)), yrecldp->rcl_fac2);
	    zevap_denom = yrecldp->rcl_cdenom1*zesatliq - yrecldp->rcl_cdenom2*ztp1[jk-1][jl-1]*zesatliq + (yrecldp->rcl_cdenom3*pow(ztp1[jk-1][jl-1], 3.0))*pap[jk-1][jl-1];                // Temperature dependent conductivity
	    zcorr2 = (pow(ztp1[jk-1][jl-1]/273.0, 1.5)*393.0)/(ztp1[jk-1][jl-1] + 120.0);
	    zka = yrecldp->rcl_ka273*zcorr2;
	    zsubsat = fmax(zzrh*zqsliq[jk-1][jl-1] - zqe, 0.0);
	    zbeta = (((((0.5/zqsliq[jk-1][jl-1])*pow(ztp1[jk-1][jl-1], 2.0))*zesatliq)*yrecldp->rcl_const1r)*(zcorr2/zevap_denom))*(0.78/pow(zlambda, yrecldp->rcl_const4r) + (yrecldp->rcl_const2r*sqrt(zrho[jl-1]*zfallcorr))/((sqrt(zcorr2)*pow(zlambda, yrecldp->rcl_const3r))));
	    zdenom = 1.0 + zbeta*ptsphy;
	    zdpevap = (((zcovpclr[jl-1]*zbeta)*ptsphy)*zsubsat)/zdenom;                //---------------------------------------------------------
	    // Add evaporation term to explicit sink.
	    // this has to be explicit since if treated in the implicit
	    // term evaporation can not reduce rain to zero and model
	    // produces small amounts of rainfall everywhere. 
	    //---------------------------------------------------------
	    // Limit rain evaporation
	    zevap = fmin(zdpevap, zqxfg[3-1][jl-1]);
	    zsolqa[3-1][5-1][jl-1] = zsolqa[3-1][5-1][jl-1] + zevap;
	    zsolqa[5-1][3-1][jl-1] = zsolqa[5-1][3-1][jl-1] - zevap;                //IF (LLCLDBUDL) ZBUDL(JL,19)=-ZEVAP*ZQTMST
	    //-------------------------------------------------------------
	    // Reduce the total precip coverage proportional to evaporation
	    // to mimic the previous scheme which had a diagnostic
	    // 2-flux treatment, abandoned due to the new prognostic precip
	    //-------------------------------------------------------------
	    zcovptot[jl-1] = fmax(yrecldp->rcovpmin, zcovptot[jl-1] - fmax(0.0, ((zcovptot[jl-1] - za[jk-1][jl-1])*zevap)/zqxfg[3-1][jl-1]));                // Update fg field 
	    zqxfg[3-1][jl-1] = zqxfg[3-1][jl-1] - zevap;
	  }

	}

      }

    }

    //----------------------------------------------------------------------
    // 4.5   EVAPORATION OF SNOW
    //----------------------------------------------------------------------
    // Snow
    if (ievapsnow == 1)
    {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zzrh = yrecldp->rprecrhmax + ((1.0 - yrecldp->rprecrhmax)*zcovpmax[jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
	zzrh = fmin(fmax(zzrh, yrecldp->rprecrhmax), 1.0);
	zqe = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsice[jk-1][jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);            //---------------------------------------------
	// humidity in moistest ZCOVPCLR part of domain
	//---------------------------------------------
	zqe = fmax(0.0, fmin(zqe, zqsice[jk-1][jl-1]));
	llo1 = zcovpclr[jl-1] > zepsec && zqxfg[4-1][jl-1] > zepsec && zqe < zzrh*zqsice[jk-1][jl-1];
	if (llo1)
	{
	  // note: zpreclr is a rain flux a
	  zpreclr = (zqxfg[4-1][jl-1]*zcovpclr[jl-1])*1.0/copysign(fmax(fabs(zcovptot[jl-1]*zdtgdp[jl-1]), zepsilon), zcovptot[jl-1]*zdtgdp[jl-1]);              //--------------------------------------
	  // actual microphysics formula in zbeta
	  //--------------------------------------
	  zbeta1 = ((sqrt(pap[jk-1][jl-1]/paph[klev+1-1][jl-1])/yrecldp->rvrfactor)*zpreclr)*1.0/fmax(zcovpclr[jl-1], zepsec);
	  zbeta = (rg*yrecldp->rpecons)*pow(zbeta1, 0.5777);
	  zdenom = 1.0 + (zbeta*ptsphy)*zcorqsice[jl-1];
	  zdpr = ((((zcovpclr[jl-1]*zbeta)*(zqsice[jk-1][jl-1] - zqe))/zdenom)*zdp[jl-1])*zrg_r;
	  zdpevap = zdpr*zdtgdp[jl-1];              //---------------------------------------------------------
	  // add evaporation term to explicit sink.
	  // this has to be explicit since if treated in the implicit
	  // term evaporation can not reduce snow to zero and model
	  // produces small amounts of snowfall everywhere. 
	  //---------------------------------------------------------
	  // Evaporate snow
	  zevap = fmin(zdpevap, zqxfg[4-1][jl-1]);
	  zsolqa[4-1][5-1][jl-1] = zsolqa[4-1][5-1][jl-1] + zevap;
	  zsolqa[5-1][4-1][jl-1] = zsolqa[5-1][4-1][jl-1] - zevap;              //IF (LLCLDBUDL) ZBUDi(JL,17)=-ZEVAP*ZQTMST
	  //-------------------------------------------------------------
	  // Reduce the total precip coverage proportional to evaporation
	  // to mimic the previous scheme which had a diagnostic
	  // 2-flux treatment, abandoned due to the new prognostic precip
	  //-------------------------------------------------------------
	  zcovptot[jl-1] = fmax(yrecldp->rcovpmin, zcovptot[jl-1] - fmax(0.0, ((zcovptot[jl-1] - za[jk-1][jl-1])*zevap)/zqxfg[4-1][jl-1]));              //Update first guess field
	  zqxfg[4-1][jl-1] = zqxfg[4-1][jl-1] - zevap;
	}

      }

      //---------------------------------------------------------
    } else {
      if (ievapsnow == 2)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  //-----------------------------------------------------------------------
	  // Calculate relative humidity limit for snow evaporation 
	  //-----------------------------------------------------------------------
	  zzrh = yrecldp->rprecrhmax + ((1.0 - yrecldp->rprecrhmax)*zcovpmax[jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);
	  zzrh = fmin(fmax(zzrh, yrecldp->rprecrhmax), 1.0);
	  zqe = (zqx[5-1][jk-1][jl-1] - za[jk-1][jl-1]*zqsice[jk-1][jl-1])*1.0/fmax(zepsec, 1.0 - za[jk-1][jl-1]);              //---------------------------------------------
	  // humidity in moistest ZCOVPCLR part of domain
	  //---------------------------------------------
	  zqe = fmax(0.0, fmin(zqe, zqsice[jk-1][jl-1]));
	  llo1 = zcovpclr[jl-1] > zepsec && zqx[4-1][jk-1][jl-1] > zepsec && zqe < zzrh*zqsice[jk-1][jl-1];
	  if (llo1)
	  {
	    // Calculate local precipitation (kg/kg)
	    zpreclr = zqx[4-1][jk-1][jl-1]/zcovptot[jl-1];
	    zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1][jl-1] - rtt))/(ztp1[jk-1][jl-1] - r4ies)))*rv)/rd;                // Particle size distribution
	    // ZTCG increases Ni with colder temperatures - essentially a 
	    // Fletcher or Meyers scheme? 
	    ztcg = 1.0;                // ZFACX1I modification is based on Andrew Barrett's results
	    zfacx1s = 1.0;
	    zaplusb = yrecldp->rcl_apb1*zvpice - yrecldp->rcl_apb2*zvpice*ztp1[jk-1][jl-1] + (pap[jk-1][jl-1]*yrecldp->rcl_apb3)*pow(ztp1[jk-1][jl-1], 3);
	    zcorrfac = sqrt(1.0/zrho[jl-1]);
	    zcorrfac2 = pow(ztp1[jk-1][jl-1]/273.0, 1.5)*(393.0/(ztp1[jk-1][jl-1] + 120.0));
	    zpr02 = ((zrho[jl-1]*zpreclr)*yrecldp->rcl_const1s)/((ztcg*zfacx1s));
	    zterm1 = (((((((zqsice[jk-1][jl-1] - zqe)*pow(ztp1[jk-1][jl-1], 2))*zvpice)*zcorrfac2)*ztcg)*yrecldp->rcl_const2s)*zfacx1s)/(((zrho[jl-1]*zaplusb)*zqsice[jk-1][jl-1]));
	    zterm2 = (0.65*yrecldp->rcl_const6s)*pow(zpr02, yrecldp->rcl_const4s) + (((yrecldp->rcl_const3s*sqrt(zcorrfac))*sqrt(zrho[jl-1]))*pow(zpr02, yrecldp->rcl_const5s))/sqrt(zcorrfac2);
	    zdpevap = fmax(((zcovpclr[jl-1]*zterm1)*zterm2)*ptsphy, 0.0);                //--------------------------------------------------------------------
	    // Limit evaporation to snow amount
	    //--------------------------------------------------------------------
	    zevap = fmin(zdpevap, zevaplimice[jl-1]);
	    zevap = fmin(zevap, zqx[4-1][jk-1][jl-1]);
	    zsolqa[4-1][5-1][jl-1] = zsolqa[4-1][5-1][jl-1] + zevap;
	    zsolqa[5-1][4-1][jl-1] = zsolqa[5-1][4-1][jl-1] - zevap;                //IF (LLCLDBUDL) ZBUDi(JL,17)=-ZEVAP*ZQTMST
	    //-------------------------------------------------------------
	    // Reduce the total precip coverage proportional to evaporation
	    // to mimic the previous scheme which had a diagnostic
	    // 2-flux treatment, abandoned due to the new prognostic precip
	    //-------------------------------------------------------------
	    zcovptot[jl-1] = fmax(yrecldp->rcovpmin, zcovptot[jl-1] - fmax(0.0, ((zcovptot[jl-1] - za[jk-1][jl-1])*zevap)/zqx[4-1][jk-1][jl-1]));                //Update first guess field
	    zqxfg[4-1][jl-1] = zqxfg[4-1][jl-1] - zevap;
	  }

	}

      }

    }

    //--------------------------------------
    // Evaporate small precipitation amounts
    //--------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      if (llfall[jm-1])
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  if (zqxfg[jm-1][jl-1] < yrecldp->rlmin)
	  {
	    zsolqa[jm-1][5-1][jl-1] = zsolqa[jm-1][5-1][jl-1] + zqxfg[jm-1][jl-1];
	    zsolqa[5-1][jm-1][jl-1] = zsolqa[5-1][jm-1][jl-1] - zqxfg[jm-1][jl-1];
	  }

	}

      }

    }

    //######################################################################
    //            5.0  *** SOLVERS FOR A AND L ***
    // now use an implicit solution rather than exact solution
    // solver is forward in time, upstream difference for advection
    //######################################################################
    //---------------------------
    // 5.1 solver for cloud cover
    //---------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zanew = (za[jk-1][jl-1] + zsolac[jl-1])/(1.0 + zsolab[jl-1]);
      zanew = fmin(zanew, 1.0);
      if (zanew < yrecldp->ramin)
      {
	zanew = 0.0;
      }

      zda[jl-1] = zanew - zaorig[jk-1][jl-1];          //---------------------------------
      // variables needed for next level
      //---------------------------------
      zanewm1[jl-1] = zanew;
    }

    //--------------------------------
    // 5.2 solver for the microphysics
    //--------------------------------
    //--------------------------------------------------------------
    // Truncate explicit sinks to avoid negatives 
    // Note: Species are treated in the order in which they run out
    // since the clipping will alter the balance for the other vars
    //--------------------------------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jn=1; jn<=5; jn+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  llindex3[jm-1][jn-1][jl-1] = false;
	}

      }

      for (jl=kidia; jl<=kfdia; jl+=1) {
	zsinksum[jm-1][jl-1] = 0.0;
      }

    }

    //----------------------------
    // collect sink terms and mark
    //----------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jn=1; jn<=5; jn+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  zsinksum[jm-1][jl-1] = zsinksum[jm-1][jl-1] - zsolqa[jn-1][jm-1][jl-1];
	}

      }

    }

    //---------------------------------------
    // calculate overshoot and scaling factor
    //---------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zmax = fmax(zqx[jm-1][jk-1][jl-1], zepsec);
	zrat = fmax(zsinksum[jm-1][jl-1], zmax);
	zratio[jm-1][jl-1] = zmax/zrat;
      }

    }

    //--------------------------------------------
    // scale the sink terms, in the correct order, 
    // recalculating the scale factor each time
    //--------------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zsinksum[jm-1][jl-1] = 0.0;
      }

    }

    //----------------
    // recalculate sum
    //----------------
    for (jm=1; jm<=5; jm+=1) {
      for (i_klon=1; i_klon<=klon; i_klon++) {
	psum_solqa[i_klon-1] = 0.0;
      }

      for (jn=1; jn<=5; jn+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  psum_solqa[jl-1] = psum_solqa[jl-1] + zsolqa[jn-1][jm-1][jl-1];
	}

      }

      for (jl=kidia; jl<=kfdia; jl+=1) {
	// ZSINKSUM(JL,JM)=ZSINKSUM(JL,JM)-SUM(ZSOLQA(JL,JM,1:NCLV))
	zsinksum[jm-1][jl-1] = zsinksum[jm-1][jl-1] - psum_solqa[jl-1];
      }

      //---------------------------
      // recalculate scaling factor
      //---------------------------
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zmm = fmax(zqx[jm-1][jk-1][jl-1], zepsec);
	zrr = fmax(zsinksum[jm-1][jl-1], zmm);
	zratio[jm-1][jl-1] = zmm/zrr;
      }

      //------
      // scale
      //------
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zzratio = zratio[jm-1][jl-1];            //DIR$ IVDEP
	//DIR$ PREFERVECTOR
	for (jn=1; jn<=5; jn+=1) {
	  if (zsolqa[jn-1][jm-1][jl-1] < 0.0)
	  {
	    zsolqa[jn-1][jm-1][jl-1] = zsolqa[jn-1][jm-1][jl-1]*zzratio;
	    zsolqa[jm-1][jn-1][jl-1] = zsolqa[jm-1][jn-1][jl-1]*zzratio;
	  }

	}

      }

    }

    //--------------------------------------------------------------
    // 5.2.2 Solver
    //------------------------
    //------------------------
    // set the LHS of equation  
    //------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jn=1; jn<=5; jn+=1) {
	//----------------------------------------------
	// diagonals: microphysical sink terms+transport
	//----------------------------------------------
	if (jn == jm)
	{
	  for (jl=kidia; jl<=kfdia; jl+=1) {
	    zqlhs[jm-1][jn-1][jl-1] = 1.0 + zfallsink[jm-1][jl-1];
	    for (jo=1; jo<=5; jo+=1) {
	      zqlhs[jm-1][jn-1][jl-1] = zqlhs[jm-1][jn-1][jl-1] + zsolqb[jn-1][jo-1][jl-1];
	    }

	  }

	  //------------------------------------------
	  // non-diagonals: microphysical source terms
	  //------------------------------------------
	} else {
	  for (jl=kidia; jl<=kfdia; jl+=1) {
	    zqlhs[jm-1][jn-1][jl-1] = -zsolqb[jm-1][jn-1][jl-1];
	  }

	}

      }

    }

    //------------------------
    // set the RHS of equation  
    //------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	//---------------------------------
	// sum the explicit source and sink
	//---------------------------------
	zexplicit = 0.0;
	for (jn=1; jn<=5; jn+=1) {
	  zexplicit = zexplicit + zsolqa[jn-1][jm-1][jl-1];
	}

	zqxn[jm-1][jl-1] = zqx[jm-1][jk-1][jl-1] + zexplicit;
      }

    }

    //-----------------------------------
    // *** solve by LU decomposition: ***
    //-----------------------------------
    // Note: This fast way of solving NCLVxNCLV system
    //       assumes a good behaviour (i.e. non-zero diagonal
    //       terms with comparable orders) of the matrix stored
    //       in ZQLHS. For the moment this is the case but
    //       be aware to preserve it when doing eventual 
    //       modifications.
    // Non pivoting recursive factorization 
    for (jn=1; jn<=4; jn+=1) {
      for (jm=jn + 1; jm<=5; jm+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  zqlhs[jn-1][jm-1][jl-1] = zqlhs[jn-1][jm-1][jl-1]/zqlhs[jn-1][jn-1][jl-1];
	}

	for (ik=jn + 1; ik<=5; ik+=1) {
	  for (jl=kidia; jl<=kfdia; jl+=1) {
	    zqlhs[ik-1][jm-1][jl-1] = zqlhs[ik-1][jm-1][jl-1] - zqlhs[jn-1][jm-1][jl-1]*zqlhs[ik-1][jn-1][jl-1];
	  }

	}

      }

    }

    // Backsubstitution 
    //  step 1 
    for (jn=2; jn<=5; jn+=1) {
      for (jm=1; jm<=jn - 1; jm+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  zqxn[jn-1][jl-1] = zqxn[jn-1][jl-1] - zqlhs[jm-1][jn-1][jl-1]*zqxn[jm-1][jl-1];
	}

      }

    }

    //  step 2
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zqxn[5-1][jl-1] = zqxn[5-1][jl-1]/zqlhs[5-1][5-1][jl-1];
    }

    for (jn=4; jn>=1; jn+=-1) {
      for (jm=jn + 1; jm<=5; jm+=1) {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  zqxn[jn-1][jl-1] = zqxn[jn-1][jl-1] - zqlhs[jm-1][jn-1][jl-1]*zqxn[jm-1][jl-1];
	}

      }

      for (jl=kidia; jl<=kfdia; jl+=1) {
	zqxn[jn-1][jl-1] = zqxn[jn-1][jl-1]/zqlhs[jn-1][jn-1][jl-1];
      }

    }

    // Ensure no small values (including negatives) remain in cloud variables nor
    // precipitation rates.
    // Evaporate l,i,r,s to water vapour. Latent heating taken into account below
    for (jn=1; jn<=4; jn+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	if (zqxn[jn-1][jl-1] < zepsec)
	{
	  zqxn[5-1][jl-1] = zqxn[5-1][jl-1] + zqxn[jn-1][jl-1];
	  zqxn[jn-1][jl-1] = 0.0;
	}

      }

    }

    //--------------------------------
    // variables needed for next level
    //--------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zqxnm1[jm-1][jl-1] = zqxn[jm-1][jl-1];
	zqxn2d[jm-1][jk-1][jl-1] = zqxn[jm-1][jl-1];
      }

    }

    //------------------------------------------------------------------------
    // 5.3 Precipitation/sedimentation fluxes to next level
    //     diagnostic precipitation fluxes
    //     It is this scaled flux that must be used for source to next layer
    //------------------------------------------------------------------------
    for (jm=1; jm<=5; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	zpfplsx[jm-1][jk+1-1][jl-1] = (zfallsink[jm-1][jl-1]*zqxn[jm-1][jl-1])*zrdtgdp[jl-1];
      }

    }

    // Ensure precipitation fraction is zero if no precipitation
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zqpretot[jl-1] = zpfplsx[4-1][jk+1-1][jl-1] + zpfplsx[3-1][jk+1-1][jl-1];
    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      if (zqpretot[jl-1] < zepsec)
      {
	zcovptot[jl-1] = 0.0;
      }

    }

    //######################################################################
    //              6  *** UPDATE TENDANCIES ***
    //######################################################################
    //--------------------------------
    // 6.1 Temperature and CLV budgets 
    //--------------------------------
    for (jm=1; jm<=4; jm+=1) {
      for (jl=kidia; jl<=kfdia; jl+=1) {
	// calculate fluxes in and out of box for conservation of TL
	zfluxq[jm-1][jl-1] = zpsupsatsrce[jm-1][jl-1] + zconvsrce[jm-1][jl-1] + zfallsrce[jm-1][jl-1] - (zfallsink[jm-1][jl-1] + zconvsink[jm-1][jl-1])*zqxn[jm-1][jl-1];
      }

      if (iphase[jm-1] == 1)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] + (ralvdcp*(zqxn[jm-1][jl-1] - zqx[jm-1][jk-1][jl-1] - zfluxq[jm-1][jl-1]))*zqtmst;
	}

      }

      if (iphase[jm-1] == 2)
      {
	for (jl=kidia; jl<=kfdia; jl+=1) {
	  tendency_loc_t[jk-1][jl-1] = tendency_loc_t[jk-1][jl-1] + (ralsdcp*(zqxn[jm-1][jl-1] - zqx[jm-1][jk-1][jl-1] - zfluxq[jm-1][jl-1]))*zqtmst;
	}

      }

      //----------------------------------------------------------------------
      // New prognostic tendencies - ice,liquid rain,snow 
      // Note: CLV arrays use PCLV in calculation of tendency while humidity
      //       uses ZQX. This is due to clipping at start of cloudsc which
      //       include the tendency already in tendency_loc%T and tendency_loc%q. ZQX was reset
      //----------------------------------------------------------------------
      for (jl=kidia; jl<=kfdia; jl+=1) {
	tendency_loc_cld[jm-1][jk-1][jl-1] = tendency_loc_cld[jm-1][jk-1][jl-1] + (zqxn[jm-1][jl-1] - zqx0[jm-1][jk-1][jl-1])*zqtmst;
      }

    }

    for (jl=kidia; jl<=kfdia; jl+=1) {
      //----------------------
      // 6.2 Humidity budget
      //----------------------
      tendency_loc_q[jk-1][jl-1] = tendency_loc_q[jk-1][jl-1] + (zqxn[5-1][jl-1] - zqx[5-1][jk-1][jl-1])*zqtmst;
      //-------------------
      // 6.3 cloud cover 
      //-----------------------
      tendency_loc_a[jk-1][jl-1] = tendency_loc_a[jk-1][jl-1] + zda[jl-1]*zqtmst;
    }

    //--------------------------------------------------
    // Copy precipitation fraction into output variable
    //-------------------------------------------------
    for (jl=kidia; jl<=kfdia; jl+=1) {
      pcovptot[jk-1][jl-1] = zcovptot[jl-1];
    }

  }

  //----------------------------------------------------------------------
  //                       END OF VERTICAL LOOP
  //----------------------------------------------------------------------
  //######################################################################
  //              8  *** FLUX/DIAGNOSTICS COMPUTATIONS ***
  //######################################################################

  //--------------------------------------------------------------------
  // Copy general precip arrays back into PFP arrays for GRIB archiving
  // Add rain and liquid fluxes, ice and snow fluxes
  //--------------------------------------------------------------------
  for (jk=1; jk<=klev + 1; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      pfplsl[jk-1][jl-1] = zpfplsx[3-1][jk-1][jl-1] + zpfplsx[1-1][jk-1][jl-1];
      pfplsn[jk-1][jl-1] = zpfplsx[4-1][jk-1][jl-1] + zpfplsx[2-1][jk-1][jl-1];
    }

  }

  //--------
  // Fluxes:
  //--------
  for (jl=kidia; jl<=kfdia; jl+=1) {
    pfsqlf[1-1][jl-1] = 0.0;
    pfsqif[1-1][jl-1] = 0.0;
    pfsqrf[1-1][jl-1] = 0.0;
    pfsqsf[1-1][jl-1] = 0.0;
    pfcqlng[1-1][jl-1] = 0.0;
    pfcqnng[1-1][jl-1] = 0.0;
    pfcqrng[1-1][jl-1] = 0.0;
    pfcqsng[1-1][jl-1] = 0.0;        // fluxes due to turbulence
    pfsqltur[1-1][jl-1] = 0.0;
    pfsqitur[1-1][jl-1] = 0.0;
  }

  for (jk=1; jk<=klev; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      zgdph_r = -zrg_r*(paph[jk+1-1][jl-1] - paph[jk-1][jl-1])*zqtmst;
      pfsqlf[jk+1-1][jl-1] = pfsqlf[jk-1][jl-1];
      pfsqif[jk+1-1][jl-1] = pfsqif[jk-1][jl-1];
      pfsqrf[jk+1-1][jl-1] = pfsqlf[jk-1][jl-1];
      pfsqsf[jk+1-1][jl-1] = pfsqif[jk-1][jl-1];
      pfcqlng[jk+1-1][jl-1] = pfcqlng[jk-1][jl-1];
      pfcqnng[jk+1-1][jl-1] = pfcqnng[jk-1][jl-1];
      pfcqrng[jk+1-1][jl-1] = pfcqlng[jk-1][jl-1];
      pfcqsng[jk+1-1][jl-1] = pfcqnng[jk-1][jl-1];
      pfsqltur[jk+1-1][jl-1] = pfsqltur[jk-1][jl-1];
      pfsqitur[jk+1-1][jl-1] = pfsqitur[jk-1][jl-1];
      zalfaw = zfoealfa[jk-1][jl-1];          // Liquid , LS scheme minus detrainment
      pfsqlf[jk+1-1][jl-1] = pfsqlf[jk+1-1][jl-1] + (zqxn2d[1-1][jk-1][jl-1] - zqx0[1-1][jk-1][jl-1] + pvfl[jk-1][jl-1]*ptsphy - zalfaw*plude[jk-1][jl-1])*zgdph_r;          // liquid, negative numbers
      pfcqlng[jk+1-1][jl-1] = pfcqlng[jk+1-1][jl-1] + zlneg[1-1][jk-1][jl-1]*zgdph_r;          // liquid, vertical diffusion
      pfsqltur[jk+1-1][jl-1] = pfsqltur[jk+1-1][jl-1] + (pvfl[jk-1][jl-1]*ptsphy)*zgdph_r;          // Rain, LS scheme 
      pfsqrf[jk+1-1][jl-1] = pfsqrf[jk+1-1][jl-1] + (zqxn2d[3-1][jk-1][jl-1] - zqx0[3-1][jk-1][jl-1])*zgdph_r;          // rain, negative numbers
      pfcqrng[jk+1-1][jl-1] = pfcqrng[jk+1-1][jl-1] + zlneg[3-1][jk-1][jl-1]*zgdph_r;          // Ice , LS scheme minus detrainment
      pfsqif[jk+1-1][jl-1] = pfsqif[jk+1-1][jl-1] + (zqxn2d[2-1][jk-1][jl-1] - zqx0[2-1][jk-1][jl-1] + pvfi[jk-1][jl-1]*ptsphy - (1.0 - zalfaw)*plude[jk-1][jl-1])*zgdph_r;          // ice, negative numbers
      pfcqnng[jk+1-1][jl-1] = pfcqnng[jk+1-1][jl-1] + zlneg[2-1][jk-1][jl-1]*zgdph_r;          // ice, vertical diffusion
      pfsqitur[jk+1-1][jl-1] = pfsqitur[jk+1-1][jl-1] + (pvfi[jk-1][jl-1]*ptsphy)*zgdph_r;          // snow, LS scheme
      pfsqsf[jk+1-1][jl-1] = pfsqsf[jk+1-1][jl-1] + (zqxn2d[4-1][jk-1][jl-1] - zqx0[4-1][jk-1][jl-1])*zgdph_r;          // snow, negative numbers
      pfcqsng[jk+1-1][jl-1] = pfcqsng[jk+1-1][jl-1] + zlneg[4-1][jk-1][jl-1]*zgdph_r;
    }

  }

  //-----------------------------------
  // enthalpy flux due to precipitation
  //-----------------------------------
  for (jk=1; jk<=klev + 1; jk+=1) {
    for (jl=kidia; jl<=kfdia; jl+=1) {
      pfhpsl[jk-1][jl-1] = -rlvtt*pfplsl[jk-1][jl-1];
      pfhpsn[jk-1][jl-1] = -rlstt*pfplsn[jk-1][jl-1];
    }

  }

  //===============================================================================
  //IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)
  return 0;
}
