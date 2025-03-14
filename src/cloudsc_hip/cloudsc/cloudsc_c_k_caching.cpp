#include "hip/hip_runtime.h"
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
#include <float.h>

__global__ void cloudsc_c(int kidia, int kfdia, int klon, dtype ptsphy,
  const dtype * __restrict__  pt,
  const dtype * __restrict__  pq, const dtype * __restrict__  tendency_tmp_t,
  const dtype * __restrict__  tendency_tmp_q, const dtype * __restrict__  tendency_tmp_a,
  const dtype * __restrict__  tendency_tmp_cld, dtype * __restrict__  tendency_loc_t,
  dtype * __restrict__  tendency_loc_q, dtype * __restrict__  tendency_loc_a,
  dtype * __restrict__  tendency_loc_cld, const dtype * __restrict__  pvfa,
  const dtype * __restrict__  pvfl, const dtype * __restrict__  pvfi, const dtype * __restrict__  pdyna,
  const dtype * __restrict__  pdynl, const dtype * __restrict__  pdyni, const dtype * __restrict__  phrsw,
  dtype * __restrict__  phrlw, const dtype * __restrict__  pvervel, const dtype * __restrict__  pap,
  const dtype * __restrict__  paph, const dtype * __restrict__  plsm,
  const int *  ktype, const dtype * __restrict__  plu, dtype * __restrict__  plude,
  const dtype * __restrict__  psnde, const dtype * __restrict__  pmfu, const dtype * __restrict__  pmfd,
  const dtype * __restrict__  pa, const dtype * __restrict__  pclv, const dtype * __restrict__  psupsat,
  const dtype * __restrict__  plcrit_aer, const dtype * __restrict__  picrit_aer,
  const dtype * __restrict__  pre_ice, const dtype * __restrict__  pccn, const dtype * __restrict__  pnice,
  dtype * __restrict__  pcovptot, dtype * __restrict__  prainfrac_toprfz,
  dtype * __restrict__  pfsqlf, dtype * __restrict__  pfsqif, dtype * __restrict__  pfcqnng,
  dtype * __restrict__  pfcqlng, dtype * __restrict__  pfsqrf, dtype * __restrict__  pfsqsf,
  dtype * __restrict__  pfcqrng, dtype * __restrict__  pfcqsng,
  dtype * __restrict__  pfsqltur, dtype * __restrict__  pfsqitur,
  dtype * __restrict__  pfplsl, dtype * __restrict__  pfplsn, dtype * __restrict__  pfhpsl,
  dtype * __restrict__  pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  dtype rg, dtype rd, dtype rcpd, dtype retv, dtype rlvtt, dtype rlstt, dtype rlmlt, dtype rtt,
  dtype rv, dtype r2es, dtype r3les, dtype r3ies, dtype r4les, dtype r4ies, dtype r5les,
  dtype r5ies, dtype r5alvcp, dtype r5alscp, dtype ralvdcp, dtype ralsdcp, dtype ralfdcp,
  dtype rtwat, dtype rtice, dtype rticecu, dtype rtwat_rtice_r, dtype rtwat_rticecu_r,
  dtype rkoop1, dtype rkoop2) {

  //-------------------------------------------------------------------------------
  //                 Declare input/output arguments
  //-------------------------------------------------------------------------------

  // PLCRIT_AER : critical liquid mmr for rain autoconversion process
  // PICRIT_AER : critical liquid mmr for snow autoconversion process
  // PRE_LIQ : liq Re
  // PRE_ICE : ice Re
  // PCCN    : liquid cloud condensation nuclei
  // PNICE   : ice number concentration (cf. CCN)

  const int klev = 137;  // Number of levels

  dtype zlcond1, zlcond2, zlevapl, zlevapi, zrainaut, zsnowaut, zliqcld, zicecld;
  dtype zlevap, zleros;
  //  condensation and evaporation terms
  // autoconversion terms
  dtype zfokoop;
  dtype zfoealfa;
  dtype zicenuclei;  // number concentration of ice nuclei

  dtype zlicld;
  dtype zacond;
  dtype zaeros;
  dtype zlfinalsum;
  dtype zdqs;
  dtype ztold;
  dtype zqold;
  dtype zdtgdp;
  dtype zrdtgdp;
  dtype ztrpaus;
  dtype zcovpclr;
  dtype zpreclr;
  dtype zcovptot;
  dtype zcovpmax;
  dtype zqpretot;
  dtype zdpevap;
  dtype zdtforc;
  dtype zdtdiab;
  dtype ztp1[2];
  dtype zldefr;
  dtype zldifdt;
  dtype zdtgdpf;
  dtype zlcust[5];
  dtype zacust;
  dtype zmf;

  dtype zrho;
  dtype ztmp1, ztmp2, ztmp3;
  dtype ztmp4, ztmp5, ztmp6, ztmp7;
  dtype zalfawm;

  // Accumulators of A,B,and C factors for cloud equations
  dtype zsolab;  // -ve implicit CC
  dtype zsolac;  // linear CC
  dtype zanew;
  dtype zanewm1;

  dtype zgdp;

  //---for flux calculation
  dtype zda;
  dtype zli;
  dtype za[2];
  dtype zaorig;

  int llflag;
  int llo1;

  int icall, ik, jk, jl, jm, jn, jo, jlen, is;

  dtype zdp, zpaphd;

  dtype zalfa;
  // & ZALFACU, ZALFALS
  dtype zalfaw;
  dtype zbeta, zbeta1;
  dtype zcfpr;
  dtype zcor;
  dtype zcdmax;
  dtype zmin;
  dtype zlcondlim;
  dtype zdenom;
  dtype zdpmxdt;
  dtype zdpr;
  dtype zdtdp;
  dtype ze;
  dtype zepsec;
  dtype zfac, zfaci, zfacw;
  dtype zgdcp;
  dtype zinew;
  dtype zlcrit;
  dtype zmfdn;
  dtype zprecip;
  dtype zqe;
  dtype zqsat, zqtmst, zrdcp;
  dtype zrhc, zsig, zsigk;
  dtype zwtot;
  dtype zzco, zzdl, zzrh, zzzdt, zqadj;
  dtype zqnew, ztnew;
  dtype zrg_r, zgdph_r, zcons1, zcond, zcons1a;
  dtype zlfinal;
  dtype zmelt;
  dtype zevap;
  dtype zfrz;
  dtype zvpliq, zvpice;
  dtype zadd, zbdd, zcvds, zice0, zdepos;
  dtype zsupsat;
  dtype zfall;
  dtype zre_ice;
  dtype zrldcp;
  dtype zqp1env;

  //----------------------------
  // Arrays for new microphysics
  //----------------------------
  int iphase[5];  // marker for water phase of each species
  // 0=vapour, 1=liquid, 2=ice

  int imelt[5];  // marks melting linkage for ice categories
  // ice->liquid, snow->rain

  int llfall[5];  // marks falling species
  // LLFALL=0, cloud cover must > 0 for zqx > 0
  // LLFALL=1, no cloud needed, zqx can evaporate

  int llindex1[5];  // index variable
  int llindex3[5*5];  // index variable
  dtype zmax;
  dtype zrat;
  int iorder[5];  // array for sorting explicit terms

  dtype zliqfrac;
  dtype zicefrac;
  dtype zqx[5];
  dtype zqx0[5];
  dtype zqxn[5];  // new values for zqx at time+1
  dtype zqxfg[5];  // first guess values including precip
  dtype zqxnm1[5];  // new values for zqx at time+1 at level above
  dtype zfluxq[5];  // fluxes convergence of species (needed?)
  // Keep the following for possible future total water variance scheme?
  //REAL(KIND=JPRB) :: ZTL(KLON,KLEV)       ! liquid water temperature
  //REAL(KIND=JPRB) :: ZABETA(KLON,KLEV)    ! cloud fraction
  //REAL(KIND=JPRB) :: ZVAR(KLON,KLEV)      ! temporary variance
  //REAL(KIND=JPRB) :: ZQTMIN(KLON,KLEV)
  //REAL(KIND=JPRB) :: ZQTMAX(KLON,KLEV)

  dtype zlneg[5];
  dtype zmeltmax;
  dtype zfrzmax;
  dtype zicetot;

  dtype zqxn2d[5];

  dtype zqsmix;
  //REAL(KIND=JPRB) :: ZQSBIN(KLON,KLEV) ! binary switched ice/liq saturation
  dtype zqsliq;
  dtype zqsice;

  //REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
  //REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
  //REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice

  dtype zfoeewmt;
  dtype zfoeew;
  dtype zfoeeliqt;

  //REAL(KIND=JPRB) :: ZFOEEICET(KLON,KLEV)

  dtype zdqsliqdt, zdqsicedt, zdqsmixdt;
  dtype zcorqsliq;
  dtype zcorqsice;
  //REAL(KIND=JPRB) :: ZCORQSBIN(KLON)
  dtype zcorqsmix;
  dtype zevaplimliq, zevaplimice, zevaplimmix;

  //-------------------------------------------------------
  // SOURCE/SINK array for implicit and explicit terms
  //-------------------------------------------------------
  // a POSITIVE value entered into the arrays is a...
  //            Source of this variable
  //            |
  //            |   Sink of this variable
  //            |   |
  //            V   V
  // ZSOLQA(JL,IQa,IQb)  = explicit terms
  // ZSOLQB(JL,IQa,IQb)  = implicit terms
  // Thus if ZSOLAB(JL,NCLDQL,IQV)=K where K>0 then this is
  // a source of NCLDQL and a sink of IQV
  // put 'magic' source terms such as PLUDE from
  // detrainment into explicit source/sink array diagnognal
  // ZSOLQA(NCLDQL,NCLDQL)= -PLUDE
  // i.e. A positive value is a sink!????? weird...
  //-------------------------------------------------------

  dtype zsolqa[5*5];  // explicit sources and sinks
  dtype zsolqb[5*5];  // implicit sources and sinks
  // e.g. microphysical pathways between ice variables.
  dtype zqlhs[5*5];  // n x n matrix storing the LHS of implicit solver
  dtype zvqx[5];  // fall speeds of three categories
  dtype zexplicit;
  dtype zratio[5], zsinksum[5];

  // for sedimentation source/sink terms
  dtype zfallsink[5];
  dtype zfallsrce[5];

  // for convection detrainment source and subsidence source/sink terms
  dtype zconvsrce[5];
  dtype zconvsink[5];

  // for supersaturation source term from previous timestep
  dtype zpsupsatsrce[5];

  // Numerical fit to wet bulb temperature
  dtype ztw1 = (dtype) 1329.31;
  dtype ztw2 = (dtype) 0.0074615;
  dtype ztw3 = (dtype) 0.85E5;
  dtype ztw4 = (dtype) 40.637;
  dtype ztw5 = (dtype) 275.0;

  dtype zsubsat;  // Subsaturation for snow melting term
  dtype ztdmtw0;  // Diff between dry-bulb temperature and
  // temperature when wet-bulb = 0degC

  // Variables for deposition term
  dtype ztcg;  // Temperature dependent function for ice PSD
  dtype zfacx1i, zfacx1s;  // PSD correction factor
  dtype zaplusb, zcorrfac, zcorrfac2, zpr02, zterm1, zterm2;  // for ice dep
  dtype zcldtopdist;  // Distance from cloud top
  dtype zinfactor;  // No. of ice nuclei factor for deposition

  // Autoconversion/accretion/riming/evaporation
  int iwarmrain;
  int ievaprain;
  int ievapsnow;
  int idepice;
  dtype zrainacc;
  dtype zraincld;
  dtype zsnowrime;
  dtype zsnowcld;
  dtype zesatliq;
  dtype zfallcorr;
  dtype zlambda;
  dtype zevap_denom;
  dtype zcorr2;
  dtype zka;
  dtype zconst;
  dtype ztemp;

  // Rain freezing
  int llrainliq;  // True if majority of raindrops are liquid (no ice core)

  //----------------------------
  // End: new microphysics
  //----------------------------

  //----------------------
  // SCM budget statistics
  //----------------------
  dtype zrain;

  dtype zhook_handle;
  dtype ztmpl, ztmpi, ztmpa;

  dtype zmm, zrr;
  dtype zrg;

  dtype zzsum, zzratio;
  dtype zepsilon;

  dtype zcond1, zqp;

  dtype psum_solqa;

  int ibl;
  int i_llfall_0;
  //dtype zqx[5 * klev];
  //dtype zqx0[5 * klev];
  dtype zpfplsx[5 * 2];
  //dtype zlneg[5 * klev];
  //dtype zqxn2d[5 * klev];

  jl = threadIdx.x; 
  ibl = blockIdx.x; 

  int jk_i;
  int jk_ip1;
  int jk_im1;

  zepsilon = (dtype) 100.*DBL_EPSILON;

  // ---------------------------------------------------------------------
  // Set version of warm-rain autoconversion/accretion
  // IWARMRAIN = 1 ! Sundquist
  // IWARMRAIN = 2 ! Khairoutdinov and Kogan (2000)
  // ---------------------------------------------------------------------
  iwarmrain = 2;
  // ---------------------------------------------------------------------
  // Set version of rain evaporation
  // IEVAPRAIN = 1 ! Sundquist
  // IEVAPRAIN = 2 ! Abel and Boutle (2013)
  // ---------------------------------------------------------------------
  ievaprain = 2;
  // ---------------------------------------------------------------------
  // Set version of snow evaporation
  // IEVAPSNOW = 1 ! Sundquist
  // IEVAPSNOW = 2 ! New
  // ---------------------------------------------------------------------
  ievapsnow = 1;
  // ---------------------------------------------------------------------
  // Set version of ice deposition
  // IDEPICE = 1 ! Rotstayn (2001)
  // IDEPICE = 2 ! New
  // ---------------------------------------------------------------------
  idepice = 1;

  // ---------------------
  // Some simple constants
  // ---------------------
  zqtmst = (dtype) 1.0 / ptsphy;
  zgdcp = rg / rcpd;
  zrdcp = rd / rcpd;
  zcons1a = rcpd / (rlmlt*rg*(*yrecldp).rtaumel);
  zepsec = (dtype) 1.E-14;
  zrg_r = (dtype) 1.0 / rg;
  zrldcp = (dtype) 1.0 / (ralsdcp - ralvdcp);

  // Note: Defined in module/yoecldp.F90
  // NCLDQL=1    ! liquid cloud water
  // NCLDQI=2    ! ice cloud water
  // NCLDQR=3    ! rain water
  // NCLDQS=4    ! snow
  // NCLDQV=5    ! vapour

  // -----------------------------------------------
  // Define species phase, 0=vapour, 1=liquid, 2=ice
  // -----------------------------------------------
  iphase[4] = 0;
  iphase[0] = 1;
  iphase[2] = 1;
  iphase[1] = 2;
  iphase[3] = 2;

  // ---------------------------------------------------
  // Set up melting/freezing index,
  // if an ice category melts/freezes, where does it go?
  // ---------------------------------------------------
  imelt[4] = -99;
  imelt[0] = 2;
  imelt[2] = 4;
  imelt[1] = 3;
  imelt[3] = 3;

  // -----------------------------------------------
  // INITIALIZATION OF OUTPUT TENDENCIES
  // -----------------------------------------------
  for (jk = 0; jk <= klev + -1; jk += 1) {
    tendency_loc_t[jl + klon*(jk + klev*ibl)] = (dtype) 0.0;
    tendency_loc_q[jl + klon*(jk + klev*ibl)] = (dtype) 0.0;
    tendency_loc_a[jl + klon*(jk + klev*ibl)] = (dtype) 0.0;
  }
  for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
    for (jk = 0; jk <= klev + -1; jk += 1) {
      tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*ibl))] = (dtype) 0.0;
    }
  }

  //-- These were uninitialized : meaningful only when we compare error differences
  for (jk = 0; jk <= klev + -1; jk += 1) {
    pcovptot[jl + klon*(jk + klev*ibl)] = (dtype) 0.0;
    tendency_loc_cld[jl + klon*(jk + klev*(5 - 1 + 5*(ibl)))] = (dtype) 0.0;
  }

  //--------
  // Fluxes:
  //--------
  pfsqlf[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfsqif[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfsqrf[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfsqsf[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfcqlng[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfcqnng[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfcqrng[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;      //rain
  pfcqsng[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;      //snow
  // fluxes due to turbulence
  pfsqltur[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;
  pfsqitur[jl + klon*(0 + (klev + 1)*ibl)] = (dtype) 0.0;

  // -------------------------
  // set up fall speeds in m/s
  // -------------------------
  zvqx[4] = (dtype) 0.0;
  zvqx[0] = (dtype) 0.0;
  zvqx[1] = (*yrecldp).rvice;
  zvqx[2] = (*yrecldp).rvrain;
  zvqx[3] = (*yrecldp).rvsnow;
  for (i_llfall_0 = 0; i_llfall_0 <= 5 + -1; i_llfall_0 += 1) {
    llfall[i_llfall_0] = false;
  }
  for (jm = 0; jm <= 5 + -1; jm += 1) {
    if (zvqx[jm] > (dtype) 0.0) {
      llfall[jm] = true;
    }
    // falling species
  }
  // Set LLFALL to false for ice (but ice still sediments!)
  // Need to rationalise this at some point
  llfall[1] = false;

  prainfrac_toprfz[jl + klon*ibl] = (dtype) 0.0;      // rain fraction at top of refreezing layer
  llrainliq = true;      // Assume all raindrops are liquid initially

  //######################################################################
  //             1.  *** INITIAL VALUES FOR VARIABLES ***
  //######################################################################

  //-----------------------------
  // Reset single level variables
  //-----------------------------

  zanewm1 = (dtype) 0.0;
  zda = (dtype) 0.0;
  zcovpclr = (dtype) 0.0;
  zcovpmax = (dtype) 0.0;
  zcovptot = (dtype) 0.0;
  zcldtopdist = (dtype) 0.0;
  //-------------
  // zero arrays
  //-------------
  for (jm = 0; jm <= 5 + -1; jm += 1) {
    zpfplsx[0 + 2*jm] = (dtype) 0.0;        // precip fluxes
    zpfplsx[1 + 2*jm] = (dtype) 0.0;
  }

  // ----------------------
  // non CLV initialization
  // ----------------------
  for (jk = 0; jk <= klev + 1 + -1; jk += 1) {

    // Fortran counting is beautiful!
    jk_i = (jk + 1) % 2;
    jk_ip1 = (jk + 2) % 2;
    jk_im1 = (jk) % 2; 

    if (1 <= jk + 1 && jk + 1 <= klev) {
      ztp1[jk_i] = pt[jl + klon*(jk + klev*ibl)] + ptsphy*
          tendency_tmp_t[jl + klon*(jk + klev*ibl)];
      zqx[4] = pq[jl + klon*(jk + klev*ibl)] + ptsphy*
          tendency_tmp_q[jl + klon*(jk + klev*ibl)];
      zqx0[4] = pq[jl + klon*(jk + klev*ibl)] + ptsphy*
          tendency_tmp_q[jl + klon*(jk + klev*ibl)];
      za[jk_i] = pa[jl + klon*(jk + klev*ibl)] + ptsphy*
          tendency_tmp_a[jl + klon*(jk + klev*ibl)];
      zaorig = pa[jl + klon*(jk + klev*ibl)] + ptsphy*
          tendency_tmp_a[jl + klon*(jk + klev*ibl)];

      // -------------------------------------
      // initialization for CLV family
      // -------------------------------------
      for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
        zqx[jm] = pclv[jl + klon*(jk + klev*(jm + 5*ibl))] + ptsphy*
            tendency_tmp_cld[jl + klon*(jk + klev*(jm + 5*ibl))];
        zqx0[jm] = pclv[jl + klon*(jk + klev*(jm + 5*ibl))] + ptsphy*
            tendency_tmp_cld[jl + klon*(jk + klev*(jm + 5*ibl))];
      }

      for (jm = 0; jm <= 5 + -1; jm += 1) {
        zqxn2d[jm] = (dtype) 0.0;            // end of timestep values in 2D
        zlneg[jm] = (dtype) 0.0;            // negative input check
      }


      // ----------------------------------------------------
      // Tidy up very small cloud cover or total cloud water
      // ----------------------------------------------------
      if (zqx[0] + zqx[1] < (*yrecldp).rlmin || za[jk_i] < (*yrecldp)
        .ramin) {

        // Evaporate small cloud liquid water amounts
        zlneg[0] = zlneg[0] + zqx[0];
        zqadj = zqx[0]*zqtmst;
        tendency_loc_q[jl + klon*(jk + klev*ibl)] =
          tendency_loc_q[jl + klon*(jk + klev*ibl)] + zqadj;
        tendency_loc_t[jl + klon*(jk + klev*ibl)] =
          tendency_loc_t[jl + klon*(jk + klev*ibl)] - ralvdcp*zqadj;
        zqx[4] = zqx[4] + zqx[0];
        zqx[0] = (dtype) 0.0;

        // Evaporate small cloud ice water amounts
        zlneg[1] = zlneg[1] + zqx[1];
        zqadj = zqx[1]*zqtmst;
        tendency_loc_q[jl + klon*(jk + klev*ibl)] =
          tendency_loc_q[jl + klon*(jk + klev*ibl)] + zqadj;
        tendency_loc_t[jl + klon*(jk + klev*ibl)] =
          tendency_loc_t[jl + klon*(jk + klev*ibl)] - ralsdcp*zqadj;
        zqx[4] = zqx[4] + zqx[1];
        zqx[1] = (dtype) 0.0;

        // Set cloud cover to zero
        za[jk_i] = (dtype) 0.0;

      }

      // ---------------------------------
      // Tidy up small CLV variables
      // ---------------------------------
      //DIR$ IVDEP
      for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
        if (zqx[jm] < (*yrecldp).rlmin) {
          zlneg[jm] = zlneg[jm] + zqx[jm];
          zqadj = zqx[jm]*zqtmst;
          tendency_loc_q[jl + klon*(jk + klev*ibl)] =
            tendency_loc_q[jl + klon*(jk + klev*ibl)] + zqadj;
          if (iphase[jm] == 1) {
            tendency_loc_t[jl + klon*(jk + klev*ibl)] =
              tendency_loc_t[jl + klon*(jk + klev*ibl)] - ralvdcp*zqadj;
          }
          if (iphase[jm] == 2) {
            tendency_loc_t[jl + klon*(jk + klev*ibl)] =
              tendency_loc_t[jl + klon*(jk + klev*ibl)] - ralsdcp*zqadj;
          }
          zqx[4] = zqx[4] + zqx[jm];
          zqx[jm] = (dtype) 0.0;
        }
      }

      // ------------------------------
      // Define saturation values
      // ------------------------------
      //----------------------------------------
      // old *diagnostic* mixed phase saturation
      //----------------------------------------
      zfoealfa = ((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))));
      zfoeewmt = MYMIN(((dtype)(r2es*((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*
        MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les)) + 
	(1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*
	MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))) / pap[jl + klon*(jk + klev*ibl)], (dtype) 0.5);
      zqsmix = zfoeewmt;
      zqsmix = zqsmix / ((dtype) 1.0 - retv*zqsmix);

      //---------------------------------------------
      // ice saturation T<273K
      // liquid water saturation for T>273K
      //---------------------------------------------
      zalfa = ((dtype)(MYMAX(0.0, copysign(1.0, ztp1[jk_i] - rtt))));
      zfoeew = MYMIN((zalfa*((dtype)(r2es*MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les)))) + 
        ((dtype) 1.0 - zalfa)*((dtype)(r2es*MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))) / 
        pap[jl + klon*(jk + klev*ibl)], (dtype) 0.5);
      zfoeew = MYMIN((dtype) 0.5, zfoeew);
      zqsice = zfoeew / ((dtype) 1.0 - retv*zfoeew);

      //----------------------------------
      // liquid water saturation
      //----------------------------------
      zfoeeliqt = MYMIN(((dtype)(r2es*MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les)))) / 
        pap[jl + klon*(jk + klev*ibl)], (dtype) 0.5);
      zqsliq = zfoeeliqt;
      zqsliq = zqsliq / ((dtype) 1.0 - retv*zqsliq);

      //   !----------------------------------
      //   ! ice water saturation
      //   !----------------------------------
      //   ZFOEEICET(JL,JK)=MIN(FOEEICE(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      //   ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
      //   ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-retv*ZQSICE(JL,JK))


      //------------------------------------------
      // Ensure cloud fraction is between 0 and 1
      //------------------------------------------
      za[jk_i] = MYMAX((dtype) 0.0, MYMIN((dtype) 1.0, za[jk_i]));

      //-------------------------------------------------------------------
      // Calculate liq/ice fractions (no longer a diagnostic relationship)
      //-------------------------------------------------------------------
      zli = zqx[0] + zqx[1];
      if (zli > (*yrecldp).rlmin) {
        zliqfrac = zqx[0] / zli;
        zicefrac = (dtype) 1.0 - zliqfrac;
      } else {
        zliqfrac = (dtype) 0.0;
        zicefrac = (dtype) 0.0;
      }

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
      //ZTRPAUS = 0.1_JPRB
      //ZPAPHD = 1.0_JPRB / PAPH(JL, KLEV + 1, IBL)
      //DO JK=1,KLEV - 1
      //  ZSIG = PAP(JL, JK, IBL)*ZPAPHD
      //  IF (ZSIG > 0.1_JPRB .and. ZSIG < 0.4_JPRB .and. ZTP1(JK_I) > ZTP1(JL, JK + 1, IBL)) THEN
      //    ZTRPAUS = ZSIG
      //  END IF
      //END DO

      //-----------------------------
      // Reset single level variables
      //-----------------------------

      //ZANEWM1 = 0.0_JPRB
      //ZDA = 0.0_JPRB
      //ZCOVPCLR = 0.0_JPRB
      //ZCOVPMAX = 0.0_JPRB
      //ZCOVPTOT = 0.0_JPRB
      //ZCLDTOPDIST = 0.0_JPRB

      //######################################################################
      //           3.       *** PHYSICS ***
      //######################################################################


      //----------------------------------------------------------------------
      //                       START OF VERTICAL LOOP
      //----------------------------------------------------------------------

      if ((*yrecldp).ncldtop <= jk + 1 && jk + 1 <= klev) {

        //----------------------------------------------------------------------
        // 3.0 INITIALIZE VARIABLES
        //----------------------------------------------------------------------

        //---------------------------------
        // First guess microphysics
        //---------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zqxfg[jm] = zqx[jm];
        }

        //---------------------------------
        // Set KLON arrays to zero
        //---------------------------------

        zlicld = (dtype) 0.0;
        zrainaut = (dtype) 0.0;            // currently needed for diags
        zrainacc = (dtype) 0.0;            // currently needed for diags
        zsnowaut = (dtype) 0.0;            // needed
        zldefr = (dtype) 0.0;
        zacust = (dtype) 0.0;            // set later when needed
        zqpretot = (dtype) 0.0;
        zlfinalsum = (dtype) 0.0;

        // Required for first guess call
        zlcond1 = (dtype) 0.0;
        zlcond2 = (dtype) 0.0;
        zsupsat = (dtype) 0.0;
        zlevapl = (dtype) 0.0;
        zlevapi = (dtype) 0.0;

        //-------------------------------------
        // solvers for cloud fraction
        //-------------------------------------
        zsolab = (dtype) 0.0;
        zsolac = (dtype) 0.0;

        zicetot = (dtype) 0.0;

        //------------------------------------------
        // reset matrix so missing pathways are set
        //------------------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            zsolqb[jn + 5*jm] = (dtype) 0.0;
            zsolqa[jn + 5*jm] = (dtype) 0.0;
          }
        }

        //----------------------------------
        // reset new microphysics variables
        //----------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zfallsrce[jm] = (dtype) 0.0;
          zfallsink[jm] = (dtype) 0.0;
          zconvsrce[jm] = (dtype) 0.0;
          zconvsink[jm] = (dtype) 0.0;
          zpsupsatsrce[jm] = (dtype) 0.0;
          zratio[jm] = (dtype) 0.0;
        }


        //-------------------------
        // derived variables needed
        //-------------------------

        zdp = paph[jl + klon*(1 + jk + (klev + 1)*ibl)] - paph[jl + klon*(jk + (klev +
          1)*ibl)];            // dp
        zgdp = rg / zdp;            // g/dp
        zrho = pap[jl + klon*(jk + klev*ibl)] / (rd*ztp1[jk_i]);            // p/RT air density

        zdtgdp = ptsphy*zgdp;            // dt g/dp
        zrdtgdp = zdp*((dtype) 1.0 / (ptsphy*rg));            // 1/(dt g/dp)

        if (jk + 1 > 1) {
          zdtgdpf = ptsphy*rg / (pap[jl + klon*(jk + klev*ibl)] - pap[jl + klon*(-1 +
            jk + klev*ibl)]);
        }

        //------------------------------------
        // Calculate dqs/dT correction factor
        //------------------------------------
        // Reminder: retv=rv/rd-1

        // liquid
        zfacw = r5les / (MYPOW((ztp1[jk_i] - r4les), 2));
        zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zfoeeliqt);
        zdqsliqdt = zfacw*zcor*zqsliq;
        zcorqsliq = (dtype) 1.0 + ralvdcp*zdqsliqdt;

        // ice
        zfaci = r5ies / (MYPOW((ztp1[jk_i] - r4ies), 2));
        zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zfoeew);
        zdqsicedt = zfaci*zcor*zqsice;
        zcorqsice = (dtype) 1.0 + ralsdcp*zdqsicedt;

        // diagnostic mixed
        zalfaw = zfoealfa;
        zalfawm = zalfaw;
        zfac = zalfaw*zfacw + ((dtype) 1.0 - zalfaw)*zfaci;
        zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zfoeewmt);
        zdqsmixdt = zfac*zcor*zqsmix;
        zcorqsmix = (dtype) 1.0 + ((dtype)((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*
          ralvdcp + (1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zdqsmixdt;

        // evaporation/sublimation limits
        zevaplimmix = MYMAX((zqsmix - zqx[4]) / zcorqsmix, (dtype) 0.0);
        zevaplimliq = MYMAX((zqsliq - zqx[4]) / zcorqsliq, (dtype) 0.0);
        zevaplimice = MYMAX((zqsice - zqx[4]) / zcorqsice, (dtype) 0.0);

        //--------------------------------
        // in-cloud consensate amount
        //--------------------------------
        ztmpa = (dtype) 1.0 / MYMAX(za[jk_i], zepsec);
        zliqcld = zqx[0]*ztmpa;
        zicecld = zqx[1]*ztmpa;
        zlicld = zliqcld + zicecld;


        //------------------------------------------------
        // Evaporate very small amounts of liquid and ice
        //------------------------------------------------

        if (zqx[0] < (*yrecldp).rlmin) {
          zsolqa[4 + 5*(0)] = zqx[0];
          zsolqa[0 + 5*(4)] = -zqx[0];
        }

        if (zqx[1] < (*yrecldp).rlmin) {
          zsolqa[4 + 5*(1)] = zqx[1];
          zsolqa[1 + 5*(4)] = -zqx[1];
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

        //-----------------------------------
        // 3.1.1 Supersaturation limit (from Koop)
        //-----------------------------------
        // Needs to be set for all temperatures
        zfokoop = ((dtype)(MYMIN(rkoop1 - rkoop2*ztp1[jk_i], (dtype)(r2es*MYEXP((r3les*(ztp1[jk_i] - rtt))/
	  (ztp1[jk_i] - r4les)))*1.0/(dtype)(r2es*MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))));

        if (ztp1[jk_i] >= rtt || (*yrecldp).nssopt == 0) {
          zfac = (dtype) 1.0;
          zfaci = (dtype) 1.0;
        } else {
          zfac = za[jk_i] + zfokoop*((dtype) 1.0 - za[jk_i]);
          zfaci = ptsphy / (*yrecldp).rkooptau;
        }

        //-------------------------------------------------------------------
        // 3.1.2 Calculate supersaturation wrt Koop including dqs/dT
        //       correction factor
        // [#Note: QSICE or QSLIQ]
        //-------------------------------------------------------------------

        // Calculate supersaturation to add to cloud
        if (za[jk_i] > (dtype) 1.0 - (*yrecldp).ramin) {
          zsupsat = MYMAX((zqx[4] - zfac*zqsice) / zcorqsice, (dtype) 0.0);
        } else {
          // Calculate environmental humidity supersaturation
          zqp1env = (zqx[4] - za[jk_i]*zqsice) / MYMAX((dtype) 1.0 - za[jk_i], zepsilon);
          //& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
          zsupsat = MYMAX(((dtype) 1.0 - za[jk_i])*(zqp1env - zfac*zqsice) /
            zcorqsice, (dtype) 0.0);
        }

        //-------------------------------------------------------------------
        // Here the supersaturation is turned into liquid water
        // However, if the temperature is below the threshold for homogeneous
        // freezing then the supersaturation is turned instantly to ice.
        //--------------------------------------------------------------------

        if (zsupsat > zepsec) {

          if (ztp1[jk_i] > (*yrecldp).rthomo) {
            // Turn supersaturation into liquid water
            zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] + zsupsat;
            zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] - zsupsat;
            // Include liquid in first guess
            zqxfg[0] = zqxfg[0] + zsupsat;
          } else {
            // Turn supersaturation into ice water
            zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] + zsupsat;
            zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] - zsupsat;
            // Add ice to first guess for deposition term
            zqxfg[1] = zqxfg[1] + zsupsat;
          }

          // Increase cloud amount using RKOOPTAU timescale
          zsolac = ((dtype) 1.0 - za[jk_i])*zfaci;

        }

        //-------------------------------------------------------
        // 3.1.3 Include supersaturation from previous timestep
        // (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
        //-------------------------------------------------------
        if (psupsat[jl + klon*(jk + klev*ibl)] > zepsec) {
          if (ztp1[jk_i] > (*yrecldp).rthomo) {
            // Turn supersaturation into liquid water
            zsolqa[0 + 5*(0)] =
              zsolqa[0 + 5*(0)] + psupsat[jl + klon*(jk + klev*ibl)];
            zpsupsatsrce[0] = psupsat[jl + klon*(jk + klev*ibl)];
            // Add liquid to first guess for deposition term
            zqxfg[0] = zqxfg[0] + psupsat[jl + klon*(jk + klev*ibl)];
            // Store cloud budget diagnostics if required
          } else {
            // Turn supersaturation into ice water
            zsolqa[1 + 5*(1)] =
              zsolqa[1 + 5*(1)] + psupsat[jl + klon*(jk + klev*ibl)];
            zpsupsatsrce[1] = psupsat[jl + klon*(jk + klev*ibl)];
            // Add ice to first guess for deposition term
            zqxfg[1] = zqxfg[1] + psupsat[jl + klon*(jk + klev*ibl)];
            // Store cloud budget diagnostics if required
          }

          // Increase cloud amount using RKOOPTAU timescale
          zsolac = ((dtype) 1.0 - za[jk_i])*zfaci;
          // Store cloud budget diagnostics if required
        }

        // on JL

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
        if (jk + 1 < klev && jk + 1 >= (*yrecldp).ncldtop) {


          plude[jl + klon*(jk + klev*ibl)] = plude[jl + klon*(jk + klev*ibl)]*zdtgdp;

          if (/*ldcum[jl + klon*ibl] &&*/ plude[jl + klon*(jk + klev*ibl)] > (*yrecldp)
            .rlmin && plu[jl + klon*(1 + jk + klev*ibl)] > zepsec) {

            zsolac = zsolac + plude[jl + klon*(jk + klev*ibl)] / plu[jl + klon*(1 + jk
              + klev*ibl)];
            // *diagnostic temperature split*
            zalfaw = zfoealfa;
            zconvsrce[0] = zalfaw*plude[jl + klon*(jk + klev*ibl)];
            zconvsrce[1] =
              ((dtype) 1.0 - zalfaw)*plude[jl + klon*(jk + klev*ibl)];
            zsolqa[0 + 5*(0)] =
              zsolqa[0 + 5*(0)] + zconvsrce[0];
            zsolqa[1 + 5*(1)] =
              zsolqa[1 + 5*(1)] + zconvsrce[1];

          } else {

            plude[jl + klon*(jk + klev*ibl)] = (dtype) 0.0;

          }
          // *convective snow detrainment source
          //if (ldcum[jl + klon*ibl]) {
            zsolqa[3 + 5*(3)] =
              zsolqa[3 + 5*(3)] + psnde[jl + klon*(jk + klev*ibl)]*zdtgdp;
          //}


        }
        // JK<KLEV

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
        if (jk + 1 > (*yrecldp).ncldtop) {

          zmf = MYMAX((dtype) 0.0, (pmfu[jl + klon*(jk + klev*ibl)] + pmfd[jl +
            klon*(jk + klev*ibl)])*zdtgdp);
          zacust = zmf*zanewm1;

          for (jm = 0; jm <= 5 + -1; jm += 1) {
            if (!llfall[jm] && iphase[jm] > 0) {
              zlcust[jm] = zmf*zqxnm1[jm];
              // record total flux for enthalpy budget:
              zconvsrce[jm] = zconvsrce[jm] + zlcust[jm];
            }
          }

          // Now have to work out how much liquid evaporates at arrival point
          // since there is no prognostic memory for in-cloud humidity, i.e.
          // we always assume cloud is saturated.

          zdtdp = zrdcp*(dtype) 0.5*(ztp1[jk_im1] + ztp1[jk_i]) / paph[jl +
            klon*(jk + (klev + 1)*ibl)];
          zdtforc = zdtdp*(pap[jl + klon*(jk + klev*ibl)] - pap[jl + klon*(-1 + jk +
            klev*ibl)]);
          //[#Note: Diagnostic mixed phase should be replaced below]
          zdqs = zanewm1*zdtforc*zdqsmixdt;

          for (jm = 0; jm <= 5 + -1; jm += 1) {
            if (!llfall[jm] && iphase[jm] > 0) {
              zlfinal = MYMAX((dtype) 0.0, zlcust[jm] - zdqs);                  //lim to zero
              // no supersaturation allowed incloud ---V
              zevap = MYMIN((zlcust[jm] - zlfinal), zevaplimmix);
              //          ZEVAP=0.0_JPRB
              zlfinal = zlcust[jm] - zevap;
              zlfinalsum = zlfinalsum + zlfinal;                  // sum

              zsolqa[jm + 5*jm] = zsolqa[jm + 5*jm] + zlcust[jm];                  // whole sum
              zsolqa[4 + 5*jm] = zsolqa[4 + 5*jm] + zevap;
              zsolqa[jm + 5*(4)] = zsolqa[jm + 5*(4)] - zevap;
            }
          }

          //  Reset the cloud contribution if no cloud water survives to this level:
          if (zlfinalsum < zepsec) {
            zacust = (dtype) 0.0;
          }
          zsolac = zsolac + zacust;

        }
        // on  JK>NCLDTOP

        //---------------------------------------------------------------------
        // Subsidence sink of cloud to the layer below
        // (Implicit - re. CFL limit on convective mass flux)
        //---------------------------------------------------------------------


        if (jk + 1 < klev) {

          zmfdn = MYMAX((dtype) 0.0, (pmfu[jl + klon*(1 + jk + klev*ibl)] + pmfd[jl +
            klon*(1 + jk + klev*ibl)])*zdtgdp);

          zsolab = zsolab + zmfdn;
          zsolqb[0 + 5*(0)] = zsolqb[0 + 5*(0)] + zmfdn;
          zsolqb[1 + 5*(1)] = zsolqb[1 + 5*(1)] + zmfdn;

          // Record sink for cloud budget and enthalpy budget diagnostics
          zconvsink[0] = zmfdn;
          zconvsink[1] = zmfdn;

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
        zldifdt = (*yrecldp).rcldiff*ptsphy;            //original version
        //Increase by factor of 5 for convective points
        if (ktype[jl + klon*ibl] > 0 && plude[jl + klon*(jk + klev*ibl)] > zepsec) {
          zldifdt = (*yrecldp).rcldiff_convi*zldifdt;
        }

        // At the moment, works on mixed RH profile and partitioned ice/liq fraction
        // so that it is similar to previous scheme
        // Should apply RHw for liquid cloud and RHi for ice cloud separately
        if (zli > zepsec) {
          // Calculate environmental humidity
          //      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
          //    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
          //      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
          ze = zldifdt*MYMAX(zqsmix - zqx[4], (dtype) 0.0);
          zleros = za[jk_i]*ze;
          zleros = MYMIN(zleros, zevaplimmix);
          zleros = MYMIN(zleros, zli);
          zaeros = zleros / zlicld;              //if linear term

          // Erosion is -ve LINEAR in L,A
          zsolac = zsolac - zaeros;              //linear

          zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] + zliqfrac*zleros;
          zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] - zliqfrac*zleros;
          zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] + zicefrac*zleros;
          zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] - zicefrac*zleros;

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

        zdtdp = zrdcp*ztp1[jk_i] / pap[jl + klon*(jk + klev*ibl)];
        zdpmxdt = zdp*zqtmst;
        zmfdn = (dtype) 0.0;
        if (jk + 1 < klev) {
          zmfdn =
            pmfu[jl + klon*(1 + jk + klev*ibl)] + pmfd[jl + klon*(1 + jk + klev*ibl)];
        }
        zwtot = pvervel[jl + klon*(jk + klev*ibl)] + (dtype) 0.5*rg*(pmfu[jl +
          klon*(jk + klev*ibl)] + pmfd[jl + klon*(jk + klev*ibl)] + zmfdn);
        zwtot = MYMIN(zdpmxdt, MYMAX(-zdpmxdt, zwtot));
        zzzdt = phrsw[jl + klon*(jk + klev*ibl)] + phrlw[jl + klon*(jk + klev*ibl)];
        zdtdiab =
          MYMIN(zdpmxdt*zdtdp, MYMAX(-zdpmxdt*zdtdp, zzzdt))*ptsphy + ralfdcp*zldefr;
        // Note: ZLDEFR should be set to the difference between the mixed phase functions
        // in the convection and cloud scheme, but this is not calculated, so is zero and
        // the functions must be the same
        zdtforc = zdtdp*zwtot*ptsphy + zdtdiab;
        zqold = zqsmix;
        ztold = ztp1[jk_i];
        ztp1[jk_i] = ztp1[jk_i] + zdtforc;
        ztp1[jk_i] = MYMAX(ztp1[jk_i], (dtype) 160.0);
        llflag = true;

        // Formerly a call to CUADJTQ(..., ICALL=5)
        zqp = (dtype) 1.0 / pap[jl + klon*(jk + klev*ibl)];
        zqsat = ((dtype)(r2es*((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*
          MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les)) + 
	  (1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*
	  MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies)))))*zqp;
        zqsat = MYMIN((dtype) 0.5, zqsat);
        zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zqsat);
        zqsat = zqsat*zcor;
        zcond = (zqsmix - zqsat) / ((dtype) 1.0 + zqsat*zcor*((dtype)(((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*
          (1.0/MYPOW(ztp1[jk_i] - r4les, 2)) + ((1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*
	  (1.0/MYPOW(ztp1[jk_i] - r4ies, 2)))));
        ztp1[jk_i] = ztp1[jk_i] + ((dtype)((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + 
          (1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zcond;
        zqsmix = zqsmix - zcond;
        zqsat = ((dtype)(r2es*((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*
          MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les)) + 
          (1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*
          MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies)))))*zqp;
        zqsat = MYMIN((dtype) 0.5, zqsat);
        zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zqsat);
        zqsat = zqsat*zcor;
        zcond1 = (zqsmix - zqsat) / ((dtype) 1.0 + zqsat*zcor*((dtype)(((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*
          (1.0/MYPOW(ztp1[jk_i] - r4les, 2)) + ((1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*
	  (1.0/MYPOW(ztp1[jk_i] - r4ies, 2)))));
        ztp1[jk_i] = ztp1[jk_i] + ((dtype)((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + 
          (1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zcond1;
        zqsmix = zqsmix - zcond1;

        zdqs = zqsmix - zqold;
        zqsmix = zqold;
        ztp1[jk_i] = ztold;

        //----------------------------------------------------------------------
        // 3.4a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
        // ----------------------------------------------------------------------
        // Erosion term is LINEAR in L
        // Changed to be uniform distribution in cloud region


        // Previous function based on DELTA DISTRIBUTION in cloud:
        if (zdqs > (dtype) 0.0) {
          //    If subsidence evaporation term is turned off, then need to use updated
          //    liquid and cloud here?
          //    ZLEVAP = MAX(ZA(JL,JK)+ZACUST(JL),1.0_JPRB)*MIN(ZDQS(JL),ZLICLD(JL)+ZLFINALSUM(JL))
          zlevap = za[jk_i]*MYMIN(zdqs, zlicld);
          zlevap = MYMIN(zlevap, zevaplimmix);
          zlevap = MYMIN(zlevap, MYMAX(zqsmix - zqx[4], (dtype) 0.0));

          // For first guess call
          zlevapl = zliqfrac*zlevap;
          zlevapi = zicefrac*zlevap;

          zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] + zliqfrac*zlevap;
          zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] - zliqfrac*zlevap;

          zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] + zicefrac*zlevap;
          zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] - zicefrac*zlevap;

        }


        //----------------------------------------------------------------------
        // 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
        //----------------------------------------------------------------------
        // (1) Increase of cloud water in existing clouds
        if (za[jk_i] > zepsec && zdqs <= -(*yrecldp).rlmin) {

          zlcond1 = MYMAX(-zdqs, (dtype) 0.0);              //new limiter

          //old limiter (significantly improves upper tropospheric humidity rms)
          if (za[jk_i] > (dtype) 0.99) {
            zcor = (dtype) 1.0 / ((dtype) 1.0 - retv*zqsmix);
            zcdmax = (zqx[4] - zqsmix) / ((dtype) 1.0 + zcor*zqsmix*((dtype)(((dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*
              (1.0/MYPOW(ztp1[jk_i] - r4les, 2)) + ((1.0 - (dtype)(MYMIN(1.0, MYPOW((MYMAX(rtice, MYMIN(rtwat, ztp1[jk_i])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*
	      (1.0/MYPOW(ztp1[jk_i] - r4ies, 2)))));
          } else {
            zcdmax = (zqx[4] - za[jk_i]*zqsmix) / za[jk_i];
          }
          zlcond1 = MYMAX(MYMIN(zlcond1, zcdmax), (dtype) 0.0);
          // end old limiter

          zlcond1 = za[jk_i]*zlcond1;
          if (zlcond1 < (*yrecldp).rlmin) {
            zlcond1 = (dtype) 0.0;
          }

          //-------------------------------------------------------------------------
          // All increase goes into liquid unless so cold cloud homogeneously freezes
          // Include new liquid formation in first guess value, otherwise liquid
          // remains at cold temperatures until next timestep.
          //-------------------------------------------------------------------------
          if (ztp1[jk_i] > (*yrecldp).rthomo) {
            zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] + zlcond1;
            zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] - zlcond1;
            zqxfg[0] = zqxfg[0] + zlcond1;
          } else {
            zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] + zlcond1;
            zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] - zlcond1;
            zqxfg[1] = zqxfg[1] + zlcond1;
          }
        }

        // (2) Generation of new clouds (da/dt>0)


        if (zdqs <= -(*yrecldp).rlmin && za[jk_i] < (dtype) 1.0 - zepsec) {

          //---------------------------
          // Critical relative humidity
          //---------------------------
          zrhc = (*yrecldp).ramid;
          zsigk =
            pap[jl + klon*(jk + klev*ibl)] / paph[jl + klon*(klev + (klev + 1)*ibl)];
          // Increase RHcrit to 1.0 towards the surface (eta>0.8)
          if (zsigk > (dtype) 0.8) {
            zrhc = (*yrecldp).ramid + ((dtype) 1.0 - (*yrecldp).ramid)*(MYPOW(((zsigk -
              (dtype) 0.8) / (dtype) 0.2), 2));
          }

          // Commented out for CY37R1 to reduce humidity in high trop and strat
          //      ! Increase RHcrit to 1.0 towards the tropopause (trop-0.2) and above
          //      ZBOTT=ZTRPAUS(JL)+0.2_JPRB
          //      IF(ZSIGK < ZBOTT) THEN
          //        ZRHC=RAMID+(1.0_JPRB-RAMID)*MIN(((ZBOTT-ZSIGK)/0.2_JPRB)**2,1.0_JPRB)
          //      ENDIF

          //---------------------------
          // Supersaturation options
          //---------------------------
          if ((*yrecldp).nssopt == 0) {
            // No scheme
            zqe = (zqx[4] - za[jk_i]*zqsice) / MYMAX(zepsec, (dtype) 1.0 -
              za[jk_i]);
            zqe = MYMAX((dtype) 0.0, zqe);
          } else if ((*yrecldp).nssopt == 1) {
            // Tompkins
            zqe = (zqx[4] - za[jk_i]*zqsice) / MYMAX(zepsec, (dtype) 1.0 -
              za[jk_i]);
            zqe = MYMAX((dtype) 0.0, zqe);
          } else if ((*yrecldp).nssopt == 2) {
            // Lohmann and Karcher
            zqe = zqx[4];
          } else if ((*yrecldp).nssopt == 3) {
            // Gierens
            zqe = zqx[4] + zli;
          }

          if (ztp1[jk_i] >= rtt || (*yrecldp).nssopt == 0) {
            // No ice supersaturation allowed
            zfac = (dtype) 1.0;
          } else {
            // Ice supersaturation
            zfac = zfokoop;
          }

          if (zqe >= zrhc*zqsice*zfac && zqe < zqsice*zfac) {
            // note: not **2 on 1-a term if ZQE is used.
            // Added correction term ZFAC to numerator 15/03/2010
            zacond = -((dtype) 1.0 - za[jk_i])*zfac*zdqs / MYMAX((dtype)
              2.0*(zfac*zqsice - zqe), zepsec);

            zacond = MYMIN(zacond, (dtype) 1.0 - za[jk_i]);                //PUT THE LIMITER BACK

            // Linear term:
            // Added correction term ZFAC 15/03/2010
            zlcond2 = -zfac*zdqs*(dtype) 0.5*zacond;                //mine linear

            // new limiter formulation
            zzdl = (dtype) 2.0*(zfac*zqsice - zqe) / MYMAX(zepsec, (dtype) 1.0 - za[jk_i]);
            // Added correction term ZFAC 15/03/2010
            if (zfac*zdqs < -zzdl) {
              // ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
              zlcondlim =
                (za[jk_i] - (dtype) 1.0)*zfac*zdqs - zfac*zqsice + zqx[4];
              zlcond2 = MYMIN(zlcond2, zlcondlim);
            }
            zlcond2 = MYMAX(zlcond2, (dtype) 0.0);

            if (zlcond2 < (*yrecldp).rlmin || ((dtype) 1.0 - za[jk_i]) < zepsec
              ) {
              zlcond2 = (dtype) 0.0;
              zacond = (dtype) 0.0;
            }
            if (zlcond2 == (dtype) 0.0) {
              zacond = (dtype) 0.0;
            }

            // Large-scale generation is LINEAR in A and LINEAR in L
            zsolac = zsolac + zacond;                //linear

            //------------------------------------------------------------------------
            // All increase goes into liquid unless so cold cloud homogeneously freezes
            // Include new liquid formation in first guess value, otherwise liquid
            // remains at cold temperatures until next timestep.
            //------------------------------------------------------------------------
            if (ztp1[jk_i] > (*yrecldp).rthomo) {
              zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] + zlcond2;
              zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] - zlcond2;
              zqxfg[0] = zqxfg[0] + zlcond2;
            } else {
              // homogeneous freezing
              zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] + zlcond2;
              zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] - zlcond2;
              zqxfg[1] = zqxfg[1] + zlcond2;
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
        if (idepice == 1) {


          //--------------------------------------------------------------
          // Calculate distance from cloud top
          // defined by cloudy layer below a layer with cloud frac <0.01
          // ZDZ = ZDP(JL)/(ZRHO(JL)*rg)
          //--------------------------------------------------------------

          if (za[jk_im1] < (*yrecldp).rcldtopcf && za[jk_i] >= (*yrecldp)
            .rcldtopcf) {
            zcldtopdist = (dtype) 0.0;
          } else {
            zcldtopdist = zcldtopdist + zdp / (zrho*rg);
          }

          //--------------------------------------------------------------
          // only treat depositional growth if liquid present. due to fact
          // that can not model ice growth from vapour without additional
          // in-cloud water vapour variable
          //--------------------------------------------------------------
          if (ztp1[jk_i] < rtt && zqxfg[0] > (*yrecldp).rlmin) {
            // T<273K

            zvpice = (((dtype)(r2es*MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))*rv) / rd;
            zvpliq = zvpice*zfokoop;
            zicenuclei = (dtype) 1000.0*MYEXP((dtype) 12.96*(zvpliq - zvpice) / zvpliq
              - (dtype) 0.639);

            //------------------------------------------------
            //   2.4e-2 is conductivity of air
            //   8.8 = 700**1/3 = density of ice to the third
            //------------------------------------------------
            zadd = rlstt*(rlstt / (rv*ztp1[jk_i]) - (dtype) 1.0) / ((dtype)
              2.4E-2*ztp1[jk_i]);
            zbdd = rv*ztp1[jk_i]*pap[jl + klon*(jk + klev*ibl)] / ((dtype)
              2.21*zvpice);
            zcvds = (dtype) 7.8*(MYPOW((zicenuclei / zrho), (dtype) 0.666))*(zvpliq -
              zvpice) / ((dtype) 8.87*(zadd + zbdd)*zvpice);

            //-----------------------------------------------------
            // RICEINIT=1.E-12_JPRB is initial mass of ice particle
            //-----------------------------------------------------
            zice0 = MYMAX(zicecld, zicenuclei*(*yrecldp).riceinit / zrho);

            //------------------
            // new value of ice:
            //------------------
            zinew = MYPOW(((dtype) 0.666*zcvds*ptsphy + (MYPOW(zice0, (dtype) 0.666))),
              (dtype) 1.5);

            //---------------------------
            // grid-mean deposition rate:
            //---------------------------
            zdepos = MYMAX(za[jk_i]*(zinew - zice0), (dtype) 0.0);

            //--------------------------------------------------------------------
            // Limit deposition to liquid water amount
            // If liquid is all frozen, ice would use up reservoir of water
            // vapour in excess of ice saturation mixing ratio - However this
            // can not be represented without a in-cloud humidity variable. Using
            // the grid-mean humidity would imply a large artificial horizontal
            // flux from the clear sky to the cloudy area. We thus rely on the
            // supersaturation check to clean up any remaining supersaturation
            //--------------------------------------------------------------------
            zdepos = MYMIN(zdepos, zqxfg[0]);                // limit to liquid water amount

            //--------------------------------------------------------------------
            // At top of cloud, reduce deposition rate near cloud top to account for
            // small scale turbulent processes, limited ice nucleation and ice fallout
            //--------------------------------------------------------------------
            //      ZDEPOS = ZDEPOS*MIN(RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH,1.0_JPRB)
            // Change to include dependence on ice nuclei concentration
            // to increase deposition rate with decreasing temperatures
            zinfactor = MYMIN(zicenuclei / (dtype) 15000., (dtype) 1.0);
            zdepos = zdepos*MYMIN(zinfactor + ((dtype) 1.0 - zinfactor)*((*yrecldp)
              .rdepliqrefrate + zcldtopdist / (*yrecldp).rdepliqrefdepth), (dtype) 1.0
              );

            //--------------
            // add to matrix
            //--------------
            zsolqa[1 + 5*(0)] = zsolqa[1 + 5*(0)] + zdepos;
            zsolqa[0 + 5*(1)] = zsolqa[0 + 5*(1)] - zdepos;
            zqxfg[1] = zqxfg[1] + zdepos;
            zqxfg[0] = zqxfg[0] - zdepos;

          }

          //--------------------------------------------------------
          //-
          //- Ice deposition assuming ice PSD
          //-
          //--------------------------------------------------------
        } else if (idepice == 2) {


          //--------------------------------------------------------------
          // Calculate distance from cloud top
          // defined by cloudy layer below a layer with cloud frac <0.01
          // ZDZ = ZDP(JL)/(ZRHO(JL)*rg)
          //--------------------------------------------------------------

          if (za[jk_im1] < (*yrecldp).rcldtopcf && za[jk_i] >= (*yrecldp)
            .rcldtopcf) {
            zcldtopdist = (dtype) 0.0;
          } else {
            zcldtopdist = zcldtopdist + zdp / (zrho*rg);
          }

          //--------------------------------------------------------------
          // only treat depositional growth if liquid present. due to fact
          // that can not model ice growth from vapour without additional
          // in-cloud water vapour variable
          //--------------------------------------------------------------
          if (ztp1[jk_i] < rtt && zqxfg[0] > (*yrecldp).rlmin) {
            // T<273K

            zvpice = (((dtype)(r2es*MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))*rv) / rd;
            zvpliq = zvpice*zfokoop;
            zicenuclei = (dtype) 1000.0*MYEXP((dtype) 12.96*(zvpliq - zvpice) / zvpliq
              - (dtype) 0.639);

            //-----------------------------------------------------
            // RICEINIT=1.E-12_JPRB is initial mass of ice particle
            //-----------------------------------------------------
            zice0 = MYMAX(zicecld, zicenuclei*(*yrecldp).riceinit / zrho);

            // Particle size distribution
            ztcg = (dtype) 1.0;
            zfacx1i = (dtype) 1.0;

            zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk_i] + 
              pap[jl + klon*(jk + klev*ibl)]*(*yrecldp).rcl_apb3*(MYPOW(ztp1[jk_i], (dtype) 3.));
            zcorrfac = MYPOW(((dtype) 1.0 / zrho), (dtype) 0.5);
            zcorrfac2 = (MYPOW((ztp1[jk_i] / (dtype) 273.0), (dtype) 1.5))
              *((dtype) 393.0 / (ztp1[jk_i] + (dtype) 120.0));

            zpr02 = zrho*zice0*(*yrecldp).rcl_const1i / (ztcg*zfacx1i);

            zterm1 = (zvpliq - zvpice)*(MYPOW(ztp1[jk_i], (dtype) 2.0))
              *zvpice*zcorrfac2*ztcg*(*yrecldp).rcl_const2i*zfacx1i /
              (zrho*zaplusb*zvpice);
            zterm2 = (dtype) 0.65*(*yrecldp).rcl_const6i*(MYPOW(zpr02, (*yrecldp)
              .rcl_const4i)) + (*yrecldp).rcl_const3i*(MYPOW(zcorrfac, (dtype) 0.5))
              *(MYPOW(zrho, (dtype) 0.5))*(MYPOW(zpr02, (*yrecldp).rcl_const5i)) /
              (MYPOW(zcorrfac2, (dtype) 0.5));

            zdepos = MYMAX(za[jk_i]*zterm1*zterm2*ptsphy, (dtype) 0.0);

            //--------------------------------------------------------------------
            // Limit deposition to liquid water amount
            // If liquid is all frozen, ice would use up reservoir of water
            // vapour in excess of ice saturation mixing ratio - However this
            // can not be represented without a in-cloud humidity variable. Using
            // the grid-mean humidity would imply a large artificial horizontal
            // flux from the clear sky to the cloudy area. We thus rely on the
            // supersaturation check to clean up any remaining supersaturation
            //--------------------------------------------------------------------
            zdepos = MYMIN(zdepos, zqxfg[0]);                // limit to liquid water amount

            //--------------------------------------------------------------------
            // At top of cloud, reduce deposition rate near cloud top to account for
            // small scale turbulent processes, limited ice nucleation and ice fallout
            //--------------------------------------------------------------------
            // Change to include dependence on ice nuclei concentration
            // to increase deposition rate with decreasing temperatures
            zinfactor = MYMIN(zicenuclei / (dtype) 15000., (dtype) 1.0);
            zdepos = zdepos*MYMIN(zinfactor + ((dtype) 1.0 - zinfactor)*((*yrecldp)
              .rdepliqrefrate + zcldtopdist / (*yrecldp).rdepliqrefdepth), (dtype) 1.0
              );

            //--------------
            // add to matrix
            //--------------
            zsolqa[1 + 5*(0)] = zsolqa[1 + 5*(0)] + zdepos;
            zsolqa[0 + 5*(1)] = zsolqa[0 + 5*(1)] - zdepos;
            zqxfg[1] = zqxfg[1] + zdepos;
            zqxfg[0] = zqxfg[0] - zdepos;
          }

        }
        // on IDEPICE

        //######################################################################
        //              4  *** PRECIPITATION PROCESSES ***
        //######################################################################

        //----------------------------------
        // revise in-cloud consensate amount
        //----------------------------------
        ztmpa = (dtype) 1.0 / MYMAX(za[jk_i], zepsec);
        zliqcld = zqxfg[0]*ztmpa;
        zicecld = zqxfg[1]*ztmpa;
        zlicld = zliqcld + zicecld;

        //----------------------------------------------------------------------
        // 4.2 SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
        //     now that rain, snow, graupel species are prognostic
        //     the precipitation flux can be defined directly level by level
        //     There is no vertical memory required from the flux variable
        //----------------------------------------------------------------------

        for (jm = 0; jm <= 5 + -1; jm += 1) {
          if (llfall[jm] || jm + 1 == 2) {
            //------------------------
            // source from layer above
            //------------------------
            if (jk + 1 > (*yrecldp).ncldtop) {
              zfallsrce[jm] = zpfplsx[jk_i + 2*jm]*zdtgdp;
              zsolqa[jm + 5*jm] = zsolqa[jm + 5*jm] + zfallsrce[jm];
              zqxfg[jm] = zqxfg[jm] + zfallsrce[jm];
              // use first guess precip----------V
              zqpretot = zqpretot + zqxfg[jm];
            }
            //-------------------------------------------------
            // sink to next layer, constant fall speed
            //-------------------------------------------------
            // if aerosol effect then override
            //  note that for T>233K this is the same as above.
            if ((*yrecldp).laericesed && jm + 1 == 2) {
              zre_ice = pre_ice[jl + klon*(jk + klev*ibl)];
              // The exponent value is from
              // Morrison et al. JAS 2005 Appendix
              zvqx[1] = (dtype) 0.002*(MYPOW(zre_ice, (dtype) 1.0));
            }
            zfall = zvqx[jm]*zrho;
            //-------------------------------------------------
            // modified by Heymsfield and Iaquinta JAS 2000
            //-------------------------------------------------
            // ZFALL = ZFALL*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
            //            &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))

            zfallsink[jm] = zdtgdp*zfall;
            // Cloud budget diagnostic stored at end as implicit
            // jl
          }
          // LLFALL
        }
        // jm

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
        if (zqpretot > zepsec) {
          zcovptot = (dtype) 1.0 - (((dtype) 1.0 - zcovptot)*((dtype) 1.0 -
            MYMAX(za[jk_i], za[jk_im1]))) / ((dtype) 1.0 - MYMIN(za[jk_im1], (dtype) 1.0 - (dtype) 1.E-06));              // here!!!
          zcovptot = MYMAX(zcovptot, (*yrecldp).rcovpmin);
          zcovpclr = MYMAX((dtype) 0.0, zcovptot - za[jk_i]);              // clear sky proportion
          zraincld = zqxfg[2] / zcovptot;
          zsnowcld = zqxfg[3] / zcovptot;
          zcovpmax = MYMAX(zcovptot, zcovpmax);
        } else {
          zraincld = (dtype) 0.0;
          zsnowcld = (dtype) 0.0;
          zcovptot = (dtype) 0.0;              // no flux - reset cover
          zcovpclr = (dtype) 0.0;              // reset clear sky proportion
          zcovpmax = (dtype) 0.0;              // reset max cover for ZZRH calc
        }

        //----------------------------------------------------------------------
        // 4.3a AUTOCONVERSION TO SNOW
        //----------------------------------------------------------------------

        if (ztp1[jk_i] <= rtt) {
          //-----------------------------------------------------
          //     Snow Autoconversion rate follow Lin et al. 1983
          //-----------------------------------------------------
          if (zicecld > zepsec) {

            zzco = ptsphy*(*yrecldp).rsnowlin1*MYEXP((*yrecldp).rsnowlin2*(ztp1[jk_i] - rtt));

            if ((*yrecldp).laericeauto) {
              zlcrit = picrit_aer[jl + klon*(jk + klev*ibl)];
              // 0.3 = N**0.333 with N=0.027
              zzco = zzco*(MYPOW(((*yrecldp).rnice / pnice[jl + klon*(jk + klev*ibl)]),
                (dtype) 0.333));
            } else {
              zlcrit = (*yrecldp).rlcritsnow;
            }

            zsnowaut = zzco*((dtype) 1.0 - MYEXP(-(MYPOW((zicecld / zlcrit), 2))));
            zsolqb[3 + 5*(1)] = zsolqb[3 + 5*(1)] + zsnowaut;

          }
        }

        //----------------------------------------------------------------------
        // 4.3b AUTOCONVERSION WARM CLOUDS
        //   Collection and accretion will require separate treatment
        //   but for now we keep this simple treatment
        //----------------------------------------------------------------------

        if (zliqcld > zepsec) {

          //--------------------------------------------------------
          //-
          //- Warm-rain process follow Sundqvist (1989)
          //-
          //--------------------------------------------------------
          if (iwarmrain == 1) {

            zzco = (*yrecldp).rkconv*ptsphy;

            if ((*yrecldp).laerliqautolsp) {
              zlcrit = plcrit_aer[jl + klon*(jk + klev*ibl)];
              // 0.3 = N**0.333 with N=125 cm-3
              zzco = zzco*(MYPOW(((*yrecldp).rccn / pccn[jl + klon*(jk + klev*ibl)]),
                (dtype) 0.333));
            } else {
              // Modify autoconversion threshold dependent on:
              //  land (polluted, high CCN, smaller droplets, higher threshold)
              //  sea  (clean, low CCN, larger droplets, lower threshold)
              if (plsm[jl + klon*ibl] > (dtype) 0.5) {
                zlcrit = (*yrecldp).rclcrit_land;                    // land
              } else {
                zlcrit = (*yrecldp).rclcrit_sea;                    // ocean
              }
            }

            //------------------------------------------------------------------
            // Parameters for cloud collection by rain and snow.
            // Note that with new prognostic variable it is now possible
            // to REPLACE this with an explicit collection parametrization
            //------------------------------------------------------------------
            zprecip = (zpfplsx[jk_i + 2*(3)] + zpfplsx[jk_i + 2*(2)
              ]) / MYMAX(zepsec, zcovptot);
            zcfpr = (dtype) 1.0 + (*yrecldp).rprc1*sqrt(MYMAX(zprecip, (dtype) 0.0));
            //      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
            //       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))

            if ((*yrecldp).laerliqcoll) {
              // 5.0 = N**0.333 with N=125 cm-3
              zcfpr = zcfpr*(MYPOW(((*yrecldp).rccn / pccn[jl + klon*(jk + klev*ibl)]),
                (dtype) 0.333));
            }

            zzco = zzco*zcfpr;
            zlcrit = zlcrit / MYMAX(zcfpr, zepsec);

            if (zliqcld / zlcrit < (dtype) 20.0) {
              // Security for exp for some compilers
              zrainaut = zzco*((dtype) 1.0 - MYEXP(-(MYPOW((zliqcld / zlcrit), 2))));
            } else {
              zrainaut = zzco;
            }

            // rain freezes instantly
            if (ztp1[jk_i] <= rtt) {
              zsolqb[3 + 5*(0)] = zsolqb[3 + 5*(0)] + zrainaut;
            } else {
              zsolqb[2 + 5*(0)] = zsolqb[2 + 5*(0)] + zrainaut;
            }

            //--------------------------------------------------------
            //-
            //- Warm-rain process follow Khairoutdinov and Kogan (2000)
            //-
            //--------------------------------------------------------
          } else if (iwarmrain == 2) {

            if (plsm[jl + klon*ibl] > (dtype) 0.5) {
              // land
              zconst = (*yrecldp).rcl_kk_cloud_num_land;
              zlcrit = (*yrecldp).rclcrit_land;
            } else {
              // ocean
              zconst = (*yrecldp).rcl_kk_cloud_num_sea;
              zlcrit = (*yrecldp).rclcrit_sea;
            }

            if (zliqcld > zlcrit) {

              zrainaut = (dtype) 1.5*za[jk_i]*ptsphy*(*yrecldp)
                .rcl_kkaau*(MYPOW(zliqcld, (*yrecldp).rcl_kkbauq))*(MYPOW(zconst, (*yrecldp
                ).rcl_kkbaun));

              zrainaut = MYMIN(zrainaut, zqxfg[0]);
              if (zrainaut < zepsec) {
                zrainaut = (dtype) 0.0;
              }

              zrainacc = (dtype) 2.0*za[jk_i]*ptsphy*(*yrecldp)
                .rcl_kkaac*(MYPOW((zliqcld*zraincld), (*yrecldp).rcl_kkbac));

              zrainacc = MYMIN(zrainacc, zqxfg[0]);
              if (zrainacc < zepsec) {
                zrainacc = (dtype) 0.0;
              }

            } else {
              zrainaut = (dtype) 0.0;
              zrainacc = (dtype) 0.0;
            }

            // If temperature < 0, then autoconversion produces snow rather than rain
            // Explicit
            if (ztp1[jk_i] <= rtt) {
              zsolqa[3 + 5*(0)] = zsolqa[3 + 5*(0)] + zrainaut;
              zsolqa[3 + 5*(0)] = zsolqa[3 + 5*(0)] + zrainacc;
              zsolqa[0 + 5*(3)] = zsolqa[0 + 5*(3)] - zrainaut;
              zsolqa[0 + 5*(3)] = zsolqa[0 + 5*(3)] - zrainacc;
            } else {
              zsolqa[2 + 5*(0)] = zsolqa[2 + 5*(0)] + zrainaut;
              zsolqa[2 + 5*(0)] = zsolqa[2 + 5*(0)] + zrainacc;
              zsolqa[0 + 5*(2)] = zsolqa[0 + 5*(2)] - zrainaut;
              zsolqa[0 + 5*(2)] = zsolqa[0 + 5*(2)] - zrainacc;
            }

          }
          // on IWARMRAIN

        }
        // on ZLIQCLD > ZEPSEC


        //----------------------------------------------------------------------
        // RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
        //      only active if T<0degC and supercooled liquid water is present
        //      AND if not Sundquist autoconversion (as this includes riming)
        //----------------------------------------------------------------------
        if (iwarmrain > 1) {

          if (ztp1[jk_i] <= rtt && zliqcld > zepsec) {

            // Fallspeed air density correction
            zfallcorr = MYPOW(((*yrecldp).rdensref / zrho), (dtype) 0.4);

            //------------------------------------------------------------------
            // Riming of snow by cloud water - implicit in lwc
            //------------------------------------------------------------------
            if (zsnowcld > zepsec && zcovptot > (dtype) 0.01) {

              // Calculate riming term
              // Factor of liq water taken out because implicit
              zsnowrime = (dtype) 0.3*zcovptot*ptsphy*(*yrecldp)
                .rcl_const7s*zfallcorr*(MYPOW((zrho*zsnowcld*(*yrecldp).rcl_const1s),
                (*yrecldp).rcl_const8s));

              // Limit snow riming term
              zsnowrime = MYMIN(zsnowrime, (dtype) 1.0);

              zsolqb[3 + 5*(0)] = zsolqb[3 + 5*(0)] + zsnowrime;

            }

            //------------------------------------------------------------------
            // Riming of ice by cloud water - implicit in lwc
            // NOT YET ACTIVE
            //------------------------------------------------------------------
            //      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
            //
            //        ! Calculate riming term
            //        ! Factor of liq water taken out because implicit
            //        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
            //     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
            //
            //        ! Limit ice riming term
            //        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
            //
            //        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
            //
            //      ENDIF
          }

        }
        // on IWARMRAIN > 1


        //----------------------------------------------------------------------
        // 4.4a  MELTING OF SNOW and ICE
        //       with new implicit solver this also has to treat snow or ice
        //       precipitating from the level above... i.e. local ice AND flux.
        //       in situ ice and snow: could arise from LS advection or warming
        //       falling ice and snow: arrives by precipitation process
        //----------------------------------------------------------------------

        zicetot = zqxfg[1] + zqxfg[3];
        zmeltmax = (dtype) 0.0;

        // If there are frozen hydrometeors present and dry-bulb temperature > 0degC
        if (zicetot > zepsec && ztp1[jk_i] > rtt) {

          // Calculate subsaturation
          zsubsat = MYMAX(zqsice - zqx[4], (dtype) 0.0);

          // Calculate difference between dry-bulb (ZTP1) and the temperature
          // at which the wet-bulb=0degC (rtt-ZSUBSAT*....) using an approx.
          // Melting only occurs if the wet-bulb temperature >0
          // i.e. warming of ice particle due to melting > cooling
          // due to evaporation.
          ztdmtw0 = ztp1[jk_i] - rtt - zsubsat*(ztw1 + ztw2*(pap[jl + klon*(jk +
            klev*ibl)] - ztw3) - ztw4*(ztp1[jk_i] - ztw5));
          // Not implicit yet...
          // Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
          zcons1 =
            MYABS(ptsphy*((dtype) 1.0 + (dtype) 0.5*ztdmtw0) / (*yrecldp).rtaumel);
          zmeltmax = MYMAX(ztdmtw0*zcons1*zrldcp, (dtype) 0.0);
        }

        // Loop over frozen hydrometeors (ice, snow)
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          if (iphase[jm] == 2) {
            jn = imelt[jm];
            if (zmeltmax > zepsec && zicetot > zepsec) {
              // Apply melting in same proportion as frozen hydrometeor fractions
              zalfa = zqxfg[jm] / zicetot;
              zmelt = MYMIN(zqxfg[jm], zalfa*zmeltmax);
              // needed in first guess
              // This implies that zqpretot has to be recalculated below
              // since is not conserved here if ice falls and liquid doesn't
              zqxfg[jm] = zqxfg[jm] - zmelt;
              zqxfg[-1 + jn] = zqxfg[-1 + jn] + zmelt;
              zsolqa[-1 + jn + 5*jm] = zsolqa[-1 + jn + 5*jm] + zmelt;
              zsolqa[jm + 5*(-1 + jn)] = zsolqa[jm + 5*(-1 + jn)] - zmelt;
            }
          }
        }

        //----------------------------------------------------------------------
        // 4.4b  FREEZING of RAIN
        //----------------------------------------------------------------------

        // If rain present
        if (zqx[2] > zepsec) {

          if (ztp1[jk_i] <= rtt && ztp1[jk_im1] > rtt) {
            // Base of melting layer/top of refreezing layer so
            // store rain/snow fraction for precip type diagnosis
            // If mostly rain, then supercooled rain slow to freeze
            // otherwise faster to freeze (snow or ice pellets)
            zqpretot = MYMAX(zqx[3] + zqx[2], zepsec);
            prainfrac_toprfz[jl + klon*ibl] = zqx[2] / zqpretot;
            if (prainfrac_toprfz[jl + klon*ibl] > 0.8) {
              llrainliq = true;
            } else {
              llrainliq = false;
            }
          }

          // If temperature less than zero
          if (ztp1[jk_i] < rtt) {

            if (prainfrac_toprfz[jl + klon*ibl] > 0.8) {

              // Majority of raindrops completely melted
              // Refreezing is by slow heterogeneous freezing

              // Slope of rain particle size distribution
              zlambda =
                MYPOW(((*yrecldp).rcl_fac1 / (zrho*zqx[2])), (*yrecldp).rcl_fac2);

              // Calculate freezing rate based on Bigg(1953) and Wisner(1972)
              ztemp = (*yrecldp).rcl_fzrab*(ztp1[jk_i] - rtt);
              zfrz = ptsphy*((*yrecldp).rcl_const5r / zrho)*(MYEXP(ztemp) - (dtype) 1.)
                *(MYPOW(zlambda, (*yrecldp).rcl_const6r));
              zfrzmax = MYMAX(zfrz, (dtype) 0.0);

            } else {

              // Majority of raindrops only partially melted
              // Refreeze with a shorter timescale (reverse of melting...for now)

              zcons1 = MYABS(ptsphy*((dtype) 1.0 + (dtype) 0.5*(rtt - ztp1[jk_i])
                ) / (*yrecldp).rtaumel);
              zfrzmax = MYMAX((rtt - ztp1[jk_i])*zcons1*zrldcp, (dtype) 0.0);

            }

            if (zfrzmax > zepsec) {
              zfrz = MYMIN(zqx[2], zfrzmax);
              zsolqa[3 + 5*(2)] = zsolqa[3 + 5*(2)] + zfrz;
              zsolqa[2 + 5*(3)] = zsolqa[2 + 5*(3)] - zfrz;
            }
          }

        }


        //----------------------------------------------------------------------
        // 4.4c  FREEZING of LIQUID
        //----------------------------------------------------------------------
        // not implicit yet...
        zfrzmax = MYMAX(((*yrecldp).rthomo - ztp1[jk_i])*zrldcp, (dtype) 0.0);

        jm = 1;
        jn = imelt[-1 + jm];
        if (zfrzmax > zepsec && zqxfg[-1 + jm] > zepsec) {
          zfrz = MYMIN(zqxfg[-1 + jm], zfrzmax);
          zsolqa[-1 + jn + 5*(-1 + jm)] = zsolqa[-1 + jn + 5*(-1 + jm)] + zfrz;
          zsolqa[-1 + jm + 5*(-1 + jn)] = zsolqa[-1 + jm + 5*(-1 + jn)] - zfrz;
        }

        //----------------------------------------------------------------------
        // 4.5   EVAPORATION OF RAIN/SNOW
        //----------------------------------------------------------------------

        //----------------------------------------
        // Rain evaporation scheme from Sundquist
        //----------------------------------------
        if (ievaprain == 1) {

          // Rain


          zzrh = (*yrecldp).rprecrhmax + ((dtype) 1.0 - (*yrecldp).rprecrhmax)
            *zcovpmax / MYMAX(zepsec, (dtype) 1.0 - za[jk_i]);
          zzrh = MYMIN(MYMAX(zzrh, (*yrecldp).rprecrhmax), (dtype) 1.0);

          zqe = (zqx[4] - za[jk_i]*zqsliq) / MYMAX(zepsec, (dtype) 1.0 -
            za[jk_i]);
          //---------------------------------------------
          // humidity in moistest ZCOVPCLR part of domain
          //---------------------------------------------
          zqe = MYMAX((dtype) 0.0, MYMIN(zqe, zqsliq));
          llo1 = zcovpclr > zepsec && zqxfg[2] > zepsec && zqe < zzrh*zqsliq;

          if (llo1) {
            // note: zpreclr is a rain flux
            zpreclr = zqxfg[2]*zcovpclr / copysign(MYMAX(MYABS(zcovptot*zdtgdp),
              zepsilon), zcovptot*zdtgdp);

            //--------------------------------------
            // actual microphysics formula in zbeta
            //--------------------------------------

            zbeta1 = sqrt(pap[jl + klon*(jk + klev*ibl)] / paph[jl + klon*(klev + (klev
               + 1)*ibl)]) / (*yrecldp).rvrfactor*zpreclr / MYMAX(zcovpclr, zepsec);

            zbeta = rg*(*yrecldp).rpecons*(dtype) 0.5*(MYPOW(zbeta1, (dtype) 0.5777));

            zdenom = (dtype) 1.0 + zbeta*ptsphy*zcorqsliq;
            zdpr = zcovpclr*zbeta*(zqsliq - zqe) / zdenom*zdp*zrg_r;
            zdpevap = zdpr*zdtgdp;

            //---------------------------------------------------------
            // add evaporation term to explicit sink.
            // this has to be explicit since if treated in the implicit
            // term evaporation can not reduce rain to zero and model
            // produces small amounts of rainfall everywhere.
            //---------------------------------------------------------

            // Evaporate rain
            zevap = MYMIN(zdpevap, zqxfg[2]);

            zsolqa[4 + 5*(2)] = zsolqa[4 + 5*(2)] + zevap;
            zsolqa[2 + 5*(4)] = zsolqa[2 + 5*(4)] - zevap;

            //-------------------------------------------------------------
            // Reduce the total precip coverage proportional to evaporation
            // to mimic the previous scheme which had a diagnostic
            // 2-flux treatment, abandoned due to the new prognostic precip
            //-------------------------------------------------------------
            zcovptot = MYMAX((*yrecldp).rcovpmin, zcovptot - MYMAX((dtype) 0.0,
              (zcovptot - za[jk_i])*zevap / zqxfg[2]));

            // Update fg field
            zqxfg[2] = zqxfg[2] - zevap;

          }


          //---------------------------------------------------------
          // Rain evaporation scheme based on Abel and Boutle (2013)
          //---------------------------------------------------------
        } else if (ievaprain == 2) {


          //-----------------------------------------------------------------------
          // Calculate relative humidity limit for rain evaporation
          // to avoid cloud formation and saturation of the grid box
          //-----------------------------------------------------------------------
          // Limit RH for rain evaporation dependent on precipitation fraction
          zzrh = (*yrecldp).rprecrhmax + ((dtype) 1.0 - (*yrecldp).rprecrhmax)
            *zcovpmax / MYMAX(zepsec, (dtype) 1.0 - za[jk_i]);
          zzrh = MYMIN(MYMAX(zzrh, (*yrecldp).rprecrhmax), (dtype) 1.0);

          // Critical relative humidity
          //ZRHC=RAMID
          //ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
          // Increase RHcrit to 1.0 towards the surface (eta>0.8)
          //IF(ZSIGK > 0.8_JPRB) THEN
          //  ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
          //ENDIF
          //ZZRH = MIN(ZRHC,ZZRH)

          // Further limit RH for rain evaporation to 80% (RHcrit in free troposphere)
          zzrh = MYMIN((dtype) 0.8, zzrh);

          zqe = MYMAX((dtype) 0.0, MYMIN(zqx[4], zqsliq));

          llo1 = zcovpclr > zepsec && zqxfg[2] > zepsec && zqe < zzrh*zqsliq;

          if (llo1) {

            //-------------------------------------------
            // Abel and Boutle (2012) evaporation
            //-------------------------------------------
            // Calculate local precipitation (kg/kg)
            zpreclr = zqxfg[2] / zcovptot;

            // Fallspeed air density correction
            zfallcorr = MYPOW(((*yrecldp).rdensref / zrho), 0.4);

            // Saturation vapour pressure with respect to liquid phase
            zesatliq = rv / rd*((dtype)(r2es*MYEXP((r3les*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4les))));

            // Slope of particle size distribution
            zlambda = MYPOW(((*yrecldp).rcl_fac1 / (zrho*zpreclr)), (*yrecldp).rcl_fac2);                // ZPRECLR=kg/kg

            zevap_denom = (*yrecldp).rcl_cdenom1*zesatliq - (*yrecldp)
              .rcl_cdenom2*ztp1[jk_i]*zesatliq + (*yrecldp)
              .rcl_cdenom3*(MYPOW(ztp1[jk_i], (dtype) 3.))*pap[jl + klon*(jk +
              klev*ibl)];

            // Temperature dependent conductivity
            zcorr2 = (MYPOW((ztp1[jk_i] / (dtype) 273.), (dtype) 1.5))*(dtype)
              393. / (ztp1[jk_i] + (dtype) 120.);
            zka = (*yrecldp).rcl_ka273*zcorr2;

            zsubsat = MYMAX(zzrh*zqsliq - zqe, (dtype) 0.0);

            zbeta = ((dtype) 0.5 / zqsliq)*(MYPOW(ztp1[jk_i], (dtype) 2.))
              *zesatliq*(*yrecldp).rcl_const1r*(zcorr2 / zevap_denom)*((dtype) 0.78 /
              (MYPOW(zlambda, (*yrecldp).rcl_const4r)) + (*yrecldp)
              .rcl_const2r*(MYPOW((zrho*zfallcorr), (dtype) 0.5)) / ((MYPOW(zcorr2,
              (dtype) 0.5))*(MYPOW(zlambda, (*yrecldp).rcl_const3r))));

            zdenom = (dtype) 1.0 + zbeta*ptsphy;                //*ZCORQSLIQ(JL)
            zdpevap = zcovpclr*zbeta*ptsphy*zsubsat / zdenom;

            //---------------------------------------------------------
            // Add evaporation term to explicit sink.
            // this has to be explicit since if treated in the implicit
            // term evaporation can not reduce rain to zero and model
            // produces small amounts of rainfall everywhere.
            //---------------------------------------------------------

            // Limit rain evaporation
            zevap = MYMIN(zdpevap, zqxfg[2]);

            zsolqa[4 + 5*(2)] = zsolqa[4 + 5*(2)] + zevap;
            zsolqa[2 + 5*(4)] = zsolqa[2 + 5*(4)] - zevap;

            //-------------------------------------------------------------
            // Reduce the total precip coverage proportional to evaporation
            // to mimic the previous scheme which had a diagnostic
            // 2-flux treatment, abandoned due to the new prognostic precip
            //-------------------------------------------------------------
            zcovptot = MYMAX((*yrecldp).rcovpmin, zcovptot - MYMAX((dtype) 0.0,
              (zcovptot - za[jk_i])*zevap / zqxfg[2]));

            // Update fg field
            zqxfg[2] = zqxfg[2] - zevap;

          }

        }
        // on IEVAPRAIN

        //----------------------------------------------------------------------
        // 4.5   EVAPORATION OF SNOW
        //----------------------------------------------------------------------
        // Snow
        if (ievapsnow == 1) {

          zzrh = (*yrecldp).rprecrhmax + ((dtype) 1.0 - (*yrecldp).rprecrhmax)
            *zcovpmax / MYMAX(zepsec, (dtype) 1.0 - za[jk_i]);
          zzrh = MYMIN(MYMAX(zzrh, (*yrecldp).rprecrhmax), (dtype) 1.0);
          zqe = (zqx[4] - za[jk_i]*zqsice) / MYMAX(zepsec, (dtype) 1.0 -
            za[jk_i]);

          //---------------------------------------------
          // humidity in moistest ZCOVPCLR part of domain
          //---------------------------------------------
          zqe = MYMAX((dtype) 0.0, MYMIN(zqe, zqsice));
          llo1 = zcovpclr > zepsec && zqxfg[3] > zepsec && zqe < zzrh*zqsice;

          if (llo1) {
            // note: zpreclr is a rain flux a
            zpreclr = zqxfg[3]*zcovpclr / copysign(MYMAX(MYABS(zcovptot*zdtgdp),
              zepsilon), zcovptot*zdtgdp);

            //--------------------------------------
            // actual microphysics formula in zbeta
            //--------------------------------------

            zbeta1 = sqrt(pap[jl + klon*(jk + klev*ibl)] / paph[jl + klon*(klev + (klev
               + 1)*ibl)]) / (*yrecldp).rvrfactor*zpreclr / MYMAX(zcovpclr, zepsec);

            zbeta = rg*(*yrecldp).rpecons*(MYPOW(zbeta1, (dtype) 0.5777));

            zdenom = (dtype) 1.0 + zbeta*ptsphy*zcorqsice;
            zdpr = zcovpclr*zbeta*(zqsice - zqe) / zdenom*zdp*zrg_r;
            zdpevap = zdpr*zdtgdp;

            //---------------------------------------------------------
            // add evaporation term to explicit sink.
            // this has to be explicit since if treated in the implicit
            // term evaporation can not reduce snow to zero and model
            // produces small amounts of snowfall everywhere.
            //---------------------------------------------------------

            // Evaporate snow
            zevap = MYMIN(zdpevap, zqxfg[3]);

            zsolqa[4 + 5*(3)] = zsolqa[4 + 5*(3)] + zevap;
            zsolqa[3 + 5*(4)] = zsolqa[3 + 5*(4)] - zevap;

            //-------------------------------------------------------------
            // Reduce the total precip coverage proportional to evaporation
            // to mimic the previous scheme which had a diagnostic
            // 2-flux treatment, abandoned due to the new prognostic precip
            //-------------------------------------------------------------
            zcovptot = MYMAX((*yrecldp).rcovpmin, zcovptot - MYMAX((dtype) 0.0,
              (zcovptot - za[jk_i])*zevap / zqxfg[3]));

            //Update first guess field
            zqxfg[3] = zqxfg[3] - zevap;

          }
          //---------------------------------------------------------
        } else if (ievapsnow == 2) {



          //-----------------------------------------------------------------------
          // Calculate relative humidity limit for snow evaporation
          //-----------------------------------------------------------------------
          zzrh = (*yrecldp).rprecrhmax + ((dtype) 1.0 - (*yrecldp).rprecrhmax)
            *zcovpmax / MYMAX(zepsec, (dtype) 1.0 - za[jk_i]);
          zzrh = MYMIN(MYMAX(zzrh, (*yrecldp).rprecrhmax), (dtype) 1.0);
          zqe = (zqx[4] - za[jk_i]*zqsice) / MYMAX(zepsec, (dtype) 1.0 -
            za[jk_i]);

          //---------------------------------------------
          // humidity in moistest ZCOVPCLR part of domain
          //---------------------------------------------
          zqe = MYMAX((dtype) 0.0, MYMIN(zqe, zqsice));
          llo1 = zcovpclr > zepsec && zqx[3] > zepsec && zqe < zzrh*zqsice;

          if (llo1) {

            // Calculate local precipitation (kg/kg)
            zpreclr = zqx[3] / zcovptot;
            zvpice = ((dtype)(r2es*MYEXP((r3ies*(ztp1[jk_i] - rtt))/(ztp1[jk_i] - r4ies))))*rv / rd;

            // Particle size distribution
            // ZTCG increases Ni with colder temperatures - essentially a
            // Fletcher or Meyers scheme?
            ztcg = (dtype) 1.0;                //v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
            // ZFACX1I modification is based on Andrew Barrett's results
            zfacx1s = (dtype) 1.0;                //v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB

            zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk_i] + 
              pap[jl + klon*(jk + klev*ibl)]*(*yrecldp).rcl_apb3*(MYPOW(ztp1[jk_i], 3));
            zcorrfac = MYPOW((1.0 / zrho), 0.5);
            zcorrfac2 =
              (MYPOW((ztp1[jk_i] / 273.0), 1.5))*(393.0 / (ztp1[jk_i] + 120.0))
              ;

            zpr02 = zrho*zpreclr*(*yrecldp).rcl_const1s / (ztcg*zfacx1s);

            zterm1 = (zqsice - zqe)*(MYPOW(ztp1[jk_i], 2))
              *zvpice*zcorrfac2*ztcg*(*yrecldp).rcl_const2s*zfacx1s /
              (zrho*zaplusb*zqsice);
            zterm2 = 0.65*(*yrecldp).rcl_const6s*(MYPOW(zpr02, (*yrecldp).rcl_const4s)) +
               (*yrecldp).rcl_const3s*(MYPOW(zcorrfac, 0.5))*(MYPOW(zrho, 0.5))*(MYPOW(zpr02,
               (*yrecldp).rcl_const5s)) / (MYPOW(zcorrfac2, 0.5));

            zdpevap = MYMAX(zcovpclr*zterm1*zterm2*ptsphy, (dtype) 0.0);

            //--------------------------------------------------------------------
            // Limit evaporation to snow amount
            //--------------------------------------------------------------------
            zevap = MYMIN(zdpevap, zevaplimice);
            zevap = MYMIN(zevap, zqx[3]);


            zsolqa[4 + 5*(3)] = zsolqa[4 + 5*(3)] + zevap;
            zsolqa[3 + 5*(4)] = zsolqa[3 + 5*(4)] - zevap;

            //-------------------------------------------------------------
            // Reduce the total precip coverage proportional to evaporation
            // to mimic the previous scheme which had a diagnostic
            // 2-flux treatment, abandoned due to the new prognostic precip
            //-------------------------------------------------------------
            zcovptot = MYMAX((*yrecldp).rcovpmin, zcovptot - MYMAX((dtype) 0.0,
              (zcovptot - za[jk_i])*zevap / zqx[3]));

            //Update first guess field
            zqxfg[3] = zqxfg[3] - zevap;

          }

        }
        // on IEVAPSNOW

        //--------------------------------------
        // Evaporate small precipitation amounts
        //--------------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          if (llfall[jm]) {
            if (zqxfg[jm] < (*yrecldp).rlmin) {
              zsolqa[4 + 5*jm] = zsolqa[4 + 5*jm] + zqxfg[jm];
              zsolqa[jm + 5*(4)] = zsolqa[jm + 5*(4)] - zqxfg[jm];
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
        zanew = (za[jk_i] + zsolac) / ((dtype) 1.0 + zsolab);
        zanew = MYMIN(zanew, (dtype) 1.0);
        if (zanew < (*yrecldp).ramin) {
          zanew = (dtype) 0.0;
        }
        zda = zanew - zaorig;
        //---------------------------------
        // variables needed for next level
        //---------------------------------
        zanewm1 = zanew;

        //--------------------------------
        // 5.2 solver for the microphysics
        //--------------------------------

        //--------------------------------------------------------------
        // Truncate explicit sinks to avoid negatives
        // Note: Species are treated in the order in which they run out
        // since the clipping will alter the balance for the other vars
        //--------------------------------------------------------------

        for (jm = 0; jm <= 5 + -1; jm += 1) {
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            llindex3[jn + 5*jm] = false;
          }
          zsinksum[jm] = (dtype) 0.0;
        }

        //----------------------------
        // collect sink terms and mark
        //----------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            zsinksum[jm] = zsinksum[jm] - zsolqa[jm + 5*jn];                // +ve total is bad
          }
        }

        //---------------------------------------
        // calculate overshoot and scaling factor
        //---------------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zmax = MYMAX(zqx[jm], zepsec);
          zrat = MYMAX(zsinksum[jm], zmax);
          zratio[jm] = zmax / zrat;
        }

        //--------------------------------------------
        // scale the sink terms, in the correct order,
        // recalculating the scale factor each time
        //--------------------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zsinksum[jm] = (dtype) 0.0;
        }

        //----------------
        // recalculate sum
        //----------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          psum_solqa = 0.0;
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            psum_solqa = psum_solqa + zsolqa[jm + 5*jn];
          }
          // ZSINKSUM(JL,JM)=ZSINKSUM(JL,JM)-SUM(ZSOLQA(JL,JM,1:NCLV))
          zsinksum[jm] = zsinksum[jm] - psum_solqa;
          //---------------------------
          // recalculate scaling factor
          //---------------------------
          zmm = MYMAX(zqx[jm], zepsec);
          zrr = MYMAX(zsinksum[jm], zmm);
          zratio[jm] = zmm / zrr;
          //------
          // scale
          //------
          zzratio = zratio[jm];
          //DIR$ IVDEP
          //DIR$ PREFERVECTOR
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            if (zsolqa[jm + 5*jn] < (dtype) 0.0) {
              zsolqa[jm + 5*jn] = zsolqa[jm + 5*jn]*zzratio;
              zsolqa[jn + 5*jm] = zsolqa[jn + 5*jm]*zzratio;
            }
          }
        }

        //--------------------------------------------------------------
        // 5.2.2 Solver
        //------------------------

        //------------------------
        // set the LHS of equation
        //------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            //----------------------------------------------
            // diagonals: microphysical sink terms+transport
            //----------------------------------------------
            if (jn + 1 == jm + 1) {
              zqlhs[jn + 5*jm] = (dtype) 1.0 + zfallsink[jm];
              for (jo = 0; jo <= 5 + -1; jo += 1) {
                zqlhs[jn + 5*jm] = zqlhs[jn + 5*jm] + zsolqb[jo + 5*jn];
              }
              //------------------------------------------
              // non-diagonals: microphysical source terms
              //------------------------------------------
            } else {
              zqlhs[jn + 5*jm] = -zsolqb[jn + 5*jm];                  // here is the delta T - missing from doc.
            }
          }
        }

        //------------------------
        // set the RHS of equation
        //------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          //---------------------------------
          // sum the explicit source and sink
          //---------------------------------
          zexplicit = (dtype) 0.0;
          for (jn = 0; jn <= 5 + -1; jn += 1) {
            zexplicit = zexplicit + zsolqa[jm + 5*jn];                // sum over middle index
          }
          zqxn[jm] = zqx[jm] + zexplicit;
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
        for (jn = 0; jn <= 5 - 1 + -1; jn += 1) {
          // number of steps
          for (jm = jn + 1; jm <= 5 + -1; jm += 1) {
            // row index
            zqlhs[jm + 5*jn] = zqlhs[jm + 5*jn] / zqlhs[jn + 5*jn];
            for (ik = jn + 1; ik <= 5 + -1; ik += 1) {
              // column index
              zqlhs[jm + 5*ik] = zqlhs[jm + 5*ik] - zqlhs[jm + 5*jn]*zqlhs[jn + 5*ik];
            }
          }
        }

        // Backsubstitution
        //  step 1
        for (jn = 1; jn <= 5 + -1; jn += 1) {
          for (jm = 0; jm <= jn + 1 - 1 + -1; jm += 1) {
            zqxn[jn] = zqxn[jn] - zqlhs[jn + 5*jm]*zqxn[jm];
          }
        }
        //  step 2
        zqxn[4] = zqxn[4] / zqlhs[4 + 5*(4)];
        for (jn = -2 + 5; jn >= 1 + -1; jn += -1) {
          for (jm = jn + 1; jm <= 5 + -1; jm += 1) {
            zqxn[jn] = zqxn[jn] - zqlhs[jn + 5*jm]*zqxn[jm];
          }
          zqxn[jn] = zqxn[jn] / zqlhs[jn + 5*jn];
        }

        // Ensure no small values (including negatives) remain in cloud variables nor
        // precipitation rates.
        // Evaporate l,i,r,s to water vapour. Latent heating taken into account below
        for (jn = 0; jn <= 5 - 1 + -1; jn += 1) {
          if (zqxn[jn] < zepsec) {
            zqxn[4] = zqxn[4] + zqxn[jn];
            zqxn[jn] = (dtype) 0.0;
          }
        }

        //--------------------------------
        // variables needed for next level
        //--------------------------------
        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zqxnm1[jm] = zqxn[jm];
          zqxn2d[jm] = zqxn[jm];
        }

        //------------------------------------------------------------------------
        // 5.3 Precipitation/sedimentation fluxes to next level
        //     diagnostic precipitation fluxes
        //     It is this scaled flux that must be used for source to next layer
        //------------------------------------------------------------------------

        for (jm = 0; jm <= 5 + -1; jm += 1) {
          zpfplsx[jk_ip1 + 2*jm] = zfallsink[jm]*zqxn[jm]*zrdtgdp;
        }

        // Ensure precipitation fraction is zero if no precipitation
        zqpretot =
          zpfplsx[jk_ip1 + 2*(3)] + zpfplsx[jk_ip1 + 2*(2)];
        if (zqpretot < zepsec) {
          zcovptot = (dtype) 0.0;
        }

        //######################################################################
        //              6  *** UPDATE TENDANCIES ***
        //######################################################################

        //--------------------------------
        // 6.1 Temperature and CLV budgets
        //--------------------------------

        for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {

          // calculate fluxes in and out of box for conservation of TL
          zfluxq[jm] = zpsupsatsrce[jm] + zconvsrce[jm] + zfallsrce[jm] -
            (zfallsink[jm] + zconvsink[jm])*zqxn[jm];

          if (iphase[jm] == 1) {
            tendency_loc_t[jl + klon*(jk + klev*ibl)] = tendency_loc_t[jl
               + klon*(jk + klev*ibl)] + ralvdcp*(zqxn[jm] - zqx[jm] -
              zfluxq[jm])*zqtmst;
          }

          if (iphase[jm] == 2) {
            tendency_loc_t[jl + klon*(jk + klev*ibl)] = tendency_loc_t[jl
               + klon*(jk + klev*ibl)] + ralsdcp*(zqxn[jm] - zqx[jm] -
              zfluxq[jm])*zqtmst;
          }

          //----------------------------------------------------------------------
          // New prognostic tendencies - ice,liquid rain,snow
          // Note: CLV arrays use PCLV in calculation of tendency while humidity
          //       uses ZQX. This is due to clipping at start of cloudsc which
          //       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
          //----------------------------------------------------------------------
          tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*ibl))] = tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*ibl))]
              + (zqxn[jm] - zqx0[jm])*zqtmst;

        }

        //----------------------
        // 6.2 Humidity budget
        //----------------------
        tendency_loc_q[jl + klon*(jk + klev*ibl)] = tendency_loc_q[jl +
          klon*(jk + klev*ibl)] + (zqxn[4] - zqx[4])*zqtmst;

        //-------------------
        // 6.3 cloud cover
        //-----------------------
        tendency_loc_a[jl + klon*(jk + klev*ibl)] =
          tendency_loc_a[jl + klon*(jk + klev*ibl)] + zda*zqtmst;

        //--------------------------------------------------
        // Copy precipitation fraction into output variable
        //-------------------------------------------------
        pcovptot[jl + klon*(jk + klev*ibl)] = zcovptot;

      }
    }

    // on vertical level JK
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
    pfplsl[jl + klon*(jk + (klev + 1)*ibl)] =
      zpfplsx[jk_i + 2*(2)] + zpfplsx[jk_i + 2*(0)];
    pfplsn[jl + klon*(jk + (klev + 1)*ibl)] =
      zpfplsx[jk_i + 2*(3)] + zpfplsx[jk_i + 2*(1)];

    if (1 <= jk + 1 && jk + 1 <= klev) {

      zgdph_r = -zrg_r*(paph[jl + klon*(1 + jk + (klev + 1)*ibl)] - paph[jl + klon*(jk
        + (klev + 1)*ibl)])*zqtmst;
      pfsqlf[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqlf[jl + klon*(jk + (klev + 1)*ibl)];
      pfsqif[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqif[jl + klon*(jk + (klev + 1)*ibl)];
      pfsqrf[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqlf[jl + klon*(jk + (klev + 1)*ibl)];
      pfsqsf[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqif[jl + klon*(jk + (klev + 1)*ibl)];
      pfcqlng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqlng[jl + klon*(jk + (klev + 1)*ibl)];
      pfcqnng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqnng[jl + klon*(jk + (klev + 1)*ibl)];
      pfcqrng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqlng[jl + klon*(jk + (klev + 1)*ibl)];
      pfcqsng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqnng[jl + klon*(jk + (klev + 1)*ibl)];
      pfsqltur[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqltur[jl + klon*(jk + (klev + 1)*ibl)];
      pfsqitur[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfsqitur[jl + klon*(jk + (klev + 1)*ibl)];

      zalfaw = zfoealfa;

      // Liquid , LS scheme minus detrainment
      pfsqlf[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqlf[jl + klon*(1 + jk + (klev +
        1)*ibl)] + (zqxn2d[0] - zqx0[0] + pvfl[jl + klon*(jk + klev*ibl)
        ]*ptsphy - zalfaw*plude[jl + klon*(jk + klev*ibl)])*zgdph_r;
      // liquid, negative numbers
      pfcqlng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqlng[jl + klon*(1 + jk + (klev + 1)*ibl)] + zlneg[0]*zgdph_r;

      // liquid, vertical diffusion
      pfsqltur[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqltur[jl + klon*(1 + jk +
        (klev + 1)*ibl)] + pvfl[jl + klon*(jk + klev*ibl)]*ptsphy*zgdph_r;

      // Rain, LS scheme
      pfsqrf[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqrf[jl + klon*(1 + jk + (klev +
        1)*ibl)] + (zqxn2d[2] - zqx0[2])*zgdph_r;
      // rain, negative numbers
      pfcqrng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqrng[jl + klon*(1 + jk + (klev + 1)*ibl)] + zlneg[2]*zgdph_r;

      // Ice , LS scheme minus detrainment
      pfsqif[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqif[jl + klon*(1 + jk + (klev +
        1)*ibl)] + (zqxn2d[1] - zqx0[1] + pvfi[jl + klon*(jk + klev*ibl)
        ]*ptsphy - ((dtype) 1.0 - zalfaw)*plude[jl + klon*(jk + klev*ibl)])*zgdph_r;
      // ice, negative numbers
      pfcqnng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqnng[jl + klon*(1 + jk + (klev + 1)*ibl)] + zlneg[1]*zgdph_r;

      // ice, vertical diffusion
      pfsqitur[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqitur[jl + klon*(1 + jk +
        (klev + 1)*ibl)] + pvfi[jl + klon*(jk + klev*ibl)]*ptsphy*zgdph_r;

      // snow, LS scheme
      pfsqsf[jl + klon*(1 + jk + (klev + 1)*ibl)] = pfsqsf[jl + klon*(1 + jk + (klev +
        1)*ibl)] + (zqxn2d[3] - zqx0[3])*zgdph_r;
      // snow, negative numbers
      pfcqsng[jl + klon*(1 + jk + (klev + 1)*ibl)] =
        pfcqsng[jl + klon*(1 + jk + (klev + 1)*ibl)] + zlneg[3]*zgdph_r;
    }

    //-----------------------------------
    // enthalpy flux due to precipitation
    //-----------------------------------
    pfhpsl[jl + klon*(jk + (klev + 1)*ibl)] =
      -rlvtt*pfplsl[jl + klon*(jk + (klev + 1)*ibl)];
    pfhpsn[jl + klon*(jk + (klev + 1)*ibl)] =
      -rlstt*pfplsn[jl + klon*(jk + (klev + 1)*ibl)];
  }
}

