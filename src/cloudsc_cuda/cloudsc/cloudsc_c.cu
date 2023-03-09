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

__global__ void cloudsc_c(int kidia, int kfdia, int klon, double ptsphy,
  const double * __restrict__  pt,
  const double * __restrict__  pq, const double * __restrict__  tendency_tmp_t,
  const double * __restrict__  tendency_tmp_q, const double * __restrict__  tendency_tmp_a,
  const double * __restrict__  tendency_tmp_cld, double * __restrict__  tendency_loc_t,
  double * __restrict__  tendency_loc_q, double * __restrict__  tendency_loc_a,
  double * __restrict__  tendency_loc_cld, const double * __restrict__  pvfa,
  const double * __restrict__  pvfl, const double * __restrict__  pvfi, const double * __restrict__  pdyna,
  const double * __restrict__  pdynl, const double * __restrict__  pdyni, const double * __restrict__  phrsw,
  double * __restrict__  phrlw, const double * __restrict__  pvervel, const double * __restrict__  pap,
  const double * __restrict__  paph, const double * __restrict__  plsm,
  const int *  ktype, const double * __restrict__  plu, double * __restrict__  plude,
  const double * __restrict__  psnde, const double * __restrict__  pmfu, const double * __restrict__  pmfd,
  const double * __restrict__  pa, const double * __restrict__  pclv, const double * __restrict__  psupsat,
  const double * __restrict__  plcrit_aer, const double * __restrict__  picrit_aer,
  const double * __restrict__  pre_ice, const double * __restrict__  pccn, const double * __restrict__  pnice,
  double * __restrict__  pcovptot, double * __restrict__  prainfrac_toprfz,
  double * __restrict__  pfsqlf, double * __restrict__  pfsqif, double * __restrict__  pfcqnng,
  double * __restrict__  pfcqlng, double * __restrict__  pfsqrf, double * __restrict__  pfsqsf,
  double * __restrict__  pfcqrng, double * __restrict__  pfcqsng,
  double * __restrict__  pfsqltur, double * __restrict__  pfsqitur,
  double * __restrict__  pfplsl, double * __restrict__  pfplsn, double * __restrict__  pfhpsl,
  double * __restrict__  pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  double rg, double rd, double rcpd, double retv, double rlvtt, double rlstt, double rlmlt, double rtt,
  double rv, double r2es, double r3les, double r3ies, double r4les, double r4ies, double r5les,
  double r5ies, double r5alvcp, double r5alscp, double ralvdcp, double ralsdcp, double ralfdcp,
  double rtwat, double rtice, double rticecu, double rtwat_rtice_r, double rtwat_rticecu_r,
  double rkoop1, double rkoop2) {

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

  double zlcond1, zlcond2, zlevapl, zlevapi, zrainaut, zsnowaut, zliqcld, zicecld;
  double zlevap, zleros;
  //  condensation and evaporation terms
  // autoconversion terms
  double zfokoop;
  double zfoealfa[klev + 1];
  double zicenuclei;  // number concentration of ice nuclei

  double zlicld;
  double zacond;
  double zaeros;
  double zlfinalsum;
  double zdqs;
  double ztold;
  double zqold;
  double zdtgdp;
  double zrdtgdp;
  double ztrpaus;
  double zcovpclr;
  double zpreclr;
  double zcovptot;
  double zcovpmax;
  double zqpretot;
  double zdpevap;
  double zdtforc;
  double zdtdiab;
  double ztp1[klev];
  double zldefr;
  double zldifdt;
  double zdtgdpf;
  double zlcust[5];
  double zacust;
  double zmf;

  double zrho;
  double ztmp1, ztmp2, ztmp3;
  double ztmp4, ztmp5, ztmp6, ztmp7;
  double zalfawm;

  // Accumulators of A,B,and C factors for cloud equations
  double zsolab;  // -ve implicit CC
  double zsolac;  // linear CC
  double zanew;
  double zanewm1;

  double zgdp;

  //---for flux calculation
  double zda;
  double zli[klev], za[klev];
  double zaorig[klev];  // start of scheme value for CC

  int llflag;
  int llo1;

  int icall, ik, jk, jl, jm, jn, jo, jlen, is;

  double zdp, zpaphd;

  double zalfa;
  // & ZALFACU, ZALFALS
  double zalfaw;
  double zbeta, zbeta1;
  //REAL(KIND=JPRB) :: ZBOTT
  double zcfpr;
  double zcor;
  double zcdmax;
  double zmin;
  double zlcondlim;
  double zdenom;
  double zdpmxdt;
  double zdpr;
  double zdtdp;
  double ze;
  double zepsec;
  double zfac, zfaci, zfacw;
  double zgdcp;
  double zinew;
  double zlcrit;
  double zmfdn;
  double zprecip;
  double zqe;
  double zqsat, zqtmst, zrdcp;
  double zrhc, zsig, zsigk;
  double zwtot;
  double zzco, zzdl, zzrh, zzzdt, zqadj;
  double zqnew, ztnew;
  double zrg_r, zgdph_r, zcons1, zcond, zcons1a;
  double zlfinal;
  double zmelt;
  double zevap;
  double zfrz;
  double zvpliq, zvpice;
  double zadd, zbdd, zcvds, zice0, zdepos;
  double zsupsat;
  double zfall;
  double zre_ice;
  double zrldcp;
  double zqp1env;

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
  int llindex3[5 * 5];  // index variable
  double zmax;
  double zrat;
  int iorder[5];  // array for sorting explicit terms

  double zliqfrac[klev];  // cloud liquid water fraction: ql/(ql+qi)
  double zicefrac[klev];  // cloud ice water fraction: qi/(ql+qi)
  double zqxn[5];  // new values for zqx at time+1
  double zqxfg[5];  // first guess values including precip
  double zqxnm1[5];  // new values for zqx at time+1 at level above
  double zfluxq[5];  // fluxes convergence of species (needed?)
  // Keep the following for possible future total water variance scheme?
  //REAL(KIND=JPRB) :: ZTL(KLON,KLEV)       ! liquid water temperature
  //REAL(KIND=JPRB) :: ZABETA(KLON,KLEV)    ! cloud fraction
  //REAL(KIND=JPRB) :: ZVAR(KLON,KLEV)      ! temporary variance
  //REAL(KIND=JPRB) :: ZQTMIN(KLON,KLEV)
  //REAL(KIND=JPRB) :: ZQTMAX(KLON,KLEV)

  double zmeltmax;
  double zfrzmax;
  double zicetot;


  double zqsmix[klev];  // diagnostic mixed phase saturation
  //REAL(KIND=JPRB) :: ZQSBIN(KLON,KLEV) ! binary switched ice/liq saturation
  double zqsliq[klev];  // liquid water saturation
  double zqsice[klev];  // ice water saturation

  //REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
  //REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
  //REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice

  double zfoeewmt[klev];
  double zfoeew[klev];
  double zfoeeliqt[klev];
  //REAL(KIND=JPRB) :: ZFOEEICET(KLON,KLEV)

  double zdqsliqdt, zdqsicedt, zdqsmixdt;
  double zcorqsliq;
  double zcorqsice;
  //REAL(KIND=JPRB) :: ZCORQSBIN(KLON)
  double zcorqsmix;
  double zevaplimliq, zevaplimice, zevaplimmix;

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

  double zsolqa[5 * 5];  // explicit sources and sinks
  double zsolqb[5 * 5];  // implicit sources and sinks
  // e.g. microphysical pathways between ice variables.
  double zqlhs[5 * 5];  // n x n matrix storing the LHS of implicit solver
  double zvqx[5];  // fall speeds of three categories
  double zexplicit;
  double zratio[5], zsinksum[5];

  // for sedimentation source/sink terms
  double zfallsink[5];
  double zfallsrce[5];

  // for convection detrainment source and subsidence source/sink terms
  double zconvsrce[5];
  double zconvsink[5];

  // for supersaturation source term from previous timestep
  double zpsupsatsrce[5];

  // Numerical fit to wet bulb temperature
  double ztw1 = (double) 1329.31;
  double ztw2 = (double) 0.0074615;
  double ztw3 = (double) 0.85E5;
  double ztw4 = (double) 40.637;
  double ztw5 = (double) 275.0;

  double zsubsat;  // Subsaturation for snow melting term
  double ztdmtw0;  // Diff between dry-bulb temperature and
  // temperature when wet-bulb = 0degC

  // Variables for deposition term
  double ztcg;  // Temperature dependent function for ice PSD
  double zfacx1i, zfacx1s;  // PSD correction factor
  double zaplusb, zcorrfac, zcorrfac2, zpr02, zterm1, zterm2;  // for ice dep
  double zcldtopdist;  // Distance from cloud top
  double zinfactor;  // No. of ice nuclei factor for deposition

  // Autoconversion/accretion/riming/evaporation
  int iwarmrain;
  int ievaprain;
  int ievapsnow;
  int idepice;
  double zrainacc;
  double zraincld;
  double zsnowrime;
  double zsnowcld;
  double zesatliq;
  double zfallcorr;
  double zlambda;
  double zevap_denom;
  double zcorr2;
  double zka;
  double zconst;
  double ztemp;

  // Rain freezing
  int llrainliq;  // True if majority of raindrops are liquid (no ice core)

  //----------------------------
  // End: new microphysics
  //----------------------------

  //----------------------
  // SCM budget statistics
  //----------------------
  double zrain;

  double zhook_handle;
  double ztmpl, ztmpi, ztmpa;

  double zmm, zrr;
  double zrg;

  double zzsum, zzratio;
  double zepsilon;

  double zcond1, zqp;

  double psum_solqa;
  int ibl;
  int i_llfall_0;
  double zqx[5 * klev];
  double zqx0[5 * klev];
  double zpfplsx[5 * (klev + 1)];
  double zlneg[5 * klev];
  double zqxn2d[5 * klev];

  jl = threadIdx.x;
  ibl = blockIdx.z; 


  //===============================================================================
  //IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)

  //===============================================================================
  //  0.0     Beginning of timestep book-keeping
  //----------------------------------------------------------------------


  //######################################################################
  //             0.  *** SET UP CONSTANTS ***
  //######################################################################

  zepsilon = (double) 100.*DBL_EPSILON;

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
  zqtmst = (double) 1.0 / ptsphy;
  zgdcp = rg / rcpd;
  zrdcp = rd / rcpd;
  zcons1a = rcpd / (rlmlt*rg*(*yrecldp).rtaumel);
  zepsec = (double) 1.E-14;
  zrg_r = (double) 1.0 / rg;
  zrldcp = (double) 1.0 / (ralsdcp - ralvdcp);

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
    tendency_loc_t[jl + klon*(jk + klev*(ibl))] = (double) 0.0;
    tendency_loc_q[jl + klon*(jk + klev*(ibl))] = (double) 0.0;
    tendency_loc_a[jl + klon*(jk + klev*(ibl))] = (double) 0.0;
  }
  for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
    for (jk = 0; jk <= klev + -1; jk += 1) {
      tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*(ibl)))] = (double) 0.0;
    }
  }

  //-- These were uninitialized : meaningful only when we compare error differences
  for (jk = 0; jk <= klev + -1; jk += 1) {
    pcovptot[jl + klon*(jk + klev*(ibl))] = (double) 0.0;
    tendency_loc_cld[jl + klon*(jk + klev*(4 + 5*(ibl)))] = (double) 0.0
      ;
  }

  // -------------------------
  // set up fall speeds in m/s
  // -------------------------
  zvqx[4] = (double) 0.0;
  zvqx[0] = (double) 0.0;
  zvqx[1] = (*yrecldp).rvice;
  zvqx[2] = (*yrecldp).rvrain;
  zvqx[3] = (*yrecldp).rvsnow;
  for (i_llfall_0 = 0; i_llfall_0 <= 5 + -1; i_llfall_0 += 1) {
    llfall[i_llfall_0] = false;
  }
  for (jm = 0; jm <= 5 + -1; jm += 1) {
    if (zvqx[jm] > (double) 0.0) {
      llfall[jm] = true;
    }
    // falling species
  }
  // Set LLFALL to false for ice (but ice still sediments!)
  // Need to rationalise this at some point
  llfall[1] = false;


  //######################################################################
  //             1.  *** INITIAL VALUES FOR VARIABLES ***
  //######################################################################


  // ----------------------
  // non CLV initialization
  // ----------------------
  for (jk = 0; jk <= klev + -1; jk += 1) {
    ztp1[jk] = pt[jl + klon*(jk + klev*(ibl))] + ptsphy*tendency_tmp_t[
      jl + klon*(jk + klev*(ibl))];
    zqx[jk + klev*(4)] = pq[jl + klon*(jk + klev*(ibl))] +
      ptsphy*tendency_tmp_q[jl + klon*(jk + klev*(ibl))];
    zqx0[jk + klev*(4)] = pq[jl + klon*(jk + klev*(ibl))] +
      ptsphy*tendency_tmp_q[jl + klon*(jk + klev*(ibl))];
    za[jk] = pa[jl + klon*(jk + klev*(ibl))] + ptsphy*tendency_tmp_a[jl
      + klon*(jk + klev*(ibl))];
    zaorig[jk] = pa[jl + klon*(jk + klev*(ibl))] + ptsphy*tendency_tmp_a[
       jl + klon*(jk + klev*(ibl))];
  }

  // -------------------------------------
  // initialization for CLV family
  // -------------------------------------
  for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
    for (jk = 0; jk <= klev + -1; jk += 1) {
      zqx[jk + klev*jm] = pclv[jl + klon*(jk + klev*(jm + 5*(ibl)))] +
        ptsphy*tendency_tmp_cld[jl + klon*(jk + klev*(jm + 5*(ibl)))];
      zqx0[jk + klev*jm] = pclv[jl + klon*(jk + klev*(jm + 5*(ibl)))] +
        ptsphy*tendency_tmp_cld[jl + klon*(jk + klev*(jm + 5*(ibl)))];
    }
  }

  //-------------
  // zero arrays
  //-------------
  for (jm = 0; jm <= 5 + -1; jm += 1) {
    for (jk = 0; jk <= klev + 1 + -1; jk += 1) {
      zpfplsx[jk + (klev + 1)*jm] = (double) 0.0;          // precip fluxes
    }
  }

  for (jm = 0; jm <= 5 + -1; jm += 1) {
    for (jk = 0; jk <= klev + -1; jk += 1) {
      zqxn2d[jk + klev*jm] = (double) 0.0;          // end of timestep values in 2D
      zlneg[jk + klev*jm] = (double) 0.0;          // negative input check
    }
  }

  prainfrac_toprfz[jl + klon*(ibl)] = (double) 0.0;      // rain fraction at top of refreezing layer
  llrainliq = true;      // Assume all raindrops are liquid initially

  // ----------------------------------------------------
  // Tidy up very small cloud cover or total cloud water
  // ----------------------------------------------------
  for (jk = 0; jk <= klev + -1; jk += 1) {
    if (zqx[jk + klev*(0)] + zqx[jk + klev*(1)] < (*yrecldp).rlmin || za[jk]
      < (*yrecldp).ramin) {

      // Evaporate small cloud liquid water amounts
      zlneg[jk + klev*(0)] = zlneg[jk + klev*(0)] + zqx[jk + klev*(0)];
      zqadj = zqx[jk + klev*(0)]*zqtmst;
      tendency_loc_q[jl + klon*(jk + klev*(ibl))] =
        tendency_loc_q[jl + klon*(jk + klev*(ibl))] + zqadj;
      tendency_loc_t[jl + klon*(jk + klev*(ibl))] =
        tendency_loc_t[jl + klon*(jk + klev*(ibl))] - ralvdcp*zqadj;
      zqx[jk + klev*(4)] = zqx[jk + klev*(4)] + zqx[jk + klev*(0)];
      zqx[jk + klev*(0)] = (double) 0.0;

      // Evaporate small cloud ice water amounts
      zlneg[jk + klev*(1)] = zlneg[jk + klev*(1)] + zqx[jk + klev*(1)];
      zqadj = zqx[jk + klev*(1)]*zqtmst;
      tendency_loc_q[jl + klon*(jk + klev*(ibl))] =
        tendency_loc_q[jl + klon*(jk + klev*(ibl))] + zqadj;
      tendency_loc_t[jl + klon*(jk + klev*(ibl))] =
        tendency_loc_t[jl + klon*(jk + klev*(ibl))] - ralsdcp*zqadj;
      zqx[jk + klev*(4)] = zqx[jk + klev*(4)] + zqx[jk + klev*(1)];
      zqx[jk + klev*(1)] = (double) 0.0;

      // Set cloud cover to zero
      za[jk] = (double) 0.0;

    }
  }

  // ---------------------------------
  // Tidy up small CLV variables
  // ---------------------------------
  //DIR$ IVDEP
  for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {
    //DIR$ IVDEP
    for (jk = 0; jk <= klev + -1; jk += 1) {
      //DIR$ IVDEP
      if (zqx[jk + klev*jm] < (*yrecldp).rlmin) {
        zlneg[jk + klev*jm] = zlneg[jk + klev*jm] + zqx[jk + klev*jm];
        zqadj = zqx[jk + klev*jm]*zqtmst;
        tendency_loc_q[jl + klon*(jk + klev*(ibl))] =
          tendency_loc_q[jl + klon*(jk + klev*(ibl))] + zqadj;
        if (iphase[jm] == 1) {
          tendency_loc_t[jl + klon*(jk + klev*(ibl))] =
            tendency_loc_t[jl + klon*(jk + klev*(ibl))] - ralvdcp*zqadj;
        }
        if (iphase[jm] == 2) {
          tendency_loc_t[jl + klon*(jk + klev*(ibl))] =
            tendency_loc_t[jl + klon*(jk + klev*(ibl))] - ralsdcp*zqadj;
        }
        zqx[jk + klev*(4)] = zqx[jk + klev*(4)] + zqx[jk + klev*jm];
        zqx[jk + klev*jm] = (double) 0.0;
      }
    }
  }


  // ------------------------------
  // Define saturation values
  // ------------------------------
  for (jk = 0; jk <= klev + -1; jk += 1) {
    //----------------------------------------
    // old *diagnostic* mixed phase saturation
    //----------------------------------------
    zfoealfa[jk] = ((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))));
    zfoeewmt[jk] =
      fmin(((double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))) / pap[jl + klon*(jk + klev*(ibl))], (double) 0.5);
    zqsmix[jk] = zfoeewmt[jk];
    zqsmix[jk] = zqsmix[jk] / ((double) 1.0 - retv*zqsmix[jk]);

    //---------------------------------------------
    // ice saturation T<273K
    // liquid water saturation for T>273K
    //---------------------------------------------
    zalfa = ((double)(fmax(0.0, copysign(1.0, ztp1[jk] - rtt))));
    zfoeew[jk] = fmin((zalfa*((double)(r2es*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)))) + ((double) 1.0 - zalfa)*((double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))) / pap[jl +
      klon*(jk + klev*(ibl))], (double) 0.5);
    zfoeew[jk] = fmin((double) 0.5, zfoeew[jk]);
    zqsice[jk] = zfoeew[jk] / ((double) 1.0 - retv*zfoeew[jk]);

    //----------------------------------
    // liquid water saturation
    //----------------------------------
    zfoeeliqt[jk] =
      fmin(((double)(r2es*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)))) / pap[jl + klon*(jk + klev*(ibl))], (double) 0.5);
    zqsliq[jk] = zfoeeliqt[jk];
    zqsliq[jk] = zqsliq[jk] / ((double) 1.0 - retv*zqsliq[jk]);

    //   !----------------------------------
    //   ! ice water saturation
    //   !----------------------------------
    //   ZFOEEICET(JL,JK)=MIN(((double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
    //   ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
    //   ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-RETV*ZQSICE(JL,JK))

  }

  for (jk = 0; jk <= klev + -1; jk += 1) {


    //------------------------------------------
    // Ensure cloud fraction is between 0 and 1
    //------------------------------------------
    za[jk] = fmax((double) 0.0, fmin((double) 1.0, za[jk]));

    //-------------------------------------------------------------------
    // Calculate liq/ice fractions (no longer a diagnostic relationship)
    //-------------------------------------------------------------------
    zli[jk] = zqx[jk + klev*(0)] + zqx[jk + klev*(1)];
    if (zli[jk] > (*yrecldp).rlmin) {
      zliqfrac[jk] = zqx[jk + klev*(0)] / zli[jk];
      zicefrac[jk] = (double) 1.0 - zliqfrac[jk];
    } else {
      zliqfrac[jk] = (double) 0.0;
      zicefrac[jk] = (double) 0.0;
    }

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
  ztrpaus = (double) 0.1;
  zpaphd = (double) 1.0 / paph[jl + klon*(klev + (klev + 1)*(ibl))];
  for (jk = 0; jk <= klev - 1 + -1; jk += 1) {
    zsig = pap[jl + klon*(jk + klev*(ibl))]*zpaphd;
    if (zsig > (double) 0.1 && zsig < (double) 0.4 && ztp1[jk] > ztp1[1 + jk]) {
      ztrpaus = zsig;
    }
  }

  //-----------------------------
  // Reset single level variables
  //-----------------------------

  zanewm1 = (double) 0.0;
  zda = (double) 0.0;
  zcovpclr = (double) 0.0;
  zcovpmax = (double) 0.0;
  zcovptot = (double) 0.0;
  zcldtopdist = (double) 0.0;

  //######################################################################
  //           3.       *** PHYSICS ***
  //######################################################################


  //----------------------------------------------------------------------
  //                       START OF VERTICAL LOOP
  //----------------------------------------------------------------------

  for (jk = -1 + (*yrecldp).ncldtop; jk <= klev + -1; jk += 1) {

    //----------------------------------------------------------------------
    // 3.0 INITIALIZE VARIABLES
    //----------------------------------------------------------------------

    //---------------------------------
    // First guess microphysics
    //---------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zqxfg[jm] = zqx[jk + klev*jm];
    }

    //---------------------------------
    // Set KLON arrays to zero
    //---------------------------------

    zlicld = (double) 0.0;
    zrainaut = (double) 0.0;        // currently needed for diags
    zrainacc = (double) 0.0;        // currently needed for diags
    zsnowaut = (double) 0.0;        // needed
    zldefr = (double) 0.0;
    zacust = (double) 0.0;        // set later when needed
    zqpretot = (double) 0.0;
    zlfinalsum = (double) 0.0;

    // Required for first guess call
    zlcond1 = (double) 0.0;
    zlcond2 = (double) 0.0;
    zsupsat = (double) 0.0;
    zlevapl = (double) 0.0;
    zlevapi = (double) 0.0;

    //-------------------------------------
    // solvers for cloud fraction
    //-------------------------------------
    zsolab = (double) 0.0;
    zsolac = (double) 0.0;

    zicetot = (double) 0.0;

    //------------------------------------------
    // reset matrix so missing pathways are set
    //------------------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      for (jn = 0; jn <= 5 + -1; jn += 1) {
        zsolqb[jn + 5*jm] = (double) 0.0;
        zsolqa[jn + 5*jm] = (double) 0.0;
      }
    }

    //----------------------------------
    // reset new microphysics variables
    //----------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zfallsrce[jm] = (double) 0.0;
      zfallsink[jm] = (double) 0.0;
      zconvsrce[jm] = (double) 0.0;
      zconvsink[jm] = (double) 0.0;
      zpsupsatsrce[jm] = (double) 0.0;
      zratio[jm] = (double) 0.0;
    }


    //-------------------------
    // derived variables needed
    //-------------------------

    zdp = paph[jl + klon*(1 + jk + (klev + 1)*(ibl))] - paph[jl +
      klon*(jk + (klev + 1)*(ibl))];        // dp
    zgdp = rg / zdp;        // g/dp
    zrho = pap[jl + klon*(jk + klev*(ibl))] / (rd*ztp1[jk]);        // p/RT air density

    zdtgdp = ptsphy*zgdp;        // dt g/dp
    zrdtgdp = zdp*((double) 1.0 / (ptsphy*rg));        // 1/(dt g/dp)

    if (jk + 1 > 1) {
      zdtgdpf = ptsphy*rg / (pap[jl + klon*(jk + klev*(ibl))] - pap[jl +
         klon*(-1 + jk + klev*(ibl))]);
    }

    //------------------------------------
    // Calculate dqs/dT correction factor
    //------------------------------------
    // Reminder: RETV=RV/RD-1

    // liquid
    zfacw = r5les / (pow((ztp1[jk] - r4les), 2));
    zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeeliqt[jk]);
    zdqsliqdt = zfacw*zcor*zqsliq[jk];
    zcorqsliq = (double) 1.0 + ralvdcp*zdqsliqdt;

    // ice
    zfaci = r5ies / (pow((ztp1[jk] - r4ies), 2));
    zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeew[jk]);
    zdqsicedt = zfaci*zcor*zqsice[jk];
    zcorqsice = (double) 1.0 + ralsdcp*zdqsicedt;

    // diagnostic mixed
    zalfaw = zfoealfa[jk];
    zalfawm = zalfaw;
    zfac = zalfaw*zfacw + ((double) 1.0 - zalfaw)*zfaci;
    zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeewmt[jk]);
    zdqsmixdt = zfac*zcor*zqsmix[jk];
    zcorqsmix = (double) 1.0 + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zdqsmixdt;

    // evaporation/sublimation limits
    zevaplimmix =
      fmax((zqsmix[jk] - zqx[jk + klev*(4)]) / zcorqsmix, (double) 0.0);
    zevaplimliq =
      fmax((zqsliq[jk] - zqx[jk + klev*(4)]) / zcorqsliq, (double) 0.0);
    zevaplimice =
      fmax((zqsice[jk] - zqx[jk + klev*(4)]) / zcorqsice, (double) 0.0);

    //--------------------------------
    // in-cloud consensate amount
    //--------------------------------
    ztmpa = (double) 1.0 / fmax(za[jk], zepsec);
    zliqcld = zqx[jk + klev*(0)]*ztmpa;
    zicecld = zqx[jk + klev*(1)]*ztmpa;
    zlicld = zliqcld + zicecld;


    //------------------------------------------------
    // Evaporate very small amounts of liquid and ice
    //------------------------------------------------

    if (zqx[jk + klev*(0)] < (*yrecldp).rlmin) {
      zsolqa[4 + 5*(0)] = zqx[jk + klev*(0)];
      zsolqa[0 + 5*(4)] = -zqx[jk + klev*(0)];
    }

    if (zqx[jk + klev*(1)] < (*yrecldp).rlmin) {
      zsolqa[4 + 5*(1)] = zqx[jk + klev*(1)];
      zsolqa[1 + 5*(4)] = -zqx[jk + klev*(1)];
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
    zfokoop = ((double)(fmin(rkoop1 - rkoop2*ztp1[jk], (double)(r2es*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)))*1.0/(double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))));

    if (ztp1[jk] >= rtt || (*yrecldp).nssopt == 0) {
      zfac = (double) 1.0;
      zfaci = (double) 1.0;
    } else {
      zfac = za[jk] + zfokoop*((double) 1.0 - za[jk]);
      zfaci = ptsphy / (*yrecldp).rkooptau;
    }

    //-------------------------------------------------------------------
    // 3.1.2 Calculate supersaturation wrt Koop including dqs/dT
    //       correction factor
    // [#Note: QSICE or QSLIQ]
    //-------------------------------------------------------------------

    // Calculate supersaturation to add to cloud
    if (za[jk] > (double) 1.0 - (*yrecldp).ramin) {
      zsupsat =
        fmax((zqx[jk + klev*(4)] - zfac*zqsice[jk]) / zcorqsice, (double) 0.0);
    } else {
      // Calculate environmental humidity supersaturation
      zqp1env = (zqx[jk + klev*(4)] - za[jk]*zqsice[jk]) / fmax((double) 1.0 -
        za[jk], zepsilon);
      //& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
      zsupsat = fmax(((double) 1.0 - za[jk])*(zqp1env - zfac*zqsice[jk]) / zcorqsice,
        (double) 0.0);
    }

    //-------------------------------------------------------------------
    // Here the supersaturation is turned into liquid water
    // However, if the temperature is below the threshold for homogeneous
    // freezing then the supersaturation is turned instantly to ice.
    //--------------------------------------------------------------------

    if (zsupsat > zepsec) {

      if (ztp1[jk] > (*yrecldp).rthomo) {
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
      zsolac = ((double) 1.0 - za[jk])*zfaci;

    }

    //-------------------------------------------------------
    // 3.1.3 Include supersaturation from previous timestep
    // (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
    //-------------------------------------------------------
    if (psupsat[jl + klon*(jk + klev*(ibl))] > zepsec) {
      if (ztp1[jk] > (*yrecldp).rthomo) {
        // Turn supersaturation into liquid water
        zsolqa[0 + 5*(0)] =
          zsolqa[0 + 5*(0)] + psupsat[jl + klon*(jk + klev*(ibl))];
        zpsupsatsrce[0] = psupsat[jl + klon*(jk + klev*(ibl))];
        // Add liquid to first guess for deposition term
        zqxfg[0] = zqxfg[0] + psupsat[jl + klon*(jk + klev*(ibl))];
        // Store cloud budget diagnostics if required
      } else {
        // Turn supersaturation into ice water
        zsolqa[1 + 5*(1)] =
          zsolqa[1 + 5*(1)] + psupsat[jl + klon*(jk + klev*(ibl))];
        zpsupsatsrce[1] = psupsat[jl + klon*(jk + klev*(ibl))];
        // Add ice to first guess for deposition term
        zqxfg[1] = zqxfg[1] + psupsat[jl + klon*(jk + klev*(ibl))];
        // Store cloud budget diagnostics if required
      }

      // Increase cloud amount using RKOOPTAU timescale
      zsolac = ((double) 1.0 - za[jk])*zfaci;
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


      plude[jl + klon*(jk + klev*(ibl))] =
        plude[jl + klon*(jk + klev*(ibl))]*zdtgdp;

      if (/*ldcum[jl + klon*(ibl)] &&*/ plude[jl + klon*(jk + klev*(ibl
        ))] > (*yrecldp).rlmin && plu[jl + klon*(1 + jk + klev*(ibl))] >
        zepsec) {

        zsolac = zsolac + plude[jl + klon*(jk + klev*(ibl))] / plu[jl +
          klon*(1 + jk + klev*(ibl))];
        // *diagnostic temperature split*
        zalfaw = zfoealfa[jk];
        zconvsrce[0] = zalfaw*plude[jl + klon*(jk + klev*(ibl))];
        zconvsrce[1] =
          ((double) 1.0 - zalfaw)*plude[jl + klon*(jk + klev*(ibl))];
        zsolqa[0 + 5*(0)] = zsolqa[0 + 5*(0)] + zconvsrce[0];
        zsolqa[1 + 5*(1)] = zsolqa[1 + 5*(1)] + zconvsrce[1];

      } else {

        plude[jl + klon*(jk + klev*(ibl))] = (double) 0.0;

      }
      // *convective snow detrainment source
      //if (ldcum[jl + klon*(ibl)]) {
        zsolqa[3 + 5*(3)] = zsolqa[3 + 5*(3)] + psnde[jl +
          klon*(jk + klev*(ibl))]*zdtgdp;
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

      zmf = fmax((double) 0.0, (pmfu[jl + klon*(jk + klev*(ibl))] + pmfd[-1 +
         jl + klon*(jk + klev*(ibl))])*zdtgdp);
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

      zdtdp = zrdcp*(double) 0.5*(ztp1[-1 + jk] + ztp1[jk]) / paph[jl + klon*(jk +
         (klev + 1)*(ibl))];
      zdtforc = zdtdp*(pap[jl + klon*(jk + klev*(ibl))] - pap[jl +
        klon*(-1 + jk + klev*(ibl))]);
      //[#Note: Diagnostic mixed phase should be replaced below]
      zdqs = zanewm1*zdtforc*zdqsmixdt;

      for (jm = 0; jm <= 5 + -1; jm += 1) {
        if (!llfall[jm] && iphase[jm] > 0) {
          zlfinal = fmax((double) 0.0, zlcust[jm] - zdqs);              //lim to zero
          // no supersaturation allowed incloud ---V
          zevap = fmin((zlcust[jm] - zlfinal), zevaplimmix);
          //          ZEVAP=0.0_JPRB
          zlfinal = zlcust[jm] - zevap;
          zlfinalsum = zlfinalsum + zlfinal;              // sum

          zsolqa[jm + 5*jm] = zsolqa[jm + 5*jm] + zlcust[jm];              // whole sum
          zsolqa[4 + 5*jm] = zsolqa[4 + 5*jm] + zevap;
          zsolqa[jm + 5*(4)] = zsolqa[jm + 5*(4)] - zevap;
        }
      }

      //  Reset the cloud contribution if no cloud water survives to this level:
      if (zlfinalsum < zepsec) {
        zacust = (double) 0.0;
      }
      zsolac = zsolac + zacust;

    }
    // on  JK>NCLDTOP

    //---------------------------------------------------------------------
    // Subsidence sink of cloud to the layer below
    // (Implicit - re. CFL limit on convective mass flux)
    //---------------------------------------------------------------------


    if (jk + 1 < klev) {

      zmfdn = fmax((double) 0.0, (pmfu[jl + klon*(1 + jk + klev*(ibl))] +
        pmfd[jl + klon*(1 + jk + klev*(ibl))])*zdtgdp);

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
    zldifdt = (*yrecldp).rcldiff*ptsphy;        //original version
    //Increase by factor of 5 for convective points
    if (ktype[jl + klon*(ibl)] > 0 && plude[jl + klon*(jk + klev*(
      ibl))] > zepsec) {
      zldifdt = (*yrecldp).rcldiff_convi*zldifdt;
    }

    // At the moment, works on mixed RH profile and partitioned ice/liq fraction
    // so that it is similar to previous scheme
    // Should apply RHw for liquid cloud and RHi for ice cloud separately
    if (zli[jk] > zepsec) {
      // Calculate environmental humidity
      //      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
      //    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
      //      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
      ze = zldifdt*fmax(zqsmix[jk] - zqx[jk + klev*(4)], (double) 0.0);
      zleros = za[jk]*ze;
      zleros = fmin(zleros, zevaplimmix);
      zleros = fmin(zleros, zli[jk]);
      zaeros = zleros / zlicld;          //if linear term

      // Erosion is -ve LINEAR in L,A
      zsolac = zsolac - zaeros;          //linear

      zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] + zliqfrac[jk]*zleros;
      zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] - zliqfrac[jk]*zleros;
      zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] + zicefrac[jk]*zleros;
      zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] - zicefrac[jk]*zleros;

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

    zdtdp = zrdcp*ztp1[jk] / pap[jl + klon*(jk + klev*(ibl))];
    zdpmxdt = zdp*zqtmst;
    zmfdn = (double) 0.0;
    if (jk + 1 < klev) {
      zmfdn = pmfu[jl + klon*(1 + jk + klev*(ibl))] + pmfd[jl + klon*(1
        + jk + klev*(ibl))];
    }
    zwtot = pvervel[jl + klon*(jk + klev*(ibl))] + (double) 0.5*rg*(pmfu[
       jl + klon*(jk + klev*(ibl))] + pmfd[jl + klon*(jk + klev*(ibl))]
      + zmfdn);
    zwtot = fmin(zdpmxdt, fmax(-zdpmxdt, zwtot));
    zzzdt = phrsw[jl + klon*(jk + klev*(ibl))] + phrlw[jl + klon*(jk +
      klev*(ibl))];
    zdtdiab = fmin(zdpmxdt*zdtdp, fmax(-zdpmxdt*zdtdp, zzzdt))*ptsphy + ralfdcp*zldefr;
    // Note: ZLDEFR should be set to the difference between the mixed phase functions
    // in the convection and cloud scheme, but this is not calculated, so is zero and
    // the functions must be the same
    zdtforc = zdtdp*zwtot*ptsphy + zdtdiab;
    zqold = zqsmix[jk];
    ztold = ztp1[jk];
    ztp1[jk] = ztp1[jk] + zdtforc;
    ztp1[jk] = fmax(ztp1[jk], (double) 160.0);
    llflag = true;

    // Formerly a call to CUADJTQ(..., ICALL=5)
    zqp = (double) 1.0 / pap[jl + klon*(jk + klev*(ibl))];
    zqsat = ((double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies)))))*zqp;
    zqsat = fmin((double) 0.5, zqsat);
    zcor = (double) 1.0 / ((double) 1.0 - retv*zqsat);
    zqsat = zqsat*zcor;
    zcond = (zqsmix[jk] - zqsat) / ((double) 1.0 + zqsat*zcor*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk] - r4ies, 2)))));
    ztp1[jk] = ztp1[jk] + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zcond;
    zqsmix[jk] = zqsmix[jk] - zcond;
    zqsat = ((double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies)))))*zqp;
    zqsat = fmin((double) 0.5, zqsat);
    zcor = (double) 1.0 / ((double) 1.0 - retv*zqsat);
    zqsat = zqsat*zcor;
    zcond1 = (zqsmix[jk] - zqsat) / ((double) 1.0 + zqsat*zcor*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk] - r4ies, 2)))));
    ztp1[jk] = ztp1[jk] + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp))*zcond1;
    zqsmix[jk] = zqsmix[jk] - zcond1;

    zdqs = zqsmix[jk] - zqold;
    zqsmix[jk] = zqold;
    ztp1[jk] = ztold;

    //----------------------------------------------------------------------
    // 3.4a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
    // ----------------------------------------------------------------------
    // Erosion term is LINEAR in L
    // Changed to be uniform distribution in cloud region


    // Previous function based on DELTA DISTRIBUTION in cloud:
    if (zdqs > (double) 0.0) {
      //    If subsidence evaporation term is turned off, then need to use updated
      //    liquid and cloud here?
      //    ZLEVAP = MAX(ZA(JL,JK)+ZACUST(JL),1.0_JPRB)*MIN(ZDQS(JL),ZLICLD(JL)+ZLFINALSUM(JL))
      zlevap = za[jk]*fmin(zdqs, zlicld);
      zlevap = fmin(zlevap, zevaplimmix);
      zlevap = fmin(zlevap, fmax(zqsmix[jk] - zqx[jk + klev*(4)], (double) 0.0));

      // For first guess call
      zlevapl = zliqfrac[jk]*zlevap;
      zlevapi = zicefrac[jk]*zlevap;

      zsolqa[4 + 5*(0)] = zsolqa[4 + 5*(0)] + zliqfrac[jk]*zlevap;
      zsolqa[0 + 5*(4)] = zsolqa[0 + 5*(4)] - zliqfrac[jk]*zlevap;

      zsolqa[4 + 5*(1)] = zsolqa[4 + 5*(1)] + zicefrac[jk]*zlevap;
      zsolqa[1 + 5*(4)] = zsolqa[1 + 5*(4)] - zicefrac[jk]*zlevap;

    }


    //----------------------------------------------------------------------
    // 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
    //----------------------------------------------------------------------
    // (1) Increase of cloud water in existing clouds
    if (za[jk] > zepsec && zdqs <= -(*yrecldp).rlmin) {

      zlcond1 = fmax(-zdqs, (double) 0.0);          //new limiter

      //old limiter (significantly improves upper tropospheric humidity rms)
      if (za[jk] > (double) 0.99) {
        zcor = (double) 1.0 / ((double) 1.0 - retv*zqsmix[jk]);
        zcdmax = (zqx[jk + klev*(4)] - zqsmix[jk]) / ((double) 1.0 +
          zcor*zqsmix[jk]*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk] - r4ies, 2)))));
      } else {
        zcdmax = (zqx[jk + klev*(4)] - za[jk]*zqsmix[jk]) / za[jk];
      }
      zlcond1 = fmax(fmin(zlcond1, zcdmax), (double) 0.0);
      // end old limiter

      zlcond1 = za[jk]*zlcond1;
      if (zlcond1 < (*yrecldp).rlmin) {
        zlcond1 = (double) 0.0;
      }

      //-------------------------------------------------------------------------
      // All increase goes into liquid unless so cold cloud homogeneously freezes
      // Include new liquid formation in first guess value, otherwise liquid
      // remains at cold temperatures until next timestep.
      //-------------------------------------------------------------------------
      if (ztp1[jk] > (*yrecldp).rthomo) {
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


    if (zdqs <= -(*yrecldp).rlmin && za[jk] < (double) 1.0 - zepsec) {

      //---------------------------
      // Critical relative humidity
      //---------------------------
      zrhc = (*yrecldp).ramid;
      zsigk = pap[jl + klon*(jk + klev*(ibl))] / paph[jl + klon*(klev +
        (klev + 1)*(ibl))];
      // Increase RHcrit to 1.0 towards the surface (eta>0.8)
      if (zsigk > (double) 0.8) {
        zrhc = (*yrecldp).ramid + ((double) 1.0 - (*yrecldp).ramid)*(pow(((zsigk -
          (double) 0.8) / (double) 0.2), 2));
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
        zqe = (zqx[jk + klev*(4)] - za[jk]*zqsice[jk]) / fmax(zepsec, (double) 1.0
           - za[jk]);
        zqe = fmax((double) 0.0, zqe);
      } else if ((*yrecldp).nssopt == 1) {
        // Tompkins
        zqe = (zqx[jk + klev*(4)] - za[jk]*zqsice[jk]) / fmax(zepsec, (double) 1.0
           - za[jk]);
        zqe = fmax((double) 0.0, zqe);
      } else if ((*yrecldp).nssopt == 2) {
        // Lohmann and Karcher
        zqe = zqx[jk + klev*(4)];
      } else if ((*yrecldp).nssopt == 3) {
        // Gierens
        zqe = zqx[jk + klev*(4)] + zli[jk];
      }

      if (ztp1[jk] >= rtt || (*yrecldp).nssopt == 0) {
        // No ice supersaturation allowed
        zfac = (double) 1.0;
      } else {
        // Ice supersaturation
        zfac = zfokoop;
      }

      if (zqe >= zrhc*zqsice[jk]*zfac && zqe < zqsice[jk]*zfac) {
        // note: not **2 on 1-a term if ZQE is used.
        // Added correction term ZFAC to numerator 15/03/2010
        zacond = -((double) 1.0 - za[jk])*zfac*zdqs / fmax((double)
          2.0*(zfac*zqsice[jk] - zqe), zepsec);

        zacond = fmin(zacond, (double) 1.0 - za[jk]);            //PUT THE LIMITER BACK

        // Linear term:
        // Added correction term ZFAC 15/03/2010
        zlcond2 = -zfac*zdqs*(double) 0.5*zacond;            //mine linear

        // new limiter formulation
        zzdl =
          (double) 2.0*(zfac*zqsice[jk] - zqe) / fmax(zepsec, (double) 1.0 - za[jk]);
        // Added correction term ZFAC 15/03/2010
        if (zfac*zdqs < -zzdl) {
          // ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
          zlcondlim = (za[jk] - (double) 1.0)*zfac*zdqs - zfac*zqsice[jk] + zqx[jk +
            klev*(4)];
          zlcond2 = fmin(zlcond2, zlcondlim);
        }
        zlcond2 = fmax(zlcond2, (double) 0.0);

        if (zlcond2 < (*yrecldp).rlmin || ((double) 1.0 - za[jk]) < zepsec) {
          zlcond2 = (double) 0.0;
          zacond = (double) 0.0;
        }
        if (zlcond2 == (double) 0.0) {
          zacond = (double) 0.0;
        }

        // Large-scale generation is LINEAR in A and LINEAR in L
        zsolac = zsolac + zacond;            //linear

        //------------------------------------------------------------------------
        // All increase goes into liquid unless so cold cloud homogeneously freezes
        // Include new liquid formation in first guess value, otherwise liquid
        // remains at cold temperatures until next timestep.
        //------------------------------------------------------------------------
        if (ztp1[jk] > (*yrecldp).rthomo) {
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
      // ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
      //--------------------------------------------------------------

      if (za[-1 + jk] < (*yrecldp).rcldtopcf && za[jk] >= (*yrecldp).rcldtopcf) {
        zcldtopdist = (double) 0.0;
      } else {
        zcldtopdist = zcldtopdist + zdp / (zrho*rg);
      }

      //--------------------------------------------------------------
      // only treat depositional growth if liquid present. due to fact
      // that can not model ice growth from vapour without additional
      // in-cloud water vapour variable
      //--------------------------------------------------------------
      if (ztp1[jk] < rtt && zqxfg[0] > (*yrecldp).rlmin) {
        // T<273K

        zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))*rv / rd;
        zvpliq = zvpice*zfokoop;
        zicenuclei = (double) 1000.0*exp((double) 12.96*(zvpliq - zvpice) / zvpliq -
          (double) 0.639);

        //------------------------------------------------
        //   2.4e-2 is conductivity of air
        //   8.8 = 700**1/3 = density of ice to the third
        //------------------------------------------------
        zadd =
          rlstt*(rlstt / (rv*ztp1[jk]) - (double) 1.0) / ((double) 2.4E-2*ztp1[jk]);
        zbdd = rv*ztp1[jk]*pap[jl + klon*(jk + klev*(ibl))] / ((double)
          2.21*zvpice);
        zcvds = (double) 7.8*(pow((zicenuclei / zrho), (double) 0.666))*(zvpliq -
          zvpice) / ((double) 8.87*(zadd + zbdd)*zvpice);

        //-----------------------------------------------------
        // RICEINIT=1.E-12_JPRB is initial mass of ice particle
        //-----------------------------------------------------
        zice0 = fmax(zicecld, zicenuclei*(*yrecldp).riceinit / zrho);

        //------------------
        // new value of ice:
        //------------------
        zinew = pow(((double) 0.666*zcvds*ptsphy + (pow(zice0, (double) 0.666))),
          (double) 1.5);

        //---------------------------
        // grid-mean deposition rate:
        //---------------------------
        zdepos = fmax(za[jk]*(zinew - zice0), (double) 0.0);

        //--------------------------------------------------------------------
        // Limit deposition to liquid water amount
        // If liquid is all frozen, ice would use up reservoir of water
        // vapour in excess of ice saturation mixing ratio - However this
        // can not be represented without a in-cloud humidity variable. Using
        // the grid-mean humidity would imply a large artificial horizontal
        // flux from the clear sky to the cloudy area. We thus rely on the
        // supersaturation check to clean up any remaining supersaturation
        //--------------------------------------------------------------------
        zdepos = fmin(zdepos, zqxfg[0]);            // limit to liquid water amount

        //--------------------------------------------------------------------
        // At top of cloud, reduce deposition rate near cloud top to account for
        // small scale turbulent processes, limited ice nucleation and ice fallout
        //--------------------------------------------------------------------
        //      ZDEPOS = ZDEPOS*MIN(RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH,1.0_JPRB)
        // Change to include dependence on ice nuclei concentration
        // to increase deposition rate with decreasing temperatures
        zinfactor = fmin(zicenuclei / (double) 15000., (double) 1.0);
        zdepos = zdepos*fmin(zinfactor + ((double) 1.0 - zinfactor)*((*yrecldp)
          .rdepliqrefrate + zcldtopdist / (*yrecldp).rdepliqrefdepth), (double) 1.0);

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
      // ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
      //--------------------------------------------------------------

      if (za[-1 + jk] < (*yrecldp).rcldtopcf && za[jk] >= (*yrecldp).rcldtopcf) {
        zcldtopdist = (double) 0.0;
      } else {
        zcldtopdist = zcldtopdist + zdp / (zrho*rg);
      }

      //--------------------------------------------------------------
      // only treat depositional growth if liquid present. due to fact
      // that can not model ice growth from vapour without additional
      // in-cloud water vapour variable
      //--------------------------------------------------------------
      if (ztp1[jk] < rtt && zqxfg[0] > (*yrecldp).rlmin) {
        // T<273K

        zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))*rv / rd;
        zvpliq = zvpice*zfokoop;
        zicenuclei = (double) 1000.0*exp((double) 12.96*(zvpliq - zvpice) / zvpliq -
          (double) 0.639);

        //-----------------------------------------------------
        // RICEINIT=1.E-12_JPRB is initial mass of ice particle
        //-----------------------------------------------------
        zice0 = fmax(zicecld, zicenuclei*(*yrecldp).riceinit / zrho);

        // Particle size distribution
        ztcg = (double) 1.0;
        zfacx1i = (double) 1.0;

        zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk] +
          pap[jl + klon*(jk + klev*(ibl))]*(*yrecldp).rcl_apb3*(pow(ztp1[jk],
           (double) 3.));
        zcorrfac = pow(((double) 1.0 / zrho), (double) 0.5);
        zcorrfac2 = (pow((ztp1[jk] / (double) 273.0), (double) 1.5))*((double) 393.0 /
          (ztp1[jk] + (double) 120.0));

        zpr02 = zrho*zice0*(*yrecldp).rcl_const1i / (ztcg*zfacx1i);

        zterm1 = (zvpliq - zvpice)*(pow(ztp1[jk], (double) 2.0))
          *zvpice*zcorrfac2*ztcg*(*yrecldp).rcl_const2i*zfacx1i / (zrho*zaplusb*zvpice)
          ;
        zterm2 = (double) 0.65*(*yrecldp).rcl_const6i*(pow(zpr02, (*yrecldp)
          .rcl_const4i)) + (*yrecldp).rcl_const3i*(pow(zcorrfac, (double) 0.5))
          *(pow(zrho, (double) 0.5))*(pow(zpr02, (*yrecldp).rcl_const5i)) /
          (pow(zcorrfac2, (double) 0.5));

        zdepos = fmax(za[jk]*zterm1*zterm2*ptsphy, (double) 0.0);

        //--------------------------------------------------------------------
        // Limit deposition to liquid water amount
        // If liquid is all frozen, ice would use up reservoir of water
        // vapour in excess of ice saturation mixing ratio - However this
        // can not be represented without a in-cloud humidity variable. Using
        // the grid-mean humidity would imply a large artificial horizontal
        // flux from the clear sky to the cloudy area. We thus rely on the
        // supersaturation check to clean up any remaining supersaturation
        //--------------------------------------------------------------------
        zdepos = fmin(zdepos, zqxfg[0]);            // limit to liquid water amount

        //--------------------------------------------------------------------
        // At top of cloud, reduce deposition rate near cloud top to account for
        // small scale turbulent processes, limited ice nucleation and ice fallout
        //--------------------------------------------------------------------
        // Change to include dependence on ice nuclei concentration
        // to increase deposition rate with decreasing temperatures
        zinfactor = fmin(zicenuclei / (double) 15000., (double) 1.0);
        zdepos = zdepos*fmin(zinfactor + ((double) 1.0 - zinfactor)*((*yrecldp)
          .rdepliqrefrate + zcldtopdist / (*yrecldp).rdepliqrefdepth), (double) 1.0);

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
    ztmpa = (double) 1.0 / fmax(za[jk], zepsec);
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
          zfallsrce[jm] = zpfplsx[jk + (klev + 1)*jm]*zdtgdp;
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
          zre_ice = pre_ice[jl + klon*(jk + klev*(ibl))];
          // The exponent value is from
          // Morrison et al. JAS 2005 Appendix
          zvqx[1] = (double) 0.002*(pow(zre_ice, (double) 1.0));
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
      zcovptot = (double) 1.0 - (((double) 1.0 - zcovptot)*((double) 1.0 - fmax(za[jk],
         za[-1 + jk])) / ((double) 1.0 - fmin(za[-1 + jk], (double) 1.0 - (double)
        1.E-06)));
      zcovptot = fmax(zcovptot, (*yrecldp).rcovpmin);
      zcovpclr = fmax((double) 0.0, zcovptot - za[jk]);          // clear sky proportion
      zraincld = zqxfg[2] / zcovptot;
      zsnowcld = zqxfg[3] / zcovptot;
      zcovpmax = fmax(zcovptot, zcovpmax);
    } else {
      zraincld = (double) 0.0;
      zsnowcld = (double) 0.0;
      zcovptot = (double) 0.0;          // no flux - reset cover
      zcovpclr = (double) 0.0;          // reset clear sky proportion
      zcovpmax = (double) 0.0;          // reset max cover for ZZRH calc
    }

    //----------------------------------------------------------------------
    // 4.3a AUTOCONVERSION TO SNOW
    //----------------------------------------------------------------------

    if (ztp1[jk] <= rtt) {
      //-----------------------------------------------------
      //     Snow Autoconversion rate follow Lin et al. 1983
      //-----------------------------------------------------
      if (zicecld > zepsec) {

        zzco = ptsphy*(*yrecldp).rsnowlin1*exp((*yrecldp).rsnowlin2*(ztp1[jk] - rtt));

        if ((*yrecldp).laericeauto) {
          zlcrit = picrit_aer[jl + klon*(jk + klev*(ibl))];
          // 0.3 = N**0.333 with N=0.027
          zzco = zzco*(pow(((*yrecldp).rnice / pnice[jl + klon*(jk + klev*(
            ibl))]), (double) 0.333));
        } else {
          zlcrit = (*yrecldp).rlcritsnow;
        }

        zsnowaut = zzco*((double) 1.0 - exp(-(pow((zicecld / zlcrit), 2))));
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
          zlcrit = plcrit_aer[jl + klon*(jk + klev*(ibl))];
          // 0.3 = N**0.333 with N=125 cm-3
          zzco = zzco*(pow(((*yrecldp).rccn / pccn[jl + klon*(jk + klev*(ibl)
            )]), (double) 0.333));
        } else {
          // Modify autoconversion threshold dependent on:
          //  land (polluted, high CCN, smaller droplets, higher threshold)
          //  sea  (clean, low CCN, larger droplets, lower threshold)
          if (plsm[jl + klon*(ibl)] > (double) 0.5) {
            zlcrit = (*yrecldp).rclcrit_land;                // land
          } else {
            zlcrit = (*yrecldp).rclcrit_sea;                // ocean
          }
        }

        //------------------------------------------------------------------
        // Parameters for cloud collection by rain and snow.
        // Note that with new prognostic variable it is now possible
        // to REPLACE this with an explicit collection parametrization
        //------------------------------------------------------------------
        zprecip = (zpfplsx[jk + (klev + 1)*(3)] + zpfplsx[jk + (klev + 1)*(2)
          ]) / fmax(zepsec, zcovptot);
        zcfpr = (double) 1.0 + (*yrecldp).rprc1*sqrt(fmax(zprecip, (double) 0.0));
        //      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
        //       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))

        if ((*yrecldp).laerliqcoll) {
          // 5.0 = N**0.333 with N=125 cm-3
          zcfpr = zcfpr*(pow(((*yrecldp).rccn / pccn[jl + klon*(jk + klev*(
            ibl))]), (double) 0.333));
        }

        zzco = zzco*zcfpr;
        zlcrit = zlcrit / fmax(zcfpr, zepsec);

        if (zliqcld / zlcrit < (double) 20.0) {
          // Security for exp for some compilers
          zrainaut = zzco*((double) 1.0 - exp(-(pow((zliqcld / zlcrit), 2))));
        } else {
          zrainaut = zzco;
        }

        // rain freezes instantly
        if (ztp1[jk] <= rtt) {
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

        if (plsm[jl + klon*(ibl)] > (double) 0.5) {
          // land
          zconst = (*yrecldp).rcl_kk_cloud_num_land;
          zlcrit = (*yrecldp).rclcrit_land;
        } else {
          // ocean
          zconst = (*yrecldp).rcl_kk_cloud_num_sea;
          zlcrit = (*yrecldp).rclcrit_sea;
        }

        if (zliqcld > zlcrit) {

          zrainaut = (double) 1.5*za[jk]*ptsphy*(*yrecldp).rcl_kkaau*(pow(zliqcld,
            (*yrecldp).rcl_kkbauq))*(pow(zconst, (*yrecldp).rcl_kkbaun));

          zrainaut = fmin(zrainaut, zqxfg[0]);
          if (zrainaut < zepsec) {
            zrainaut = (double) 0.0;
          }

          zrainacc = (double) 2.0*za[jk]*ptsphy*(*yrecldp)
            .rcl_kkaac*(pow((zliqcld*zraincld), (*yrecldp).rcl_kkbac));

          zrainacc = fmin(zrainacc, zqxfg[0]);
          if (zrainacc < zepsec) {
            zrainacc = (double) 0.0;
          }

        } else {
          zrainaut = (double) 0.0;
          zrainacc = (double) 0.0;
        }

        // If temperature < 0, then autoconversion produces snow rather than rain
        // Explicit
        if (ztp1[jk] <= rtt) {
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

      if (ztp1[jk] <= rtt && zliqcld > zepsec) {

        // Fallspeed air density correction
        zfallcorr = pow(((*yrecldp).rdensref / zrho), (double) 0.4);

        //------------------------------------------------------------------
        // Riming of snow by cloud water - implicit in lwc
        //------------------------------------------------------------------
        if (zsnowcld > zepsec && zcovptot > (double) 0.01) {

          // Calculate riming term
          // Factor of liq water taken out because implicit
          zsnowrime = (double) 0.3*zcovptot*ptsphy*(*yrecldp)
            .rcl_const7s*zfallcorr*(pow((zrho*zsnowcld*(*yrecldp).rcl_const1s),
            (*yrecldp).rcl_const8s));

          // Limit snow riming term
          zsnowrime = fmin(zsnowrime, (double) 1.0);

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
    zmeltmax = (double) 0.0;

    // If there are frozen hydrometeors present and dry-bulb temperature > 0degC
    if (zicetot > zepsec && ztp1[jk] > rtt) {

      // Calculate subsaturation
      zsubsat = fmax(zqsice[jk] - zqx[jk + klev*(4)], (double) 0.0);

      // Calculate difference between dry-bulb (ZTP1) and the temperature
      // at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
      // Melting only occurs if the wet-bulb temperature >0
      // i.e. warming of ice particle due to melting > cooling
      // due to evaporation.
      ztdmtw0 = ztp1[jk] - rtt - zsubsat*(ztw1 + ztw2*(pap[jl + klon*(jk +
        klev*(ibl))] - ztw3) - ztw4*(ztp1[jk] - ztw5));
      // Not implicit yet...
      // Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
      zcons1 = fabs(ptsphy*((double) 1.0 + (double) 0.5*ztdmtw0) / (*yrecldp).rtaumel);
      zmeltmax = fmax(ztdmtw0*zcons1*zrldcp, (double) 0.0);
    }

    // Loop over frozen hydrometeors (ice, snow)
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      if (iphase[jm] == 2) {
        jn = imelt[jm];
        if (zmeltmax > zepsec && zicetot > zepsec) {
          // Apply melting in same proportion as frozen hydrometeor fractions
          zalfa = zqxfg[jm] / zicetot;
          zmelt = fmin(zqxfg[jm], zalfa*zmeltmax);
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
    if (zqx[jk + klev*(2)] > zepsec) {

      if (ztp1[jk] <= rtt && ztp1[-1 + jk] > rtt) {
        // Base of melting layer/top of refreezing layer so
        // store rain/snow fraction for precip type diagnosis
        // If mostly rain, then supercooled rain slow to freeze
        // otherwise faster to freeze (snow or ice pellets)
        zqpretot = fmax(zqx[jk + klev*(3)] + zqx[jk + klev*(2)], zepsec);
        prainfrac_toprfz[jl + klon*(ibl)] =
          zqx[jk + klev*(2)] / zqpretot;
        if (prainfrac_toprfz[jl + klon*(ibl)] > 0.8) {
          llrainliq = true;
        } else {
          llrainliq = false;
        }
      }

      // If temperature less than zero
      if (ztp1[jk] < rtt) {

        if (prainfrac_toprfz[jl + klon*(ibl)] > 0.8) {

          // Majority of raindrops completely melted
          // Refreezing is by slow heterogeneous freezing

          // Slope of rain particle size distribution
          zlambda = pow(((*yrecldp).rcl_fac1 / (zrho*zqx[jk + klev*(2)])),
            (*yrecldp).rcl_fac2);

          // Calculate freezing rate based on Bigg(1953) and Wisner(1972)
          ztemp = (*yrecldp).rcl_fzrab*(ztp1[jk] - rtt);
          zfrz = ptsphy*((*yrecldp).rcl_const5r / zrho)*(exp(ztemp) - (double) 1.)
            *(pow(zlambda, (*yrecldp).rcl_const6r));
          zfrzmax = fmax(zfrz, (double) 0.0);

        } else {

          // Majority of raindrops only partially melted
          // Refreeze with a shorter timescale (reverse of melting...for now)

          zcons1 = fabs(ptsphy*((double) 1.0 + (double) 0.5*(rtt - ztp1[jk])) /
            (*yrecldp).rtaumel);
          zfrzmax = fmax((rtt - ztp1[jk])*zcons1*zrldcp, (double) 0.0);

        }

        if (zfrzmax > zepsec) {
          zfrz = fmin(zqx[jk + klev*(2)], zfrzmax);
          zsolqa[3 + 5*(2)] = zsolqa[3 + 5*(2)] + zfrz;
          zsolqa[2 + 5*(3)] = zsolqa[2 + 5*(3)] - zfrz;
        }
      }

    }


    //----------------------------------------------------------------------
    // 4.4c  FREEZING of LIQUID
    //----------------------------------------------------------------------
    // not implicit yet...
    zfrzmax = fmax(((*yrecldp).rthomo - ztp1[jk])*zrldcp, (double) 0.0);

    jm = 1;
    jn = imelt[-1 + jm];
    if (zfrzmax > zepsec && zqxfg[-1 + jm] > zepsec) {
      zfrz = fmin(zqxfg[-1 + jm], zfrzmax);
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


      zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax /
        fmax(zepsec, (double) 1.0 - za[jk]);
      zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);

      zqe = (zqx[jk + klev*(4)] - za[jk]*zqsliq[jk]) / fmax(zepsec, (double) 1.0 -
         za[jk]);
      //---------------------------------------------
      // humidity in moistest ZCOVPCLR part of domain
      //---------------------------------------------
      zqe = fmax((double) 0.0, fmin(zqe, zqsliq[jk]));
      llo1 = zcovpclr > zepsec && zqxfg[2] > zepsec && zqe < zzrh*zqsliq[jk];

      if (llo1) {
        // note: zpreclr is a rain flux
        zpreclr = zqxfg[2]*zcovpclr / copysign(fmax(fabs(zcovptot*zdtgdp),
          zepsilon), zcovptot*zdtgdp);

        //--------------------------------------
        // actual microphysics formula in zbeta
        //--------------------------------------

        zbeta1 = sqrt(pap[jl + klon*(jk + klev*(ibl))] / paph[jl +
          klon*(klev + (klev + 1)*(ibl))]) / (*yrecldp).rvrfactor*zpreclr /
          fmax(zcovpclr, zepsec);

        zbeta = rg*(*yrecldp).rpecons*(double) 0.5*(pow(zbeta1, (double) 0.5777));

        zdenom = (double) 1.0 + zbeta*ptsphy*zcorqsliq;
        zdpr = zcovpclr*zbeta*(zqsliq[jk] - zqe) / zdenom*zdp*zrg_r;
        zdpevap = zdpr*zdtgdp;

        //---------------------------------------------------------
        // add evaporation term to explicit sink.
        // this has to be explicit since if treated in the implicit
        // term evaporation can not reduce rain to zero and model
        // produces small amounts of rainfall everywhere.
        //---------------------------------------------------------

        // Evaporate rain
        zevap = fmin(zdpevap, zqxfg[2]);

        zsolqa[4 + 5*(2)] = zsolqa[4 + 5*(2)] + zevap;
        zsolqa[2 + 5*(4)] = zsolqa[2 + 5*(4)] - zevap;

        //-------------------------------------------------------------
        // Reduce the total precip coverage proportional to evaporation
        // to mimic the previous scheme which had a diagnostic
        // 2-flux treatment, abandoned due to the new prognostic precip
        //-------------------------------------------------------------
        zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot -
          za[jk])*zevap / zqxfg[2]));

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
      zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax /
        fmax(zepsec, (double) 1.0 - za[jk]);
      zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);

      // Critical relative humidity
      //ZRHC=RAMID
      //ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
      // Increase RHcrit to 1.0 towards the surface (eta>0.8)
      //IF(ZSIGK > 0.8_JPRB) THEN
      //  ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
      //ENDIF
      //ZZRH = MIN(ZRHC,ZZRH)

      // Further limit RH for rain evaporation to 80% (RHcrit in free troposphere)
      zzrh = fmin((double) 0.8, zzrh);

      zqe = fmax((double) 0.0, fmin(zqx[jk + klev*(4)], zqsliq[jk]));

      llo1 = zcovpclr > zepsec && zqxfg[2] > zepsec && zqe < zzrh*zqsliq[jk];

      if (llo1) {

        //-------------------------------------------
        // Abel and Boutle (2012) evaporation
        //-------------------------------------------
        // Calculate local precipitation (kg/kg)
        zpreclr = zqxfg[2] / zcovptot;

        // Fallspeed air density correction
        zfallcorr = pow(((*yrecldp).rdensref / zrho), 0.4);

        // Saturation vapour pressure with respect to liquid phase
        zesatliq = rv / rd*((double)(r2es*exp((r3les*(ztp1[jk] - rtt))/(ztp1[jk] - r4les))));

        // Slope of particle size distribution
        zlambda = pow(((*yrecldp).rcl_fac1 / (zrho*zpreclr)), (*yrecldp).rcl_fac2);            // ZPRECLR=kg/kg

        zevap_denom = (*yrecldp).rcl_cdenom1*zesatliq - (*yrecldp)
          .rcl_cdenom2*ztp1[jk]*zesatliq + (*yrecldp).rcl_cdenom3*(pow(ztp1[jk],
          (double) 3.))*pap[jl + klon*(jk + klev*(ibl))];

        // Temperature dependent conductivity
        zcorr2 = (pow((ztp1[jk] / (double) 273.), (double) 1.5))*(double) 393. /
          (ztp1[jk] + (double) 120.);
        zka = (*yrecldp).rcl_ka273*zcorr2;

        zsubsat = fmax(zzrh*zqsliq[jk] - zqe, (double) 0.0);

        zbeta = ((double) 0.5 / zqsliq[jk])*(pow(ztp1[jk], (double) 2.))
          *zesatliq*(*yrecldp).rcl_const1r*(zcorr2 / zevap_denom)*((double) 0.78 /
          (pow(zlambda, (*yrecldp).rcl_const4r)) + (*yrecldp)
          .rcl_const2r*(pow((zrho*zfallcorr), (double) 0.5)) / ((pow(zcorr2, (double)
          0.5))*(pow(zlambda, (*yrecldp).rcl_const3r))));

        zdenom = (double) 1.0 + zbeta*ptsphy;            //*ZCORQSLIQ(JL)
        zdpevap = zcovpclr*zbeta*ptsphy*zsubsat / zdenom;

        //---------------------------------------------------------
        // Add evaporation term to explicit sink.
        // this has to be explicit since if treated in the implicit
        // term evaporation can not reduce rain to zero and model
        // produces small amounts of rainfall everywhere.
        //---------------------------------------------------------

        // Limit rain evaporation
        zevap = fmin(zdpevap, zqxfg[2]);

        zsolqa[4 + 5*(2)] = zsolqa[4 + 5*(2)] + zevap;
        zsolqa[2 + 5*(4)] = zsolqa[2 + 5*(4)] - zevap;

        //-------------------------------------------------------------
        // Reduce the total precip coverage proportional to evaporation
        // to mimic the previous scheme which had a diagnostic
        // 2-flux treatment, abandoned due to the new prognostic precip
        //-------------------------------------------------------------
        zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot -
          za[jk])*zevap / zqxfg[2]));

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

      zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax /
        fmax(zepsec, (double) 1.0 - za[jk]);
      zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);
      zqe = (zqx[jk + klev*(4)] - za[jk]*zqsice[jk]) / fmax(zepsec, (double) 1.0 -
         za[jk]);

      //---------------------------------------------
      // humidity in moistest ZCOVPCLR part of domain
      //---------------------------------------------
      zqe = fmax((double) 0.0, fmin(zqe, zqsice[jk]));
      llo1 = zcovpclr > zepsec && zqxfg[3] > zepsec && zqe < zzrh*zqsice[jk];

      if (llo1) {
        // note: zpreclr is a rain flux a
        zpreclr = zqxfg[3]*zcovpclr / copysign(fmax(fabs(zcovptot*zdtgdp),
          zepsilon), zcovptot*zdtgdp);

        //--------------------------------------
        // actual microphysics formula in zbeta
        //--------------------------------------

        zbeta1 = sqrt(pap[jl + klon*(jk + klev*(ibl))] / paph[jl +
          klon*(klev + (klev + 1)*(ibl))]) / (*yrecldp).rvrfactor*zpreclr /
          fmax(zcovpclr, zepsec);

        zbeta = rg*(*yrecldp).rpecons*(pow(zbeta1, (double) 0.5777));

        zdenom = (double) 1.0 + zbeta*ptsphy*zcorqsice;
        zdpr = zcovpclr*zbeta*(zqsice[jk] - zqe) / zdenom*zdp*zrg_r;
        zdpevap = zdpr*zdtgdp;

        //---------------------------------------------------------
        // add evaporation term to explicit sink.
        // this has to be explicit since if treated in the implicit
        // term evaporation can not reduce snow to zero and model
        // produces small amounts of snowfall everywhere.
        //---------------------------------------------------------

        // Evaporate snow
        zevap = fmin(zdpevap, zqxfg[3]);

        zsolqa[4 + 5*(3)] = zsolqa[4 + 5*(3)] + zevap;
        zsolqa[3 + 5*(4)] = zsolqa[3 + 5*(4)] - zevap;

        //-------------------------------------------------------------
        // Reduce the total precip coverage proportional to evaporation
        // to mimic the previous scheme which had a diagnostic
        // 2-flux treatment, abandoned due to the new prognostic precip
        //-------------------------------------------------------------
        zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot -
          za[jk])*zevap / zqxfg[3]));

        //Update first guess field
        zqxfg[3] = zqxfg[3] - zevap;

      }
      //---------------------------------------------------------
    } else if (ievapsnow == 2) {



      //-----------------------------------------------------------------------
      // Calculate relative humidity limit for snow evaporation
      //-----------------------------------------------------------------------
      zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax /
        fmax(zepsec, (double) 1.0 - za[jk]);
      zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);
      zqe = (zqx[jk + klev*(4)] - za[jk]*zqsice[jk]) / fmax(zepsec, (double) 1.0 -
         za[jk]);

      //---------------------------------------------
      // humidity in moistest ZCOVPCLR part of domain
      //---------------------------------------------
      zqe = fmax((double) 0.0, fmin(zqe, zqsice[jk]));
      llo1 =
        zcovpclr > zepsec && zqx[jk + klev*(3)] > zepsec && zqe < zzrh*zqsice[jk];

      if (llo1) {

        // Calculate local precipitation (kg/kg)
        zpreclr = zqx[jk + klev*(3)] / zcovptot;
        zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk] - rtt))/(ztp1[jk] - r4ies))))*rv / rd;

        // Particle size distribution
        // ZTCG increases Ni with colder temperatures - essentially a
        // Fletcher or Meyers scheme?
        ztcg = (double) 1.0;            //v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
        // ZFACX1I modification is based on Andrew Barrett's results
        zfacx1s = (double) 1.0;            //v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB

        zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk] +
          pap[jl + klon*(jk + klev*(ibl))]*(*yrecldp).rcl_apb3*(pow(ztp1[jk],
           3));
        zcorrfac = pow((1.0 / zrho), 0.5);
        zcorrfac2 = (pow((ztp1[jk] / 273.0), 1.5))*(393.0 / (ztp1[jk] + 120.0));

        zpr02 = zrho*zpreclr*(*yrecldp).rcl_const1s / (ztcg*zfacx1s);

        zterm1 = (zqsice[jk] - zqe)*(pow(ztp1[jk], 2))*zvpice*zcorrfac2*ztcg*(*yrecldp)
          .rcl_const2s*zfacx1s / (zrho*zaplusb*zqsice[jk]);
        zterm2 = 0.65*(*yrecldp).rcl_const6s*(pow(zpr02, (*yrecldp).rcl_const4s)) +
          (*yrecldp).rcl_const3s*(pow(zcorrfac, 0.5))*(pow(zrho, 0.5))*(pow(zpr02,
          (*yrecldp).rcl_const5s)) / (pow(zcorrfac2, 0.5));

        zdpevap = fmax(zcovpclr*zterm1*zterm2*ptsphy, (double) 0.0);

        //--------------------------------------------------------------------
        // Limit evaporation to snow amount
        //--------------------------------------------------------------------
        zevap = fmin(zdpevap, zevaplimice);
        zevap = fmin(zevap, zqx[jk + klev*(3)]);


        zsolqa[4 + 5*(3)] = zsolqa[4 + 5*(3)] + zevap;
        zsolqa[3 + 5*(4)] = zsolqa[3 + 5*(4)] - zevap;

        //-------------------------------------------------------------
        // Reduce the total precip coverage proportional to evaporation
        // to mimic the previous scheme which had a diagnostic
        // 2-flux treatment, abandoned due to the new prognostic precip
        //-------------------------------------------------------------
        zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot -
          za[jk])*zevap / zqx[jk + klev*(3)]));

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
    zanew = (za[jk] + zsolac) / ((double) 1.0 + zsolab);
    zanew = fmin(zanew, (double) 1.0);
    if (zanew < (*yrecldp).ramin) {
      zanew = (double) 0.0;
    }
    zda = zanew - zaorig[jk];
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
      zsinksum[jm] = (double) 0.0;
    }

    //----------------------------
    // collect sink terms and mark
    //----------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      for (jn = 0; jn <= 5 + -1; jn += 1) {
        zsinksum[jm] = zsinksum[jm] - zsolqa[jm + 5*jn];            // +ve total is bad
      }
    }

    //---------------------------------------
    // calculate overshoot and scaling factor
    //---------------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zmax = fmax(zqx[jk + klev*jm], zepsec);
      zrat = fmax(zsinksum[jm], zmax);
      zratio[jm] = zmax / zrat;
    }

    //--------------------------------------------
    // scale the sink terms, in the correct order,
    // recalculating the scale factor each time
    //--------------------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zsinksum[jm] = (double) 0.0;
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
      zmm = fmax(zqx[jk + klev*jm], zepsec);
      zrr = fmax(zsinksum[jm], zmm);
      zratio[jm] = zmm / zrr;
      //------
      // scale
      //------
      zzratio = zratio[jm];
      //DIR$ IVDEP
      //DIR$ PREFERVECTOR
      for (jn = 0; jn <= 5 + -1; jn += 1) {
        if (zsolqa[jm + 5*jn] < (double) 0.0) {
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
          zqlhs[jn + 5*jm] = (double) 1.0 + zfallsink[jm];
          for (jo = 0; jo <= 5 + -1; jo += 1) {
            zqlhs[jn + 5*jm] = zqlhs[jn + 5*jm] + zsolqb[jo + 5*jn];
          }
          //------------------------------------------
          // non-diagonals: microphysical source terms
          //------------------------------------------
        } else {
          zqlhs[jn + 5*jm] = -zsolqb[jn + 5*jm];              // here is the delta T - missing from doc.
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
      zexplicit = (double) 0.0;
      for (jn = 0; jn <= 5 + -1; jn += 1) {
        zexplicit = zexplicit + zsolqa[jm + 5*jn];            // sum over middle index
      }
      zqxn[jm] = zqx[jk + klev*jm] + zexplicit;
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
        zqxn[jn] = (double) 0.0;
      }
    }

    //--------------------------------
    // variables needed for next level
    //--------------------------------
    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zqxnm1[jm] = zqxn[jm];
      zqxn2d[jk + klev*jm] = zqxn[jm];
    }

    //------------------------------------------------------------------------
    // 5.3 Precipitation/sedimentation fluxes to next level
    //     diagnostic precipitation fluxes
    //     It is this scaled flux that must be used for source to next layer
    //------------------------------------------------------------------------

    for (jm = 0; jm <= 5 + -1; jm += 1) {
      zpfplsx[1 + jk + (klev + 1)*jm] = zfallsink[jm]*zqxn[jm]*zrdtgdp;
    }

    // Ensure precipitation fraction is zero if no precipitation
    zqpretot =
      zpfplsx[1 + jk + (klev + 1)*(3)] + zpfplsx[1 + jk + (klev + 1)*(2)];
    if (zqpretot < zepsec) {
      zcovptot = (double) 0.0;
    }

    //######################################################################
    //              6  *** UPDATE TENDANCIES ***
    //######################################################################

    //--------------------------------
    // 6.1 Temperature and CLV budgets
    //--------------------------------

    for (jm = 0; jm <= 5 - 1 + -1; jm += 1) {

      // calculate fluxes in and out of box for conservation of TL
      zfluxq[jm] = zpsupsatsrce[jm] + zconvsrce[jm] + zfallsrce[jm] - (zfallsink[jm] +
        zconvsink[jm])*zqxn[jm];

      if (iphase[jm] == 1) {
        tendency_loc_t[jl + klon*(jk + klev*(ibl))] = tendency_loc_t[jl
          + klon*(jk + klev*(ibl))] + ralvdcp*(zqxn[jm] - zqx[jk + klev*jm] -
          zfluxq[jm])*zqtmst;
      }

      if (iphase[jm] == 2) {
        tendency_loc_t[jl + klon*(jk + klev*(ibl))] = tendency_loc_t[jl
          + klon*(jk + klev*(ibl))] + ralsdcp*(zqxn[jm] - zqx[jk + klev*jm] -
          zfluxq[jm])*zqtmst;
      }

      //----------------------------------------------------------------------
      // New prognostic tendencies - ice,liquid rain,snow
      // Note: CLV arrays use PCLV in calculation of tendency while humidity
      //       uses ZQX. This is due to clipping at start of cloudsc which
      //       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
      //----------------------------------------------------------------------
      tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*(ibl)))] =
        tendency_loc_cld[jl + klon*(jk + klev*(jm + 5*(ibl)))] + (zqxn[jm] -
        zqx0[jk + klev*jm])*zqtmst;

    }

    //----------------------
    // 6.2 Humidity budget
    //----------------------
    tendency_loc_q[jl + klon*(jk + klev*(ibl))] = tendency_loc_q[jl +
      klon*(jk + klev*(ibl))] + (zqxn[4] - zqx[jk + klev*(4)])*zqtmst;

    //-------------------
    // 6.3 cloud cover
    //-----------------------
    tendency_loc_a[jl + klon*(jk + klev*(ibl))] =
      tendency_loc_a[jl + klon*(jk + klev*(ibl))] + zda*zqtmst;

    //--------------------------------------------------
    // Copy precipitation fraction into output variable
    //-------------------------------------------------
    pcovptot[jl + klon*(jk + klev*(ibl))] = zcovptot;

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
  for (jk = 0; jk <= klev + 1 + -1; jk += 1) {
    pfplsl[jl + klon*(jk + (klev + 1)*(ibl))] =
      zpfplsx[jk + (klev + 1)*(2)] + zpfplsx[jk + (klev + 1)*(0)];
    pfplsn[jl + klon*(jk + (klev + 1)*(ibl))] =
      zpfplsx[jk + (klev + 1)*(3)] + zpfplsx[jk + (klev + 1)*(1)];
  }

  //--------
  // Fluxes:
  //--------
  pfsqlf[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfsqif[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfsqrf[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfsqsf[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfcqlng[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfcqnng[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfcqrng[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;      //rain
  pfcqsng[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;      //snow
  // fluxes due to turbulence
  pfsqltur[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;
  pfsqitur[jl + klon*(0 + (klev + 1)*(ibl))] = (double) 0.0;

  for (jk = 0; jk <= klev + -1; jk += 1) {

    zgdph_r = -zrg_r*(paph[jl + klon*(1 + jk + (klev + 1)*(ibl))] - paph[
       jl + klon*(jk + (klev + 1)*(ibl))])*zqtmst;
    pfsqlf[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqlf[jl + klon*(jk + (klev + 1)*(ibl))];
    pfsqif[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqif[jl + klon*(jk + (klev + 1)*(ibl))];
    pfsqrf[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqlf[jl + klon*(jk + (klev + 1)*(ibl))];
    pfsqsf[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqif[jl + klon*(jk + (klev + 1)*(ibl))];
    pfcqlng[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfcqlng[jl + klon*(jk + (klev + 1)*(ibl))];
    pfcqnng[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfcqnng[jl + klon*(jk + (klev + 1)*(ibl))];
    pfcqrng[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfcqlng[jl + klon*(jk + (klev + 1)*(ibl))];
    pfcqsng[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfcqnng[jl + klon*(jk + (klev + 1)*(ibl))];
    pfsqltur[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqltur[jl + klon*(jk + (klev + 1)*(ibl))];
    pfsqitur[jl + klon*(1 + jk + (klev + 1)*(ibl))] =
      pfsqitur[jl + klon*(jk + (klev + 1)*(ibl))];

    zalfaw = zfoealfa[jk];

    // Liquid , LS scheme minus detrainment
    pfsqlf[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqlf[jl + klon*(1
      + jk + (klev + 1)*(ibl))] + (zqxn2d[jk + klev*(0)] - zqx0[jk + klev*(-1
       + 1)] + pvfl[jl + klon*(jk + klev*(ibl))]*ptsphy - zalfaw*plude[
      jl + klon*(jk + klev*(ibl))])*zgdph_r;
    // liquid, negative numbers
    pfcqlng[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfcqlng[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + zlneg[jk + klev*(0)]*zgdph_r;

    // liquid, vertical diffusion
    pfsqltur[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqltur[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + pvfl[jl + klon*(jk + klev*(ibl
      ))]*ptsphy*zgdph_r;

    // Rain, LS scheme
    pfsqrf[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqrf[jl + klon*(1
      + jk + (klev + 1)*(ibl))] + (zqxn2d[jk + klev*(2)] - zqx0[jk + klev*(-1
       + 3)])*zgdph_r;
    // rain, negative numbers
    pfcqrng[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfcqrng[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + zlneg[jk + klev*(2)]*zgdph_r;

    // Ice , LS scheme minus detrainment
    pfsqif[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqif[jl + klon*(1
      + jk + (klev + 1)*(ibl))] + (zqxn2d[jk + klev*(1)] - zqx0[jk + klev*(-1
       + 2)] + pvfi[jl + klon*(jk + klev*(ibl))]*ptsphy - ((double) 1.0 -
      zalfaw)*plude[jl + klon*(jk + klev*(ibl))])*zgdph_r;
    // ice, negative numbers
    pfcqnng[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfcqnng[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + zlneg[jk + klev*(1)]*zgdph_r;

    // ice, vertical diffusion
    pfsqitur[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqitur[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + pvfi[jl + klon*(jk + klev*(ibl
      ))]*ptsphy*zgdph_r;

    // snow, LS scheme
    pfsqsf[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfsqsf[jl + klon*(1
      + jk + (klev + 1)*(ibl))] + (zqxn2d[jk + klev*(3)] - zqx0[jk + klev*(-1
       + 4)])*zgdph_r;
    // snow, negative numbers
    pfcqsng[jl + klon*(1 + jk + (klev + 1)*(ibl))] = pfcqsng[jl +
      klon*(1 + jk + (klev + 1)*(ibl))] + zlneg[jk + klev*(3)]*zgdph_r;
  }

  //-----------------------------------
  // enthalpy flux due to precipitation
  //-----------------------------------
  for (jk = 0; jk <= klev + 1 + -1; jk += 1) {
    pfhpsl[jl + klon*(jk + (klev + 1)*(ibl))] =
      -rlvtt*pfplsl[jl + klon*(jk + (klev + 1)*(ibl))];
    pfhpsn[jl + klon*(jk + (klev + 1)*(ibl))] =
      -rlstt*pfplsn[jl + klon*(jk + (klev + 1)*(ibl))];
  }
}

