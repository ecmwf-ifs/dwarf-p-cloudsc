#include "cloudsc_c.h"
//#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
//#include "parkind1_c.h"
//#include "yomphyder_c.h"
//#include "yomcst_c.h"
//#include "yoethf_c.h"
//#include "yoecldp_c.h"
//#include "fcttre_mod_c.h"
//#include "fccld_mod_c.h"

void cloudsc_c(int kidia, int kfdia, int klon/*, int klev*/, double ptsphy, double * restrict v_pt, 
  double * restrict v_pq, double * restrict v_tendency_tmp_t, 
  double * restrict v_tendency_tmp_q, double * restrict v_tendency_tmp_a, 
  double * restrict v_tendency_tmp_cld, double * restrict v_tendency_loc_t, 
  double * restrict v_tendency_loc_q, double * restrict v_tendency_loc_a, 
  double * restrict v_tendency_loc_cld, double * restrict v_pvfa, 
  double * restrict v_pvfl, double * restrict v_pvfi, double * restrict v_pdyna, 
  double * restrict v_pdynl, double * restrict v_pdyni, double * restrict v_phrsw, 
  double * restrict v_phrlw, double * restrict v_pvervel, double * restrict v_pap, 
  double * restrict v_paph, double * restrict v_plsm, double * restrict v_ldcum, 
  int * restrict v_ktype, double * restrict v_plu, double * restrict v_plude, 
  double * restrict v_psnde, double * restrict v_pmfu, double * restrict v_pmfd, 
  double * restrict v_pa, double * restrict v_pclv, double * restrict v_psupsat, 
  double * restrict v_plcrit_aer, double * restrict v_picrit_aer, 
  double * restrict v_pre_ice, double * restrict v_pccn, double * restrict v_pnice, 
  double * restrict v_pcovptot, double * restrict v_prainfrac_toprfz, 
  double * restrict v_pfsqlf, double * restrict v_pfsqif, double * restrict v_pfcqnng, 
  double * restrict v_pfcqlng, double * restrict v_pfsqrf, double * restrict v_pfsqsf, 
  double * restrict v_pfcqrng, double * restrict v_pfcqsng, 
  double * restrict v_pfsqltur, double * restrict v_pfsqitur, 
  double * restrict v_pfplsl, double * restrict v_pfplsn, double * restrict v_pfhpsl, 
  double * restrict v_pfhpsn, struct TECLDP *yrecldp, int ngpblks,
  double rg, double rd, double rcpd, double retv, double rlvtt, double rlstt, double rlmlt, double rtt,
  double rv, double r2es, double r3les, double r3ies, double r4les, double r4ies, double r5les,
  double r5ies, double r5alvcp, double r5alscp, double ralvdcp, double ralsdcp, double ralfdcp,
  double rtwat, double rtice, double rticecu, double rtwat_rtice_r, double rtwat_rticecu_r,
  double rkoop1, double rkoop2) {

  /*
  double rg, double rd, double rcpd, double retv, double rlvtt, double rlstt, double rlmlt, double rtt,
  double rv, double r2es, double r3les, double r3ies, double r4les, double r4ies, double r5les,
  double r5ies, double r5alvcp, double r5alscp, double ralvdcp, double ralsdcp, double ralfdcp,
  double rtwat, double rtice, double rticecu, double rtwat_rtice_r, double rtwat_rticecu_r,
  double rkoop1, double rkoop2
  */

  
  //-------------------------------------------------------------------------------
  //                 Declare input/output arguments
  //-------------------------------------------------------------------------------
  
  // PLCRIT_AER : critical liquid mmr for rain autoconversion process
  // PICRIT_AER : critical liquid mmr for snow autoconversion process
  // PRE_LIQ : liq Re
  // PRE_ICE : ice Re
  // PCCN    : liquid cloud condensation nuclei
  // PNICE   : ice number concentration (cf. CCN)
  
  int klev = 137;  // Number of levels
  
  
  
  // Supersat clipped at previous time level in SLTEND
  // Flux diagnostics for DDH budget
  
  
  //-------------------------------------------------------------------------------
  //                       Declare local variables
  //-------------------------------------------------------------------------------
  
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
  int llindex3[5][5];  // index variable
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
  
  double zsolqa[5][5];  // explicit sources and sinks
  double zsolqb[5][5];  // implicit sources and sinks
  // e.g. microphysical pathways between ice variables.
  double zqlhs[5][5];  // n x n matrix storing the LHS of implicit solver
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
  double zqx[5][klev];
  double zqx0[5][klev];
  double zpfplsx[5][klev + 1];
  double zlneg[5][klev];
  double zqxn2d[5][klev];
  /* Array casts for pointer arguments */
  double (*pt)[klev][klon] = (double (*)[klev][klon]) v_pt;
  double (*pq)[klev][klon] = (double (*)[klev][klon]) v_pq;
  double (*tendency_tmp_t)[klev][klon] = (double (*)[klev][klon]) v_tendency_tmp_t;
  double (*tendency_tmp_q)[klev][klon] = (double (*)[klev][klon]) v_tendency_tmp_q;
  double (*tendency_tmp_a)[klev][klon] = (double (*)[klev][klon]) v_tendency_tmp_a;
  double (*tendency_tmp_cld)[5][klev][klon] = (double (*)[5][klev][klon]) v_tendency_tmp_cld;
  double (*tendency_loc_t)[klev][klon] = (double (*)[klev][klon]) v_tendency_loc_t;
  double (*tendency_loc_q)[klev][klon] = (double (*)[klev][klon]) v_tendency_loc_q;
  double (*tendency_loc_a)[klev][klon] = (double (*)[klev][klon]) v_tendency_loc_a;
  double (*tendency_loc_cld)[5][klev][klon] = (double (*)[5][klev][klon]) v_tendency_loc_cld;
  double (*pvfa)[klev][klon] = (double (*)[klev][klon]) v_pvfa;
  double (*pvfl)[klev][klon] = (double (*)[klev][klon]) v_pvfl;
  double (*pvfi)[klev][klon] = (double (*)[klev][klon]) v_pvfi;
  double (*pdyna)[klev][klon] = (double (*)[klev][klon]) v_pdyna;
  double (*pdynl)[klev][klon] = (double (*)[klev][klon]) v_pdynl;
  double (*pdyni)[klev][klon] = (double (*)[klev][klon]) v_pdyni;
  double (*phrsw)[klev][klon] = (double (*)[klev][klon]) v_phrsw;
  double (*phrlw)[klev][klon] = (double (*)[klev][klon]) v_phrlw;
  double (*pvervel)[klev][klon] = (double (*)[klev][klon]) v_pvervel;
  double (*pap)[klev][klon] = (double (*)[klev][klon]) v_pap;
  double (*paph)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_paph;
  double (*plsm)[klon] = (double (*)[klon]) v_plsm;
  int (*ldcum)[klon] = (int (*)[klon]) v_ldcum;
  int (*ktype)[klon] = (int (*)[klon]) v_ktype;
  double (*plu)[klev][klon] = (double (*)[klev][klon]) v_plu;
  double (*plude)[klev][klon] = (double (*)[klev][klon]) v_plude;
  double (*psnde)[klev][klon] = (double (*)[klev][klon]) v_psnde;
  double (*pmfu)[klev][klon] = (double (*)[klev][klon]) v_pmfu;
  double (*pmfd)[klev][klon] = (double (*)[klev][klon]) v_pmfd;
  double (*pa)[klev][klon] = (double (*)[klev][klon]) v_pa;
  double (*pclv)[5][klev][klon] = (double (*)[5][klev][klon]) v_pclv;
  double (*psupsat)[klev][klon] = (double (*)[klev][klon]) v_psupsat;
  double (*plcrit_aer)[klev][klon] = (double (*)[klev][klon]) v_plcrit_aer;
  double (*picrit_aer)[klev][klon] = (double (*)[klev][klon]) v_picrit_aer;
  double (*pre_ice)[klev][klon] = (double (*)[klev][klon]) v_pre_ice;
  double (*pccn)[klev][klon] = (double (*)[klev][klon]) v_pccn;
  double (*pnice)[klev][klon] = (double (*)[klev][klon]) v_pnice;
  double (*pcovptot)[klev][klon] = (double (*)[klev][klon]) v_pcovptot;
  double (*prainfrac_toprfz)[klon] = (double (*)[klon]) v_prainfrac_toprfz;
  double (*pfsqlf)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqlf;
  double (*pfsqif)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqif;
  double (*pfcqnng)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfcqnng;
  double (*pfcqlng)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfcqlng;
  double (*pfsqrf)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqrf;
  double (*pfsqsf)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqsf;
  double (*pfcqrng)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfcqrng;
  double (*pfcqsng)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfcqsng;
  double (*pfsqltur)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqltur;
  double (*pfsqitur)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfsqitur;
  double (*pfplsl)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfplsl;
  double (*pfplsn)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfplsn;
  double (*pfhpsl)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfhpsl;
  double (*pfhpsn)[klev + 1][klon] = (double (*)[klev + 1][klon]) v_pfhpsn;
  /*
  double rg;
  rg = yomcst__get__rg();
  double rd;
  rd = yomcst__get__rd();
  double rcpd;
  rcpd = yomcst__get__rcpd();
  double retv;
  retv = yomcst__get__retv();
  double rlvtt;
  rlvtt = yomcst__get__rlvtt();
  double rlstt;
  rlstt = yomcst__get__rlstt();
  double rlmlt;
  rlmlt = yomcst__get__rlmlt();
  double rtt;
  rtt = yomcst__get__rtt();
  double rv;
  rv = yomcst__get__rv();
  double r2es;
  r2es = yoethf__get__r2es();
  double r3les;
  r3les = yoethf__get__r3les();
  double r3ies;
  r3ies = yoethf__get__r3ies();
  double r4les;
  r4les = yoethf__get__r4les();
  double r4ies;
  r4ies = yoethf__get__r4ies();
  double r5les;
  r5les = yoethf__get__r5les();
  double r5ies;
  r5ies = yoethf__get__r5ies();
  double r5alvcp;
  r5alvcp = yoethf__get__r5alvcp();
  double r5alscp;
  r5alscp = yoethf__get__r5alscp();
  double ralvdcp;
  ralvdcp = yoethf__get__ralvdcp();
  double ralsdcp;
  ralsdcp = yoethf__get__ralsdcp();
  double ralfdcp;
  ralfdcp = yoethf__get__ralfdcp();
  double rtwat;
  rtwat = yoethf__get__rtwat();
  double rtice;
  rtice = yoethf__get__rtice();
  double rticecu;
  rticecu = yoethf__get__rticecu();
  double rtwat_rtice_r;
  rtwat_rtice_r = yoethf__get__rtwat_rtice_r();
  double rtwat_rticecu_r;
  rtwat_rticecu_r = yoethf__get__rtwat_rticecu_r();
  double rkoop1;
  rkoop1 = yoethf__get__rkoop1();
  double rkoop2;
  rkoop2 = yoethf__get__rkoop2();
  */
  jl = 0; //threadIdx.x; //THREADIDX%X;
  ibl = 0; //blockIdx.z; //BLOCKIDX%Z;
  
  if (ibl <= ngpblks && jl <= klon) {
    
    
    
    
    
    
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
    iphase[5 - 1] = 0;
    iphase[1 - 1] = 1;
    iphase[3 - 1] = 1;
    iphase[2 - 1] = 2;
    iphase[4 - 1] = 2;
    
    // ---------------------------------------------------
    // Set up melting/freezing index,
    // if an ice category melts/freezes, where does it go?
    // ---------------------------------------------------
    imelt[5 - 1] = -99;
    imelt[1 - 1] = 2;
    imelt[3 - 1] = 4;
    imelt[2 - 1] = 3;
    imelt[4 - 1] = 3;
    
    // -----------------------------------------------
    // INITIALIZATION OF OUTPUT TENDENCIES
    // -----------------------------------------------
    for (jk = 1; jk <= klev; jk += 1) {
      tendency_loc_t[ibl - 1][jk - 1][jl - 1] = (double) 0.0;
      tendency_loc_q[ibl - 1][jk - 1][jl - 1] = (double) 0.0;
      tendency_loc_a[ibl - 1][jk - 1][jl - 1] = (double) 0.0;
    }
    for (jm = 1; jm <= 5 - 1; jm += 1) {
      for (jk = 1; jk <= klev; jk += 1) {
        tendency_loc_cld[ibl - 1][jm - 1][jk - 1][jl - 1] = (double) 0.0;
      }
    }
    
    //-- These were uninitialized : meaningful only when we compare error differences
    for (jk = 1; jk <= klev; jk += 1) {
      pcovptot[ibl - 1][jk - 1][jl - 1] = (double) 0.0;
      tendency_loc_cld[ibl - 1][5 - 1][jk - 1][jl - 1] = (double) 0.0;
    }
    
    // -------------------------
    // set up fall speeds in m/s
    // -------------------------
    zvqx[5 - 1] = (double) 0.0;
    zvqx[1 - 1] = (double) 0.0;
    zvqx[2 - 1] = (*yrecldp).rvice;
    zvqx[3 - 1] = (*yrecldp).rvrain;
    zvqx[4 - 1] = (*yrecldp).rvsnow;
    for (i_llfall_0 = 1; i_llfall_0 <= 5; i_llfall_0 += 1) {
      llfall[i_llfall_0 - 1] = false;
    }
    for (jm = 1; jm <= 5; jm += 1) {
      if (zvqx[jm - 1] > (double) 0.0) {
        llfall[jm - 1] = true;
      }
      // falling species
    }
    // Set LLFALL to false for ice (but ice still sediments!)
    // Need to rationalise this at some point
    llfall[2 - 1] = false;
    
    
    //######################################################################
    //             1.  *** INITIAL VALUES FOR VARIABLES ***
    //######################################################################
    
    
    // ----------------------
    // non CLV initialization
    // ----------------------
    for (jk = 1; jk <= klev; jk += 1) {
      ztp1[jk - 1] = 
        pt[ibl - 1][jk - 1][jl - 1] + ptsphy*tendency_tmp_t[ibl - 1][jk - 1][jl - 1];
      zqx[5 - 1][jk - 1] = 
        pq[ibl - 1][jk - 1][jl - 1] + ptsphy*tendency_tmp_q[ibl - 1][jk - 1][jl - 1];
      zqx0[5 - 1][jk - 1] = 
        pq[ibl - 1][jk - 1][jl - 1] + ptsphy*tendency_tmp_q[ibl - 1][jk - 1][jl - 1];
      za[jk - 1] = 
        pa[ibl - 1][jk - 1][jl - 1] + ptsphy*tendency_tmp_a[ibl - 1][jk - 1][jl - 1];
      zaorig[jk - 1] = 
        pa[ibl - 1][jk - 1][jl - 1] + ptsphy*tendency_tmp_a[ibl - 1][jk - 1][jl - 1];
    }
    
    // -------------------------------------
    // initialization for CLV family
    // -------------------------------------
    for (jm = 1; jm <= 5 - 1; jm += 1) {
      for (jk = 1; jk <= klev; jk += 1) {
        zqx[jm - 1][jk - 1] = pclv[ibl - 1][jm - 1][jk - 1][jl - 1] + 
          ptsphy*tendency_tmp_cld[ibl - 1][jm - 1][jk - 1][jl - 1];
        zqx0[jm - 1][jk - 1] = pclv[ibl - 1][jm - 1][jk - 1][jl - 1] + 
          ptsphy*tendency_tmp_cld[ibl - 1][jm - 1][jk - 1][jl - 1];
      }
    }
    
    //-------------
    // zero arrays
    //-------------
    for (jm = 1; jm <= 5; jm += 1) {
      for (jk = 1; jk <= klev + 1; jk += 1) {
        zpfplsx[jm - 1][jk - 1] = (double) 0.0;          // precip fluxes
      }
    }
    
    for (jm = 1; jm <= 5; jm += 1) {
      for (jk = 1; jk <= klev; jk += 1) {
        zqxn2d[jm - 1][jk - 1] = (double) 0.0;          // end of timestep values in 2D
        zlneg[jm - 1][jk - 1] = (double) 0.0;          // negative input check
      }
    }
    
    prainfrac_toprfz[ibl - 1][jl - 1] = (double) 0.0;      // rain fraction at top of refreezing layer
    llrainliq = true;      // Assume all raindrops are liquid initially
    
    // ----------------------------------------------------
    // Tidy up very small cloud cover or total cloud water
    // ----------------------------------------------------
    for (jk = 1; jk <= klev; jk += 1) {
      if (zqx[1 - 1][jk - 1] + zqx[2 - 1][jk - 1] < (*yrecldp).rlmin || za[jk - 1] < 
        (*yrecldp).ramin) {
        
        // Evaporate small cloud liquid water amounts
        zlneg[1 - 1][jk - 1] = zlneg[1 - 1][jk - 1] + zqx[1 - 1][jk - 1];
        zqadj = zqx[1 - 1][jk - 1]*zqtmst;
        tendency_loc_q[ibl - 1][jk - 1][jl - 1] = 
          tendency_loc_q[ibl - 1][jk - 1][jl - 1] + zqadj;
        tendency_loc_t[ibl - 1][jk - 1][jl - 1] = 
          tendency_loc_t[ibl - 1][jk - 1][jl - 1] - ralvdcp*zqadj;
        zqx[5 - 1][jk - 1] = zqx[5 - 1][jk - 1] + zqx[1 - 1][jk - 1];
        zqx[1 - 1][jk - 1] = (double) 0.0;
        
        // Evaporate small cloud ice water amounts
        zlneg[2 - 1][jk - 1] = zlneg[2 - 1][jk - 1] + zqx[2 - 1][jk - 1];
        zqadj = zqx[2 - 1][jk - 1]*zqtmst;
        tendency_loc_q[ibl - 1][jk - 1][jl - 1] = 
          tendency_loc_q[ibl - 1][jk - 1][jl - 1] + zqadj;
        tendency_loc_t[ibl - 1][jk - 1][jl - 1] = 
          tendency_loc_t[ibl - 1][jk - 1][jl - 1] - ralsdcp*zqadj;
        zqx[5 - 1][jk - 1] = zqx[5 - 1][jk - 1] + zqx[2 - 1][jk - 1];
        zqx[2 - 1][jk - 1] = (double) 0.0;
        
        // Set cloud cover to zero
        za[jk - 1] = (double) 0.0;
        
      }
    }
    
    // ---------------------------------
    // Tidy up small CLV variables
    // ---------------------------------
    //DIR$ IVDEP
    for (jm = 1; jm <= 5 - 1; jm += 1) {
      //DIR$ IVDEP
      for (jk = 1; jk <= klev; jk += 1) {
        //DIR$ IVDEP
        if (zqx[jm - 1][jk - 1] < (*yrecldp).rlmin) {
          zlneg[jm - 1][jk - 1] = zlneg[jm - 1][jk - 1] + zqx[jm - 1][jk - 1];
          zqadj = zqx[jm - 1][jk - 1]*zqtmst;
          tendency_loc_q[ibl - 1][jk - 1][jl - 1] = 
            tendency_loc_q[ibl - 1][jk - 1][jl - 1] + zqadj;
          if (iphase[jm - 1] == 1) {
            tendency_loc_t[ibl - 1][jk - 1][jl - 1] = 
              tendency_loc_t[ibl - 1][jk - 1][jl - 1] - ralvdcp*zqadj;
          }
          if (iphase[jm - 1] == 2) {
            tendency_loc_t[ibl - 1][jk - 1][jl - 1] = 
              tendency_loc_t[ibl - 1][jk - 1][jl - 1] - ralsdcp*zqadj;
          }
          zqx[5 - 1][jk - 1] = zqx[5 - 1][jk - 1] + zqx[jm - 1][jk - 1];
          zqx[jm - 1][jk - 1] = (double) 0.0;
        }
      }
    }
    
    
    // ------------------------------
    // Define saturation values
    // ------------------------------
    for (jk = 1; jk <= klev; jk += 1) {
      //----------------------------------------
      // old *diagnostic* mixed phase saturation
      //----------------------------------------
      zfoealfa[jk - 1] =  (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))); //foealfa;
      zfoeewmt[jk - 1] = fmin((double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies))))/pap[ibl - 1][jk-1][jl-1], 0.5); //fmin(foeewm / pap[ibl - 1][jk - 1][jl - 1], (double) 0.5);
      zqsmix[jk - 1] = zfoeewmt[jk - 1];
      zqsmix[jk - 1] = zqsmix[jk - 1] / ((double) 1.0 - retv*zqsmix[jk - 1]);
      
      //---------------------------------------------
      // ice saturation T<273K
      // liquid water saturation for T>273K
      //---------------------------------------------
      zalfa = (double)(fmax(0.0, copysign(1.0, ztp1[jk-1] - rtt))); //foedelta;
      zfoeew[jk - 1] = fmin((zalfa*((double)(r2es*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)))) + ((double) 1.0 - zalfa)*((double)(r2es*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies))))) / pap[ibl -
        1][jk - 1][jl - 1], (double) 0.5);
      zfoeew[jk - 1] = fmin((double) 0.5, zfoeew[jk - 1]);
      zqsice[jk - 1] = zfoeew[jk - 1] / ((double) 1.0 - retv*zfoeew[jk - 1]);
      
      //----------------------------------
      // liquid water saturation
      //----------------------------------
      zfoeeliqt[jk - 1] = fmin(((double)(r2es*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)))) / pap[ibl - 1][jk - 1][jl - 1], (double) 0.5);
      zqsliq[jk - 1] = zfoeeliqt[jk - 1];
      zqsliq[jk - 1] = zqsliq[jk - 1] / ((double) 1.0 - retv*zqsliq[jk - 1]);
      
      //   !----------------------------------
      //   ! ice water saturation
      //   !----------------------------------
      //   ZFOEEICET(JL,JK)=MIN(FOEEICE(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      //   ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
      //   ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-RETV*ZQSICE(JL,JK))
      
    }
    
    for (jk = 1; jk <= klev; jk += 1) {
      
      
      //------------------------------------------
      // Ensure cloud fraction is between 0 and 1
      //------------------------------------------
      za[jk - 1] = fmax((double) 0.0, fmin((double) 1.0, za[jk - 1]));
      
      //-------------------------------------------------------------------
      // Calculate liq/ice fractions (no longer a diagnostic relationship)
      //-------------------------------------------------------------------
      zli[jk - 1] = zqx[1 - 1][jk - 1] + zqx[2 - 1][jk - 1];
      if (zli[jk - 1] > (*yrecldp).rlmin) {
        zliqfrac[jk - 1] = zqx[1 - 1][jk - 1] / zli[jk - 1];
        zicefrac[jk - 1] = (double) 1.0 - zliqfrac[jk - 1];
      } else {
        zliqfrac[jk - 1] = (double) 0.0;
        zicefrac[jk - 1] = (double) 0.0;
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
    zpaphd = (double) 1.0 / paph[ibl - 1][klev + 1 - 1][jl - 1];
    for (jk = 1; jk <= klev - 1; jk += 1) {
      zsig = pap[ibl - 1][jk - 1][jl - 1]*zpaphd;
      if (zsig > (double) 0.1 && zsig < (double) 0.4 && ztp1[jk - 1] > ztp1[jk + 1 - 1]
        ) {
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
    
    for (jk = (*yrecldp).ncldtop; jk <= klev; jk += 1) {
      
      //----------------------------------------------------------------------
      // 3.0 INITIALIZE VARIABLES
      //----------------------------------------------------------------------
      
      //---------------------------------
      // First guess microphysics
      //---------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        zqxfg[jm - 1] = zqx[jm - 1][jk - 1];
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
      for (jm = 1; jm <= 5; jm += 1) {
        for (jn = 1; jn <= 5; jn += 1) {
          zsolqb[jm - 1][jn - 1] = (double) 0.0;
          zsolqa[jm - 1][jn - 1] = (double) 0.0;
        }
      }
      
      //----------------------------------
      // reset new microphysics variables
      //----------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        zfallsrce[jm - 1] = (double) 0.0;
        zfallsink[jm - 1] = (double) 0.0;
        zconvsrce[jm - 1] = (double) 0.0;
        zconvsink[jm - 1] = (double) 0.0;
        zpsupsatsrce[jm - 1] = (double) 0.0;
        zratio[jm - 1] = (double) 0.0;
      }
      
      
      //-------------------------
      // derived variables needed
      //-------------------------
      
      zdp = paph[ibl - 1][jk + 1 - 1][jl - 1] - paph[ibl - 1][jk - 1][jl - 1];        // dp
      zgdp = rg / zdp;        // g/dp
      zrho = pap[ibl - 1][jk - 1][jl - 1] / (rd*ztp1[jk - 1]);        // p/RT air density
      
      zdtgdp = ptsphy*zgdp;        // dt g/dp
      zrdtgdp = zdp*((double) 1.0 / (ptsphy*rg));        // 1/(dt g/dp)
      
      if (jk > 1) {
        zdtgdpf = 
          ptsphy*rg / (pap[ibl - 1][jk - 1][jl - 1] - pap[ibl - 1][jk - 1 - 1][jl - 1]);
      }
      
      //------------------------------------
      // Calculate dqs/dT correction factor
      //------------------------------------
      // Reminder: RETV=RV/RD-1
      
      // liquid
      zfacw = r5les / (pow((ztp1[jk - 1] - r4les), 2));
      zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeeliqt[jk - 1]);
      zdqsliqdt = zfacw*zcor*zqsliq[jk - 1];
      zcorqsliq = (double) 1.0 + ralvdcp*zdqsliqdt;
      
      // ice
      zfaci = r5ies / (pow((ztp1[jk - 1] - r4ies), 2));
      zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeew[jk - 1]);
      zdqsicedt = zfaci*zcor*zqsice[jk - 1];
      zcorqsice = (double) 1.0 + ralsdcp*zdqsicedt;
      
      // diagnostic mixed
      zalfaw = zfoealfa[jk - 1];
      zalfawm = zalfaw;
      zfac = zalfaw*zfacw + ((double) 1.0 - zalfaw)*zfaci;
      zcor = (double) 1.0 / ((double) 1.0 - retv*zfoeewmt[jk - 1]);
      zdqsmixdt = zfac*zcor*zqsmix[jk - 1];
      zcorqsmix = (double) 1.0 + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)) *zdqsmixdt;
      
      // evaporation/sublimation limits
      zevaplimmix = fmax((zqsmix[jk - 1] - zqx[5 - 1][jk - 1]) / zcorqsmix, (double) 0.0)
        ;
      zevaplimliq = fmax((zqsliq[jk - 1] - zqx[5 - 1][jk - 1]) / zcorqsliq, (double) 0.0)
        ;
      zevaplimice = fmax((zqsice[jk - 1] - zqx[5 - 1][jk - 1]) / zcorqsice, (double) 0.0)
        ;
      
      //--------------------------------
      // in-cloud consensate amount
      //--------------------------------
      ztmpa = (double) 1.0 / fmax(za[jk - 1], zepsec);
      zliqcld = zqx[1 - 1][jk - 1]*ztmpa;
      zicecld = zqx[2 - 1][jk - 1]*ztmpa;
      zlicld = zliqcld + zicecld;
      
      
      //------------------------------------------------
      // Evaporate very small amounts of liquid and ice
      //------------------------------------------------
      
      if (zqx[1 - 1][jk - 1] < (*yrecldp).rlmin) {
        zsolqa[1 - 1][5 - 1] = zqx[1 - 1][jk - 1];
        zsolqa[5 - 1][1 - 1] = -zqx[1 - 1][jk - 1];
      }
      
      if (zqx[2 - 1][jk - 1] < (*yrecldp).rlmin) {
        zsolqa[2 - 1][5 - 1] = zqx[2 - 1][jk - 1];
        zsolqa[5 - 1][2 - 1] = -zqx[2 - 1][jk - 1];
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
      zfokoop = (double)(fmin(rkoop1 - rkoop2*ztp1[jk-1], (double)(r2es*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)))*1.0/(double)(r2es*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))));//fokoop;
      
      if (ztp1[jk - 1] >= rtt || (*yrecldp).nssopt == 0) {
        zfac = (double) 1.0;
        zfaci = (double) 1.0;
      } else {
        zfac = za[jk - 1] + zfokoop*((double) 1.0 - za[jk - 1]);
        zfaci = ptsphy / (*yrecldp).rkooptau;
      }
      
      //-------------------------------------------------------------------
      // 3.1.2 Calculate supersaturation wrt Koop including dqs/dT
      //       correction factor
      // [#Note: QSICE or QSLIQ]
      //-------------------------------------------------------------------
      
      // Calculate supersaturation to add to cloud
      if (za[jk - 1] > (double) 1.0 - (*yrecldp).ramin) {
        zsupsat = 
          fmax((zqx[5 - 1][jk - 1] - zfac*zqsice[jk - 1]) / zcorqsice, (double) 0.0);
      } else {
        // Calculate environmental humidity supersaturation
        zqp1env = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsice[jk - 1]) / fmax((double) 1.0 - 
          za[jk - 1], zepsilon);
        //& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
        zsupsat = fmax(((double) 1.0 - za[jk - 1])*(zqp1env - zfac*zqsice[jk - 1]) / 
          zcorqsice, (double) 0.0);
      }
      
      //-------------------------------------------------------------------
      // Here the supersaturation is turned into liquid water
      // However, if the temperature is below the threshold for homogeneous
      // freezing then the supersaturation is turned instantly to ice.
      //--------------------------------------------------------------------
      
      if (zsupsat > zepsec) {
        
        if (ztp1[jk - 1] > (*yrecldp).rthomo) {
          // Turn supersaturation into liquid water
          zsolqa[5 - 1][1 - 1] = zsolqa[5 - 1][1 - 1] + zsupsat;
          zsolqa[1 - 1][5 - 1] = zsolqa[1 - 1][5 - 1] - zsupsat;
          // Include liquid in first guess
          zqxfg[1 - 1] = zqxfg[1 - 1] + zsupsat;
        } else {
          // Turn supersaturation into ice water
          zsolqa[5 - 1][2 - 1] = zsolqa[5 - 1][2 - 1] + zsupsat;
          zsolqa[2 - 1][5 - 1] = zsolqa[2 - 1][5 - 1] - zsupsat;
          // Add ice to first guess for deposition term
          zqxfg[2 - 1] = zqxfg[2 - 1] + zsupsat;
        }
        
        // Increase cloud amount using RKOOPTAU timescale
        zsolac = ((double) 1.0 - za[jk - 1])*zfaci;
        
      }
      
      //-------------------------------------------------------
      // 3.1.3 Include supersaturation from previous timestep
      // (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
      //-------------------------------------------------------
      if (psupsat[ibl - 1][jk - 1][jl - 1] > zepsec) {
        if (ztp1[jk - 1] > (*yrecldp).rthomo) {
          // Turn supersaturation into liquid water
          zsolqa[1 - 1][1 - 1] = zsolqa[1 - 1][1 - 1] + psupsat[ibl - 1][jk - 1][jl - 1];
          zpsupsatsrce[1 - 1] = psupsat[ibl - 1][jk - 1][jl - 1];
          // Add liquid to first guess for deposition term
          zqxfg[1 - 1] = zqxfg[1 - 1] + psupsat[ibl - 1][jk - 1][jl - 1];
          // Store cloud budget diagnostics if required
        } else {
          // Turn supersaturation into ice water
          zsolqa[2 - 1][2 - 1] = zsolqa[2 - 1][2 - 1] + psupsat[ibl - 1][jk - 1][jl - 1];
          zpsupsatsrce[2 - 1] = psupsat[ibl - 1][jk - 1][jl - 1];
          // Add ice to first guess for deposition term
          zqxfg[2 - 1] = zqxfg[2 - 1] + psupsat[ibl - 1][jk - 1][jl - 1];
          // Store cloud budget diagnostics if required
        }
        
        // Increase cloud amount using RKOOPTAU timescale
        zsolac = ((double) 1.0 - za[jk - 1])*zfaci;
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
      if (jk < klev && jk >= (*yrecldp).ncldtop) {
        
        
        plude[ibl - 1][jk - 1][jl - 1] = plude[ibl - 1][jk - 1][jl - 1]*zdtgdp;
        
        if (ldcum[ibl - 1][jl - 1] && plude[ibl - 1][jk - 1][jl - 1] > (*yrecldp).rlmin 
          && plu[ibl - 1][jk + 1 - 1][jl - 1] > zepsec) {
          
          zsolac = 
            zsolac + plude[ibl - 1][jk - 1][jl - 1] / plu[ibl - 1][jk + 1 - 1][jl - 1];
          // *diagnostic temperature split*
          zalfaw = zfoealfa[jk - 1];
          zconvsrce[1 - 1] = zalfaw*plude[ibl - 1][jk - 1][jl - 1];
          zconvsrce[2 - 1] = ((double) 1.0 - zalfaw)*plude[ibl - 1][jk - 1][jl - 1];
          zsolqa[1 - 1][1 - 1] = zsolqa[1 - 1][1 - 1] + zconvsrce[1 - 1];
          zsolqa[2 - 1][2 - 1] = zsolqa[2 - 1][2 - 1] + zconvsrce[2 - 1];
          
        } else {
          
          plude[ibl - 1][jk - 1][jl - 1] = (double) 0.0;
          
        }
        // *convective snow detrainment source
        if (ldcum[ibl - 1][jl - 1]) {
          zsolqa[4 - 1][4 - 1] = 
            zsolqa[4 - 1][4 - 1] + psnde[ibl - 1][jk - 1][jl - 1]*zdtgdp;
        }
        
        
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
      if (jk > (*yrecldp).ncldtop) {
        
        zmf = fmax((double) 0.0, (pmfu[ibl - 1][jk - 1][jl - 1] + pmfd[ibl - 1][jk - 
          1][jl - 1])*zdtgdp);
        zacust = zmf*zanewm1;
        
        for (jm = 1; jm <= 5; jm += 1) {
          if (!llfall[jm - 1] && iphase[jm - 1] > 0) {
            zlcust[jm - 1] = zmf*zqxnm1[jm - 1];
            // record total flux for enthalpy budget:
            zconvsrce[jm - 1] = zconvsrce[jm - 1] + zlcust[jm - 1];
          }
        }
        
        // Now have to work out how much liquid evaporates at arrival point
        // since there is no prognostic memory for in-cloud humidity, i.e.
        // we always assume cloud is saturated.
        
        zdtdp = zrdcp*(double) 0.5*(ztp1[jk - 1 - 1] + ztp1[jk - 1]) / paph[ibl - 1][jk -
           1][jl - 1];
        zdtforc = zdtdp*(pap[ibl - 1][jk - 1][jl - 1] - pap[ibl - 1][jk - 1 - 1][jl - 1])
          ;
        //[#Note: Diagnostic mixed phase should be replaced below]
        zdqs = zanewm1*zdtforc*zdqsmixdt;
        
        for (jm = 1; jm <= 5; jm += 1) {
          if (!llfall[jm - 1] && iphase[jm - 1] > 0) {
            zlfinal = fmax((double) 0.0, zlcust[jm - 1] - zdqs);              //lim to zero
            // no supersaturation allowed incloud ---V
            zevap = fmin((zlcust[jm - 1] - zlfinal), zevaplimmix);
            //          ZEVAP=0.0_JPRB
            zlfinal = zlcust[jm - 1] - zevap;
            zlfinalsum = zlfinalsum + zlfinal;              // sum
            
            zsolqa[jm - 1][jm - 1] = zsolqa[jm - 1][jm - 1] + zlcust[jm - 1];              // whole sum
            zsolqa[jm - 1][5 - 1] = zsolqa[jm - 1][5 - 1] + zevap;
            zsolqa[5 - 1][jm - 1] = zsolqa[5 - 1][jm - 1] - zevap;
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
      
      
      if (jk < klev) {
        
        zmfdn = fmax((double) 0.0, (pmfu[ibl - 1][jk + 1 - 1][jl - 1] + pmfd[ibl - 1][jk 
          + 1 - 1][jl - 1])*zdtgdp);
        
        zsolab = zsolab + zmfdn;
        zsolqb[1 - 1][1 - 1] = zsolqb[1 - 1][1 - 1] + zmfdn;
        zsolqb[2 - 1][2 - 1] = zsolqb[2 - 1][2 - 1] + zmfdn;
        
        // Record sink for cloud budget and enthalpy budget diagnostics
        zconvsink[1 - 1] = zmfdn;
        zconvsink[2 - 1] = zmfdn;
        
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
      if (ktype[ibl - 1][jl - 1] > 0 && plude[ibl - 1][jk - 1][jl - 1] > zepsec) {
        zldifdt = (*yrecldp).rcldiff_convi*zldifdt;
      }
      
      // At the moment, works on mixed RH profile and partitioned ice/liq fraction
      // so that it is similar to previous scheme
      // Should apply RHw for liquid cloud and RHi for ice cloud separately
      if (zli[jk - 1] > zepsec) {
        // Calculate environmental humidity
        //      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
        //    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
        //      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
        ze = zldifdt*fmax(zqsmix[jk - 1] - zqx[5 - 1][jk - 1], (double) 0.0);
        zleros = za[jk - 1]*ze;
        zleros = fmin(zleros, zevaplimmix);
        zleros = fmin(zleros, zli[jk - 1]);
        zaeros = zleros / zlicld;          //if linear term
        
        // Erosion is -ve LINEAR in L,A
        zsolac = zsolac - zaeros;          //linear
        
        zsolqa[1 - 1][5 - 1] = zsolqa[1 - 1][5 - 1] + zliqfrac[jk - 1]*zleros;
        zsolqa[5 - 1][1 - 1] = zsolqa[5 - 1][1 - 1] - zliqfrac[jk - 1]*zleros;
        zsolqa[2 - 1][5 - 1] = zsolqa[2 - 1][5 - 1] + zicefrac[jk - 1]*zleros;
        zsolqa[5 - 1][2 - 1] = zsolqa[5 - 1][2 - 1] - zicefrac[jk - 1]*zleros;
        
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
      
      zdtdp = zrdcp*ztp1[jk - 1] / pap[ibl - 1][jk - 1][jl - 1];
      zdpmxdt = zdp*zqtmst;
      zmfdn = (double) 0.0;
      if (jk < klev) {
        zmfdn = pmfu[ibl - 1][jk + 1 - 1][jl - 1] + pmfd[ibl - 1][jk + 1 - 1][jl - 1];
      }
      zwtot = pvervel[ibl - 1][jk - 1][jl - 1] + (double) 0.5*rg*(pmfu[ibl - 1][jk - 
        1][jl - 1] + pmfd[ibl - 1][jk - 1][jl - 1] + zmfdn);
      zwtot = fmin(zdpmxdt, fmax(-zdpmxdt, zwtot));
      zzzdt = phrsw[ibl - 1][jk - 1][jl - 1] + phrlw[ibl - 1][jk - 1][jl - 1];
      zdtdiab = fmin(zdpmxdt*zdtdp, fmax(-zdpmxdt*zdtdp, zzzdt))*ptsphy + ralfdcp*zldefr;
      // Note: ZLDEFR should be set to the difference between the mixed phase functions
      // in the convection and cloud scheme, but this is not calculated, so is zero and
      // the functions must be the same
      zdtforc = zdtdp*zwtot*ptsphy + zdtdiab;
      zqold = zqsmix[jk - 1];
      ztold = ztp1[jk - 1];
      ztp1[jk - 1] = ztp1[jk - 1] + zdtforc;
      ztp1[jk - 1] = fmax(ztp1[jk - 1], (double) 160.0);
      llflag = true;
      
      // Formerly a call to CUADJTQ(..., ICALL=5)
      zqp = (double) 1.0 / pap[ibl - 1][jk - 1][jl - 1];
      zqsat = (double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))) * zqp;//foeewm*zqp;
      zqsat = fmin((double) 0.5, zqsat);
      zcor = (double) 1.0 / ((double) 1.0 - retv*zqsat);
      zqsat = zqsat*zcor;
      // (double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1] - r4ies, 2))
      zcond = (zqsmix[jk - 1] - zqsat) / ((double) 1.0 + zqsat*zcor*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1] - r4ies, 2)))));
      ztp1[jk - 1] = ztp1[jk - 1] + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)) *zcond;
      zqsmix[jk - 1] = zqsmix[jk - 1] - zcond;
      zqsat = (double)(r2es*((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)) + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))) * zqp;//foeewm*zqp;
      zqsat = fmin((double) 0.5, zqsat);
      zcor = (double) 1.0 / ((double) 1.0 - retv*zqsat);
      zqsat = zqsat*zcor;
      zcond1 = (zqsmix[jk - 1] - zqsat) / ((double) 1.0 + zqsat*zcor*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1] - r4ies, 2)))));
      ztp1[jk - 1] = ztp1[jk - 1] + ((double)((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*ralvdcp + (1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*ralsdcp)) *zcond1;
      zqsmix[jk - 1] = zqsmix[jk - 1] - zcond1;
      
      zdqs = zqsmix[jk - 1] - zqold;
      zqsmix[jk - 1] = zqold;
      ztp1[jk - 1] = ztold;
      
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
        zlevap = za[jk - 1]*fmin(zdqs, zlicld);
        zlevap = fmin(zlevap, zevaplimmix);
        zlevap = fmin(zlevap, fmax(zqsmix[jk - 1] - zqx[5 - 1][jk - 1], (double) 0.0));
        
        // For first guess call
        zlevapl = zliqfrac[jk - 1]*zlevap;
        zlevapi = zicefrac[jk - 1]*zlevap;
        
        zsolqa[1 - 1][5 - 1] = zsolqa[1 - 1][5 - 1] + zliqfrac[jk - 1]*zlevap;
        zsolqa[5 - 1][1 - 1] = zsolqa[5 - 1][1 - 1] - zliqfrac[jk - 1]*zlevap;
        
        zsolqa[2 - 1][5 - 1] = zsolqa[2 - 1][5 - 1] + zicefrac[jk - 1]*zlevap;
        zsolqa[5 - 1][2 - 1] = zsolqa[5 - 1][2 - 1] - zicefrac[jk - 1]*zlevap;
        
      }
      
      
      //----------------------------------------------------------------------
      // 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
      //----------------------------------------------------------------------
      // (1) Increase of cloud water in existing clouds
      if (za[jk - 1] > zepsec && zdqs <= -(*yrecldp).rlmin) {
        
        zlcond1 = fmax(-zdqs, (double) 0.0);          //new limiter
        
        //old limiter (significantly improves upper tropospheric humidity rms)
        if (za[jk - 1] > (double) 0.99) {
          zcor = (double) 1.0 / ((double) 1.0 - retv*zqsmix[jk - 1]);
          zcdmax = (zqx[5 - 1][jk - 1] - zqsmix[jk - 1]) / ((double) 1.0 + zcor*zqsmix[jk
             - 1]*((double)(((double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2)))*r5alvcp)*(1.0/pow(ztp1[jk-1] - r4les, 2)) + ((1.0 - (double)(fmin(1.0, pow((fmax(rtice, fmin(rtwat, ztp1[jk-1])) - rtice)*rtwat_rtice_r, 2))))*r5alscp)*(1.0/pow(ztp1[jk-1] - r4ies, 2)))));
        } else {
          zcdmax = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsmix[jk - 1]) / za[jk - 1];
        }
        zlcond1 = fmax(fmin(zlcond1, zcdmax), (double) 0.0);
        // end old limiter
        
        zlcond1 = za[jk - 1]*zlcond1;
        if (zlcond1 < (*yrecldp).rlmin) {
          zlcond1 = (double) 0.0;
        }
        
        //-------------------------------------------------------------------------
        // All increase goes into liquid unless so cold cloud homogeneously freezes
        // Include new liquid formation in first guess value, otherwise liquid
        // remains at cold temperatures until next timestep.
        //-------------------------------------------------------------------------
        if (ztp1[jk - 1] > (*yrecldp).rthomo) {
          zsolqa[5 - 1][1 - 1] = zsolqa[5 - 1][1 - 1] + zlcond1;
          zsolqa[1 - 1][5 - 1] = zsolqa[1 - 1][5 - 1] - zlcond1;
          zqxfg[1 - 1] = zqxfg[1 - 1] + zlcond1;
        } else {
          zsolqa[5 - 1][2 - 1] = zsolqa[5 - 1][2 - 1] + zlcond1;
          zsolqa[2 - 1][5 - 1] = zsolqa[2 - 1][5 - 1] - zlcond1;
          zqxfg[2 - 1] = zqxfg[2 - 1] + zlcond1;
        }
      }
      
      // (2) Generation of new clouds (da/dt>0)
      
      
      if (zdqs <= -(*yrecldp).rlmin && za[jk - 1] < (double) 1.0 - zepsec) {
        
        //---------------------------
        // Critical relative humidity
        //---------------------------
        zrhc = (*yrecldp).ramid;
        zsigk = pap[ibl - 1][jk - 1][jl - 1] / paph[ibl - 1][klev + 1 - 1][jl - 1];
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
          zqe = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsice[jk - 1]) / fmax(zepsec, (double) 
            1.0 - za[jk - 1]);
          zqe = fmax((double) 0.0, zqe);
        } else if ((*yrecldp).nssopt == 1) {
          // Tompkins
          zqe = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsice[jk - 1]) / fmax(zepsec, (double) 
            1.0 - za[jk - 1]);
          zqe = fmax((double) 0.0, zqe);
        } else if ((*yrecldp).nssopt == 2) {
          // Lohmann and Karcher
          zqe = zqx[5 - 1][jk - 1];
        } else if ((*yrecldp).nssopt == 3) {
          // Gierens
          zqe = zqx[5 - 1][jk - 1] + zli[jk - 1];
        }
        
        if (ztp1[jk - 1] >= rtt || (*yrecldp).nssopt == 0) {
          // No ice supersaturation allowed
          zfac = (double) 1.0;
        } else {
          // Ice supersaturation
          zfac = zfokoop;
        }
        
        if (zqe >= zrhc*zqsice[jk - 1]*zfac && zqe < zqsice[jk - 1]*zfac) {
          // note: not **2 on 1-a term if ZQE is used.
          // Added correction term ZFAC to numerator 15/03/2010
          zacond = -((double) 1.0 - za[jk - 1])*zfac*zdqs / fmax((double) 
            2.0*(zfac*zqsice[jk - 1] - zqe), zepsec);
          
          zacond = fmin(zacond, (double) 1.0 - za[jk - 1]);            //PUT THE LIMITER BACK
          
          // Linear term:
          // Added correction term ZFAC 15/03/2010
          zlcond2 = -zfac*zdqs*(double) 0.5*zacond;            //mine linear
          
          // new limiter formulation
          zzdl = (double) 2.0*(zfac*zqsice[jk - 1] - zqe) / fmax(zepsec, (double) 1.0 - 
            za[jk - 1]);
          // Added correction term ZFAC 15/03/2010
          if (zfac*zdqs < -zzdl) {
            // ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
            zlcondlim = (za[jk - 1] - (double) 1.0)*zfac*zdqs - zfac*zqsice[jk - 1] + 
              zqx[5 - 1][jk - 1];
            zlcond2 = fmin(zlcond2, zlcondlim);
          }
          zlcond2 = fmax(zlcond2, (double) 0.0);
          
          if (zlcond2 < (*yrecldp).rlmin || ((double) 1.0 - za[jk - 1]) < zepsec) {
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
          if (ztp1[jk - 1] > (*yrecldp).rthomo) {
            zsolqa[5 - 1][1 - 1] = zsolqa[5 - 1][1 - 1] + zlcond2;
            zsolqa[1 - 1][5 - 1] = zsolqa[1 - 1][5 - 1] - zlcond2;
            zqxfg[1 - 1] = zqxfg[1 - 1] + zlcond2;
          } else {
            // homogeneous freezing
            zsolqa[5 - 1][2 - 1] = zsolqa[5 - 1][2 - 1] + zlcond2;
            zsolqa[2 - 1][5 - 1] = zsolqa[2 - 1][5 - 1] - zlcond2;
            zqxfg[2 - 1] = zqxfg[2 - 1] + zlcond2;
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
        
        if (za[jk - 1 - 1] < (*yrecldp).rcldtopcf && za[jk - 1] >= (*yrecldp).rcldtopcf
          ) {
          zcldtopdist = (double) 0.0;
        } else {
          zcldtopdist = zcldtopdist + zdp / (zrho*rg);
        }
        
        //--------------------------------------------------------------
        // only treat depositional growth if liquid present. due to fact
        // that can not model ice growth from vapour without additional
        // in-cloud water vapour variable
        //--------------------------------------------------------------
        if (ztp1[jk - 1] < rtt && zqxfg[1 - 1] > (*yrecldp).rlmin) {
          // T<273K
          
          zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))*rv) / rd;
          zvpliq = zvpice*zfokoop;
          zicenuclei = (double) 1000.0*exp((double) 12.96*(zvpliq - zvpice) / zvpliq - 
            (double) 0.639);
          
          //------------------------------------------------
          //   2.4e-2 is conductivity of air
          //   8.8 = 700**1/3 = density of ice to the third
          //------------------------------------------------
          zadd = rlstt*(rlstt / (rv*ztp1[jk - 1]) - (double) 1.0) / ((double) 
            2.4E-2*ztp1[jk - 1]);
          zbdd = rv*ztp1[jk - 1]*pap[ibl - 1][jk - 1][jl - 1] / ((double) 2.21*zvpice);
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
          zdepos = fmax(za[jk - 1]*(zinew - zice0), (double) 0.0);
          
          //--------------------------------------------------------------------
          // Limit deposition to liquid water amount
          // If liquid is all frozen, ice would use up reservoir of water
          // vapour in excess of ice saturation mixing ratio - However this
          // can not be represented without a in-cloud humidity variable. Using
          // the grid-mean humidity would imply a large artificial horizontal
          // flux from the clear sky to the cloudy area. We thus rely on the
          // supersaturation check to clean up any remaining supersaturation
          //--------------------------------------------------------------------
          zdepos = fmin(zdepos, zqxfg[1 - 1]);            // limit to liquid water amount
          
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
          zsolqa[1 - 1][2 - 1] = zsolqa[1 - 1][2 - 1] + zdepos;
          zsolqa[2 - 1][1 - 1] = zsolqa[2 - 1][1 - 1] - zdepos;
          zqxfg[2 - 1] = zqxfg[2 - 1] + zdepos;
          zqxfg[1 - 1] = zqxfg[1 - 1] - zdepos;
          
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
        
        if (za[jk - 1 - 1] < (*yrecldp).rcldtopcf && za[jk - 1] >= (*yrecldp).rcldtopcf
          ) {
          zcldtopdist = (double) 0.0;
        } else {
          zcldtopdist = zcldtopdist + zdp / (zrho*rg);
        }
        
        //--------------------------------------------------------------
        // only treat depositional growth if liquid present. due to fact
        // that can not model ice growth from vapour without additional
        // in-cloud water vapour variable
        //--------------------------------------------------------------
        if (ztp1[jk - 1] < rtt && zqxfg[1 - 1] > (*yrecldp).rlmin) {
          // T<273K
          
          zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))*rv) / rd;
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
          
          zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk - 1] 
            + pap[ibl - 1][jk - 1][jl - 1]*(*yrecldp).rcl_apb3*(pow(ztp1[jk - 1], (double
            ) 3.));
          zcorrfac = pow(((double) 1.0 / zrho), (double) 0.5);
          zcorrfac2 = (pow((ztp1[jk - 1] / (double) 273.0), (double) 1.5))*((double) 
            393.0 / (ztp1[jk - 1] + (double) 120.0));
          
          zpr02 = zrho*zice0*(*yrecldp).rcl_const1i / (ztcg*zfacx1i);
          
          zterm1 = (zvpliq - zvpice)*(pow(ztp1[jk - 1], (double) 2.0))
            *zvpice*zcorrfac2*ztcg*(*yrecldp).rcl_const2i*zfacx1i / (zrho*zaplusb*zvpice)
            ;
          zterm2 = (double) 0.65*(*yrecldp).rcl_const6i*(pow(zpr02, (*yrecldp)
            .rcl_const4i)) + (*yrecldp).rcl_const3i*(pow(zcorrfac, (double) 0.5))
            *(pow(zrho, (double) 0.5))*(pow(zpr02, (*yrecldp).rcl_const5i)) / 
            (pow(zcorrfac2, (double) 0.5));
          
          zdepos = fmax(za[jk - 1]*zterm1*zterm2*ptsphy, (double) 0.0);
          
          //--------------------------------------------------------------------
          // Limit deposition to liquid water amount
          // If liquid is all frozen, ice would use up reservoir of water
          // vapour in excess of ice saturation mixing ratio - However this
          // can not be represented without a in-cloud humidity variable. Using
          // the grid-mean humidity would imply a large artificial horizontal
          // flux from the clear sky to the cloudy area. We thus rely on the
          // supersaturation check to clean up any remaining supersaturation
          //--------------------------------------------------------------------
          zdepos = fmin(zdepos, zqxfg[1 - 1]);            // limit to liquid water amount
          
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
          zsolqa[1 - 1][2 - 1] = zsolqa[1 - 1][2 - 1] + zdepos;
          zsolqa[2 - 1][1 - 1] = zsolqa[2 - 1][1 - 1] - zdepos;
          zqxfg[2 - 1] = zqxfg[2 - 1] + zdepos;
          zqxfg[1 - 1] = zqxfg[1 - 1] - zdepos;
        }
        
      }
      // on IDEPICE
      
      //######################################################################
      //              4  *** PRECIPITATION PROCESSES ***
      //######################################################################
      
      //----------------------------------
      // revise in-cloud consensate amount
      //----------------------------------
      ztmpa = (double) 1.0 / fmax(za[jk - 1], zepsec);
      zliqcld = zqxfg[1 - 1]*ztmpa;
      zicecld = zqxfg[2 - 1]*ztmpa;
      zlicld = zliqcld + zicecld;
      
      //----------------------------------------------------------------------
      // 4.2 SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
      //     now that rain, snow, graupel species are prognostic
      //     the precipitation flux can be defined directly level by level
      //     There is no vertical memory required from the flux variable
      //----------------------------------------------------------------------
      
      for (jm = 1; jm <= 5; jm += 1) {
        if (llfall[jm - 1] || jm == 2) {
          //------------------------
          // source from layer above
          //------------------------
          if (jk > (*yrecldp).ncldtop) {
            zfallsrce[jm - 1] = zpfplsx[jm - 1][jk - 1]*zdtgdp;
            zsolqa[jm - 1][jm - 1] = zsolqa[jm - 1][jm - 1] + zfallsrce[jm - 1];
            zqxfg[jm - 1] = zqxfg[jm - 1] + zfallsrce[jm - 1];
            // use first guess precip----------V
            zqpretot = zqpretot + zqxfg[jm - 1];
          }
          //-------------------------------------------------
          // sink to next layer, constant fall speed
          //-------------------------------------------------
          // if aerosol effect then override
          //  note that for T>233K this is the same as above.
          if ((*yrecldp).laericesed && jm == 2) {
            zre_ice = pre_ice[ibl - 1][jk - 1][jl - 1];
            // The exponent value is from
            // Morrison et al. JAS 2005 Appendix
            zvqx[2 - 1] = (double) 0.002*(pow(zre_ice, (double) 1.0));
          }
          zfall = zvqx[jm - 1]*zrho;
          //-------------------------------------------------
          // modified by Heymsfield and Iaquinta JAS 2000
          //-------------------------------------------------
          // ZFALL = ZFALL*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
          //            &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))
          
          zfallsink[jm - 1] = zdtgdp*zfall;
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
        zcovptot = (double) 1.0 - (((double) 1.0 - zcovptot)*((double) 1.0 - fmax(za[jk -
           1], za[jk - 1 - 1])) / ((double) 1.0 - fmin(za[jk - 1 - 1], (double) 1.0 - 
          (double) 1.E-06)));
        zcovptot = fmax(zcovptot, (*yrecldp).rcovpmin);
        zcovpclr = fmax((double) 0.0, zcovptot - za[jk - 1]);          // clear sky proportion
        zraincld = zqxfg[3 - 1] / zcovptot;
        zsnowcld = zqxfg[4 - 1] / zcovptot;
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
      
      if (ztp1[jk - 1] <= rtt) {
        //-----------------------------------------------------
        //     Snow Autoconversion rate follow Lin et al. 1983
        //-----------------------------------------------------
        if (zicecld > zepsec) {
          
          zzco = 
            ptsphy*(*yrecldp).rsnowlin1*exp((*yrecldp).rsnowlin2*(ztp1[jk - 1] - rtt));
          
          if ((*yrecldp).laericeauto) {
            zlcrit = picrit_aer[ibl - 1][jk - 1][jl - 1];
            // 0.3 = N**0.333 with N=0.027
            zzco = zzco*(pow(((*yrecldp).rnice / pnice[ibl - 1][jk - 1][jl - 1]), (double
              ) 0.333));
          } else {
            zlcrit = (*yrecldp).rlcritsnow;
          }
          
          zsnowaut = zzco*((double) 1.0 - exp(-(pow((zicecld / zlcrit), 2))));
          zsolqb[2 - 1][4 - 1] = zsolqb[2 - 1][4 - 1] + zsnowaut;
          
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
            zlcrit = plcrit_aer[ibl - 1][jk - 1][jl - 1];
            // 0.3 = N**0.333 with N=125 cm-3
            zzco = zzco*(pow(((*yrecldp).rccn / pccn[ibl - 1][jk - 1][jl - 1]), (double) 
              0.333));
          } else {
            // Modify autoconversion threshold dependent on:
            //  land (polluted, high CCN, smaller droplets, higher threshold)
            //  sea  (clean, low CCN, larger droplets, lower threshold)
            if (plsm[ibl - 1][jl - 1] > (double) 0.5) {
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
          zprecip = 
            (zpfplsx[4 - 1][jk - 1] + zpfplsx[3 - 1][jk - 1]) / fmax(zepsec, zcovptot);
          zcfpr = (double) 1.0 + (*yrecldp).rprc1*sqrt(fmax(zprecip, (double) 0.0));
          //      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
          //       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))
          
          if ((*yrecldp).laerliqcoll) {
            // 5.0 = N**0.333 with N=125 cm-3
            zcfpr = zcfpr*(pow(((*yrecldp).rccn / pccn[ibl - 1][jk - 1][jl - 1]), (double
              ) 0.333));
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
          if (ztp1[jk - 1] <= rtt) {
            zsolqb[1 - 1][4 - 1] = zsolqb[1 - 1][4 - 1] + zrainaut;
          } else {
            zsolqb[1 - 1][3 - 1] = zsolqb[1 - 1][3 - 1] + zrainaut;
          }
          
          //--------------------------------------------------------
          //-
          //- Warm-rain process follow Khairoutdinov and Kogan (2000)
          //-
          //--------------------------------------------------------
        } else if (iwarmrain == 2) {
          
          if (plsm[ibl - 1][jl - 1] > (double) 0.5) {
            // land
            zconst = (*yrecldp).rcl_kk_cloud_num_land;
            zlcrit = (*yrecldp).rclcrit_land;
          } else {
            // ocean
            zconst = (*yrecldp).rcl_kk_cloud_num_sea;
            zlcrit = (*yrecldp).rclcrit_sea;
          }
          
          if (zliqcld > zlcrit) {
            
            zrainaut = (double) 1.5*za[jk - 1]*ptsphy*(*yrecldp).rcl_kkaau*(pow(zliqcld, 
              (*yrecldp).rcl_kkbauq))*(pow(zconst, (*yrecldp).rcl_kkbaun));
            
            zrainaut = fmin(zrainaut, zqxfg[1 - 1]);
            if (zrainaut < zepsec) {
              zrainaut = (double) 0.0;
            }
            
            zrainacc = (double) 2.0*za[jk - 1]*ptsphy*(*yrecldp)
              .rcl_kkaac*(pow((zliqcld*zraincld), (*yrecldp).rcl_kkbac));
            
            zrainacc = fmin(zrainacc, zqxfg[1 - 1]);
            if (zrainacc < zepsec) {
              zrainacc = (double) 0.0;
            }
            
          } else {
            zrainaut = (double) 0.0;
            zrainacc = (double) 0.0;
          }
          
          // If temperature < 0, then autoconversion produces snow rather than rain
          // Explicit
          if (ztp1[jk - 1] <= rtt) {
            zsolqa[1 - 1][4 - 1] = zsolqa[1 - 1][4 - 1] + zrainaut;
            zsolqa[1 - 1][4 - 1] = zsolqa[1 - 1][4 - 1] + zrainacc;
            zsolqa[4 - 1][1 - 1] = zsolqa[4 - 1][1 - 1] - zrainaut;
            zsolqa[4 - 1][1 - 1] = zsolqa[4 - 1][1 - 1] - zrainacc;
          } else {
            zsolqa[1 - 1][3 - 1] = zsolqa[1 - 1][3 - 1] + zrainaut;
            zsolqa[1 - 1][3 - 1] = zsolqa[1 - 1][3 - 1] + zrainacc;
            zsolqa[3 - 1][1 - 1] = zsolqa[3 - 1][1 - 1] - zrainaut;
            zsolqa[3 - 1][1 - 1] = zsolqa[3 - 1][1 - 1] - zrainacc;
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
        
        if (ztp1[jk - 1] <= rtt && zliqcld > zepsec) {
          
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
            
            zsolqb[1 - 1][4 - 1] = zsolqb[1 - 1][4 - 1] + zsnowrime;
            
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
      
      zicetot = zqxfg[2 - 1] + zqxfg[4 - 1];
      zmeltmax = (double) 0.0;
      
      // If there are frozen hydrometeors present and dry-bulb temperature > 0degC
      if (zicetot > zepsec && ztp1[jk - 1] > rtt) {
        
        // Calculate subsaturation
        zsubsat = fmax(zqsice[jk - 1] - zqx[5 - 1][jk - 1], (double) 0.0);
        
        // Calculate difference between dry-bulb (ZTP1) and the temperature
        // at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
        // Melting only occurs if the wet-bulb temperature >0
        // i.e. warming of ice particle due to melting > cooling
        // due to evaporation.
        ztdmtw0 = ztp1[jk - 1] - rtt - zsubsat*(ztw1 + ztw2*(pap[ibl - 1][jk - 1][jl - 1]
           - ztw3) - ztw4*(ztp1[jk - 1] - ztw5));
        // Not implicit yet...
        // Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
        zcons1 = fabs(ptsphy*((double) 1.0 + (double) 0.5*ztdmtw0) / (*yrecldp).rtaumel);
        zmeltmax = fmax(ztdmtw0*zcons1*zrldcp, (double) 0.0);
      }
      
      // Loop over frozen hydrometeors (ice, snow)
      for (jm = 1; jm <= 5; jm += 1) {
        if (iphase[jm - 1] == 2) {
          jn = imelt[jm - 1];
          if (zmeltmax > zepsec && zicetot > zepsec) {
            // Apply melting in same proportion as frozen hydrometeor fractions
            zalfa = zqxfg[jm - 1] / zicetot;
            zmelt = fmin(zqxfg[jm - 1], zalfa*zmeltmax);
            // needed in first guess
            // This implies that zqpretot has to be recalculated below
            // since is not conserved here if ice falls and liquid doesn't
            zqxfg[jm - 1] = zqxfg[jm - 1] - zmelt;
            zqxfg[jn - 1] = zqxfg[jn - 1] + zmelt;
            zsolqa[jm - 1][jn - 1] = zsolqa[jm - 1][jn - 1] + zmelt;
            zsolqa[jn - 1][jm - 1] = zsolqa[jn - 1][jm - 1] - zmelt;
          }
        }
      }
      
      //----------------------------------------------------------------------
      // 4.4b  FREEZING of RAIN
      //----------------------------------------------------------------------
      
      // If rain present
      if (zqx[3 - 1][jk - 1] > zepsec) {
        
        if (ztp1[jk - 1] <= rtt && ztp1[jk - 1 - 1] > rtt) {
          // Base of melting layer/top of refreezing layer so
          // store rain/snow fraction for precip type diagnosis
          // If mostly rain, then supercooled rain slow to freeze
          // otherwise faster to freeze (snow or ice pellets)
          zqpretot = fmax(zqx[4 - 1][jk - 1] + zqx[3 - 1][jk - 1], zepsec);
          prainfrac_toprfz[ibl - 1][jl - 1] = zqx[3 - 1][jk - 1] / zqpretot;
          if (prainfrac_toprfz[ibl - 1][jl - 1] > 0.8) {
            llrainliq = true;
          } else {
            llrainliq = false;
          }
        }
        
        // If temperature less than zero
        if (ztp1[jk - 1] < rtt) {
          
          if (prainfrac_toprfz[ibl - 1][jl - 1] > 0.8) {
            
            // Majority of raindrops completely melted
            // Refreezing is by slow heterogeneous freezing
            
            // Slope of rain particle size distribution
            zlambda = 
              pow(((*yrecldp).rcl_fac1 / (zrho*zqx[3 - 1][jk - 1])), (*yrecldp).rcl_fac2)
              ;
            
            // Calculate freezing rate based on Bigg(1953) and Wisner(1972)
            ztemp = (*yrecldp).rcl_fzrab*(ztp1[jk - 1] - rtt);
            zfrz = ptsphy*((*yrecldp).rcl_const5r / zrho)*(exp(ztemp) - (double) 1.)
              *(pow(zlambda, (*yrecldp).rcl_const6r));
            zfrzmax = fmax(zfrz, (double) 0.0);
            
          } else {
            
            // Majority of raindrops only partially melted
            // Refreeze with a shorter timescale (reverse of melting...for now)
            
            zcons1 = fabs(ptsphy*((double) 1.0 + (double) 0.5*(rtt - ztp1[jk - 1])) / 
              (*yrecldp).rtaumel);
            zfrzmax = fmax((rtt - ztp1[jk - 1])*zcons1*zrldcp, (double) 0.0);
            
          }
          
          if (zfrzmax > zepsec) {
            zfrz = fmin(zqx[3 - 1][jk - 1], zfrzmax);
            zsolqa[3 - 1][4 - 1] = zsolqa[3 - 1][4 - 1] + zfrz;
            zsolqa[4 - 1][3 - 1] = zsolqa[4 - 1][3 - 1] - zfrz;
          }
        }
        
      }
      
      
      //----------------------------------------------------------------------
      // 4.4c  FREEZING of LIQUID
      //----------------------------------------------------------------------
      // not implicit yet...
      zfrzmax = fmax(((*yrecldp).rthomo - ztp1[jk - 1])*zrldcp, (double) 0.0);
      
      jm = 1;
      jn = imelt[jm - 1];
      if (zfrzmax > zepsec && zqxfg[jm - 1] > zepsec) {
        zfrz = fmin(zqxfg[jm - 1], zfrzmax);
        zsolqa[jm - 1][jn - 1] = zsolqa[jm - 1][jn - 1] + zfrz;
        zsolqa[jn - 1][jm - 1] = zsolqa[jn - 1][jm - 1] - zfrz;
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
          fmax(zepsec, (double) 1.0 - za[jk - 1]);
        zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);
        
        zqe = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsliq[jk - 1]) / fmax(zepsec, (double) 
          1.0 - za[jk - 1]);
        //---------------------------------------------
        // humidity in moistest ZCOVPCLR part of domain
        //---------------------------------------------
        zqe = fmax((double) 0.0, fmin(zqe, zqsliq[jk - 1]));
        llo1 = zcovpclr > zepsec && zqxfg[3 - 1] > zepsec && zqe < zzrh*zqsliq[jk - 1];
        
        if (llo1) {
          // note: zpreclr is a rain flux
          zpreclr = zqxfg[3 - 1]*zcovpclr / copysign(fmax(fabs(zcovptot*zdtgdp), zepsilon
            ), zcovptot*zdtgdp);
          
          //--------------------------------------
          // actual microphysics formula in zbeta
          //--------------------------------------
          
          zbeta1 = sqrt(pap[ibl - 1][jk - 1][jl - 1] / paph[ibl - 1][klev + 1 - 1][jl - 
            1]) / (*yrecldp).rvrfactor*zpreclr / fmax(zcovpclr, zepsec);
          
          zbeta = rg*(*yrecldp).rpecons*(double) 0.5*(pow(zbeta1, (double) 0.5777));
          
          zdenom = (double) 1.0 + zbeta*ptsphy*zcorqsliq;
          zdpr = zcovpclr*zbeta*(zqsliq[jk - 1] - zqe) / zdenom*zdp*zrg_r;
          zdpevap = zdpr*zdtgdp;
          
          //---------------------------------------------------------
          // add evaporation term to explicit sink.
          // this has to be explicit since if treated in the implicit
          // term evaporation can not reduce rain to zero and model
          // produces small amounts of rainfall everywhere.
          //---------------------------------------------------------
          
          // Evaporate rain
          zevap = fmin(zdpevap, zqxfg[3 - 1]);
          
          zsolqa[3 - 1][5 - 1] = zsolqa[3 - 1][5 - 1] + zevap;
          zsolqa[5 - 1][3 - 1] = zsolqa[5 - 1][3 - 1] - zevap;
          
          //-------------------------------------------------------------
          // Reduce the total precip coverage proportional to evaporation
          // to mimic the previous scheme which had a diagnostic
          // 2-flux treatment, abandoned due to the new prognostic precip
          //-------------------------------------------------------------
          zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot - 
            za[jk - 1])*zevap / zqxfg[3 - 1]));
          
          // Update fg field
          zqxfg[3 - 1] = zqxfg[3 - 1] - zevap;
          
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
          fmax(zepsec, (double) 1.0 - za[jk - 1]);
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
        
        zqe = fmax((double) 0.0, fmin(zqx[5 - 1][jk - 1], zqsliq[jk - 1]));
        
        llo1 = zcovpclr > zepsec && zqxfg[3 - 1] > zepsec && zqe < zzrh*zqsliq[jk - 1];
        
        if (llo1) {
          
          //-------------------------------------------
          // Abel and Boutle (2012) evaporation
          //-------------------------------------------
          // Calculate local precipitation (kg/kg)
          zpreclr = zqxfg[3 - 1] / zcovptot;
          
          // Fallspeed air density correction
          zfallcorr = pow(((*yrecldp).rdensref / zrho), 0.4);
          
          // Saturation vapour pressure with respect to liquid phase
          zesatliq = (rv/rd)*(double)(r2es*exp((r3les*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4les)));  //rv / rd*foeeliq;
          
          // Slope of particle size distribution
          zlambda = pow(((*yrecldp).rcl_fac1 / (zrho*zpreclr)), (*yrecldp).rcl_fac2);            // ZPRECLR=kg/kg
          
          zevap_denom = (*yrecldp).rcl_cdenom1*zesatliq - (*yrecldp).rcl_cdenom2*ztp1[jk 
            - 1]*zesatliq + (*yrecldp).rcl_cdenom3*(pow(ztp1[jk - 1], (double) 3.))
            *pap[ibl - 1][jk - 1][jl - 1];
          
          // Temperature dependent conductivity
          zcorr2 = (pow((ztp1[jk - 1] / (double) 273.), (double) 1.5))*(double) 393. / 
            (ztp1[jk - 1] + (double) 120.);
          zka = (*yrecldp).rcl_ka273*zcorr2;
          
          zsubsat = fmax(zzrh*zqsliq[jk - 1] - zqe, (double) 0.0);
          
          zbeta = ((double) 0.5 / zqsliq[jk - 1])*(pow(ztp1[jk - 1], (double) 2.))
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
          zevap = fmin(zdpevap, zqxfg[3 - 1]);
          
          zsolqa[3 - 1][5 - 1] = zsolqa[3 - 1][5 - 1] + zevap;
          zsolqa[5 - 1][3 - 1] = zsolqa[5 - 1][3 - 1] - zevap;
          
          //-------------------------------------------------------------
          // Reduce the total precip coverage proportional to evaporation
          // to mimic the previous scheme which had a diagnostic
          // 2-flux treatment, abandoned due to the new prognostic precip
          //-------------------------------------------------------------
          zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot - 
            za[jk - 1])*zevap / zqxfg[3 - 1]));
          
          // Update fg field
          zqxfg[3 - 1] = zqxfg[3 - 1] - zevap;
          
        }
        
      }
      // on IEVAPRAIN
      
      //----------------------------------------------------------------------
      // 4.5   EVAPORATION OF SNOW
      //----------------------------------------------------------------------
      // Snow
      if (ievapsnow == 1) {
        
        zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax / 
          fmax(zepsec, (double) 1.0 - za[jk - 1]);
        zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);
        zqe = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsice[jk - 1]) / fmax(zepsec, (double) 
          1.0 - za[jk - 1]);
        
        //---------------------------------------------
        // humidity in moistest ZCOVPCLR part of domain
        //---------------------------------------------
        zqe = fmax((double) 0.0, fmin(zqe, zqsice[jk - 1]));
        llo1 = zcovpclr > zepsec && zqxfg[4 - 1] > zepsec && zqe < zzrh*zqsice[jk - 1];
        
        if (llo1) {
          // note: zpreclr is a rain flux a
          zpreclr = zqxfg[4 - 1]*zcovpclr / copysign(fmax(fabs(zcovptot*zdtgdp), zepsilon
            ), zcovptot*zdtgdp);
          
          //--------------------------------------
          // actual microphysics formula in zbeta
          //--------------------------------------
          
          zbeta1 = sqrt(pap[ibl - 1][jk - 1][jl - 1] / paph[ibl - 1][klev + 1 - 1][jl - 
            1]) / (*yrecldp).rvrfactor*zpreclr / fmax(zcovpclr, zepsec);
          
          zbeta = rg*(*yrecldp).rpecons*(pow(zbeta1, (double) 0.5777));
          
          zdenom = (double) 1.0 + zbeta*ptsphy*zcorqsice;
          zdpr = zcovpclr*zbeta*(zqsice[jk - 1] - zqe) / zdenom*zdp*zrg_r;
          zdpevap = zdpr*zdtgdp;
          
          //---------------------------------------------------------
          // add evaporation term to explicit sink.
          // this has to be explicit since if treated in the implicit
          // term evaporation can not reduce snow to zero and model
          // produces small amounts of snowfall everywhere.
          //---------------------------------------------------------
          
          // Evaporate snow
          zevap = fmin(zdpevap, zqxfg[4 - 1]);
          
          zsolqa[4 - 1][5 - 1] = zsolqa[4 - 1][5 - 1] + zevap;
          zsolqa[5 - 1][4 - 1] = zsolqa[5 - 1][4 - 1] - zevap;
          
          //-------------------------------------------------------------
          // Reduce the total precip coverage proportional to evaporation
          // to mimic the previous scheme which had a diagnostic
          // 2-flux treatment, abandoned due to the new prognostic precip
          //-------------------------------------------------------------
          zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot - 
            za[jk - 1])*zevap / zqxfg[4 - 1]));
          
          //Update first guess field
          zqxfg[4 - 1] = zqxfg[4 - 1] - zevap;
          
        }
        //---------------------------------------------------------
      } else if (ievapsnow == 2) {
        
        
        
        //-----------------------------------------------------------------------
        // Calculate relative humidity limit for snow evaporation
        //-----------------------------------------------------------------------
        zzrh = (*yrecldp).rprecrhmax + ((double) 1.0 - (*yrecldp).rprecrhmax)*zcovpmax / 
          fmax(zepsec, (double) 1.0 - za[jk - 1]);
        zzrh = fmin(fmax(zzrh, (*yrecldp).rprecrhmax), (double) 1.0);
        zqe = (zqx[5 - 1][jk - 1] - za[jk - 1]*zqsice[jk - 1]) / fmax(zepsec, (double) 
          1.0 - za[jk - 1]);
        
        //---------------------------------------------
        // humidity in moistest ZCOVPCLR part of domain
        //---------------------------------------------
        zqe = fmax((double) 0.0, fmin(zqe, zqsice[jk - 1]));
        llo1 = 
          zcovpclr > zepsec && zqx[4 - 1][jk - 1] > zepsec && zqe < zzrh*zqsice[jk - 1];
        
        if (llo1) {
          
          // Calculate local precipitation (kg/kg)
          zpreclr = zqx[4 - 1][jk - 1] / zcovptot;
          zvpice = ((double)(r2es*exp((r3ies*(ztp1[jk-1] - rtt))/(ztp1[jk-1] - r4ies)))*rv) / rd;
          
          // Particle size distribution
          // ZTCG increases Ni with colder temperatures - essentially a
          // Fletcher or Meyers scheme?
          ztcg = (double) 1.0;            //v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
          // ZFACX1I modification is based on Andrew Barrett's results
          zfacx1s = (double) 1.0;            //v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB
          
          zaplusb = (*yrecldp).rcl_apb1*zvpice - (*yrecldp).rcl_apb2*zvpice*ztp1[jk - 1] 
            + pap[ibl - 1][jk - 1][jl - 1]*(*yrecldp).rcl_apb3*(pow(ztp1[jk - 1], 3));
          zcorrfac = pow((1.0 / zrho), 0.5);
          zcorrfac2 = (pow((ztp1[jk - 1] / 273.0), 1.5))*(393.0 / (ztp1[jk - 1] + 120.0))
            ;
          
          zpr02 = zrho*zpreclr*(*yrecldp).rcl_const1s / (ztcg*zfacx1s);
          
          zterm1 = (zqsice[jk - 1] - zqe)*(pow(ztp1[jk - 1], 2))
            *zvpice*zcorrfac2*ztcg*(*yrecldp).rcl_const2s*zfacx1s / 
            (zrho*zaplusb*zqsice[jk - 1]);
          zterm2 = 0.65*(*yrecldp).rcl_const6s*(pow(zpr02, (*yrecldp).rcl_const4s)) + 
            (*yrecldp).rcl_const3s*(pow(zcorrfac, 0.5))*(pow(zrho, 0.5))*(pow(zpr02, 
            (*yrecldp).rcl_const5s)) / (pow(zcorrfac2, 0.5));
          
          zdpevap = fmax(zcovpclr*zterm1*zterm2*ptsphy, (double) 0.0);
          
          //--------------------------------------------------------------------
          // Limit evaporation to snow amount
          //--------------------------------------------------------------------
          zevap = fmin(zdpevap, zevaplimice);
          zevap = fmin(zevap, zqx[4 - 1][jk - 1]);
          
          
          zsolqa[4 - 1][5 - 1] = zsolqa[4 - 1][5 - 1] + zevap;
          zsolqa[5 - 1][4 - 1] = zsolqa[5 - 1][4 - 1] - zevap;
          
          //-------------------------------------------------------------
          // Reduce the total precip coverage proportional to evaporation
          // to mimic the previous scheme which had a diagnostic
          // 2-flux treatment, abandoned due to the new prognostic precip
          //-------------------------------------------------------------
          zcovptot = fmax((*yrecldp).rcovpmin, zcovptot - fmax((double) 0.0, (zcovptot - 
            za[jk - 1])*zevap / zqx[4 - 1][jk - 1]));
          
          //Update first guess field
          zqxfg[4 - 1] = zqxfg[4 - 1] - zevap;
          
        }
        
      }
      // on IEVAPSNOW
      
      //--------------------------------------
      // Evaporate small precipitation amounts
      //--------------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        if (llfall[jm - 1]) {
          if (zqxfg[jm - 1] < (*yrecldp).rlmin) {
            zsolqa[jm - 1][5 - 1] = zsolqa[jm - 1][5 - 1] + zqxfg[jm - 1];
            zsolqa[5 - 1][jm - 1] = zsolqa[5 - 1][jm - 1] - zqxfg[jm - 1];
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
      zanew = (za[jk - 1] + zsolac) / ((double) 1.0 + zsolab);
      zanew = fmin(zanew, (double) 1.0);
      if (zanew < (*yrecldp).ramin) {
        zanew = (double) 0.0;
      }
      zda = zanew - zaorig[jk - 1];
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
      
      for (jm = 1; jm <= 5; jm += 1) {
        for (jn = 1; jn <= 5; jn += 1) {
          llindex3[jm - 1][jn - 1] = false;
        }
        zsinksum[jm - 1] = (double) 0.0;
      }
      
      //----------------------------
      // collect sink terms and mark
      //----------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        for (jn = 1; jn <= 5; jn += 1) {
          zsinksum[jm - 1] = zsinksum[jm - 1] - zsolqa[jn - 1][jm - 1];            // +ve total is bad
        }
      }
      
      //---------------------------------------
      // calculate overshoot and scaling factor
      //---------------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        zmax = fmax(zqx[jm - 1][jk - 1], zepsec);
        zrat = fmax(zsinksum[jm - 1], zmax);
        zratio[jm - 1] = zmax / zrat;
      }
      
      //--------------------------------------------
      // scale the sink terms, in the correct order,
      // recalculating the scale factor each time
      //--------------------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        zsinksum[jm - 1] = (double) 0.0;
      }
      
      //----------------
      // recalculate sum
      //----------------
      for (jm = 1; jm <= 5; jm += 1) {
        psum_solqa = 0.0;
        for (jn = 1; jn <= 5; jn += 1) {
          psum_solqa = psum_solqa + zsolqa[jn - 1][jm - 1];
        }
        // ZSINKSUM(JL,JM)=ZSINKSUM(JL,JM)-SUM(ZSOLQA(JL,JM,1:NCLV))
        zsinksum[jm - 1] = zsinksum[jm - 1] - psum_solqa;
        //---------------------------
        // recalculate scaling factor
        //---------------------------
        zmm = fmax(zqx[jm - 1][jk - 1], zepsec);
        zrr = fmax(zsinksum[jm - 1], zmm);
        zratio[jm - 1] = zmm / zrr;
        //------
        // scale
        //------
        zzratio = zratio[jm - 1];
        //DIR$ IVDEP
        //DIR$ PREFERVECTOR
        for (jn = 1; jn <= 5; jn += 1) {
          if (zsolqa[jn - 1][jm - 1] < (double) 0.0) {
            zsolqa[jn - 1][jm - 1] = zsolqa[jn - 1][jm - 1]*zzratio;
            zsolqa[jm - 1][jn - 1] = zsolqa[jm - 1][jn - 1]*zzratio;
          }
        }
      }
      
      //--------------------------------------------------------------
      // 5.2.2 Solver
      //------------------------
      
      //------------------------
      // set the LHS of equation
      //------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        for (jn = 1; jn <= 5; jn += 1) {
          //----------------------------------------------
          // diagonals: microphysical sink terms+transport
          //----------------------------------------------
          if (jn == jm) {
            zqlhs[jm - 1][jn - 1] = (double) 1.0 + zfallsink[jm - 1];
            for (jo = 1; jo <= 5; jo += 1) {
              zqlhs[jm - 1][jn - 1] = zqlhs[jm - 1][jn - 1] + zsolqb[jn - 1][jo - 1];
            }
            //------------------------------------------
            // non-diagonals: microphysical source terms
            //------------------------------------------
          } else {
            zqlhs[jm - 1][jn - 1] = -zsolqb[jm - 1][jn - 1];              // here is the delta T - missing from doc.
          }
        }
      }
      
      //------------------------
      // set the RHS of equation
      //------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        //---------------------------------
        // sum the explicit source and sink
        //---------------------------------
        zexplicit = (double) 0.0;
        for (jn = 1; jn <= 5; jn += 1) {
          zexplicit = zexplicit + zsolqa[jn - 1][jm - 1];            // sum over middle index
        }
        zqxn[jm - 1] = zqx[jm - 1][jk - 1] + zexplicit;
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
      for (jn = 1; jn <= 5 - 1; jn += 1) {
        // number of steps
        for (jm = jn + 1; jm <= 5; jm += 1) {
          // row index
          zqlhs[jn - 1][jm - 1] = zqlhs[jn - 1][jm - 1] / zqlhs[jn - 1][jn - 1];
          for (ik = jn + 1; ik <= 5; ik += 1) {
            // column index
            zqlhs[ik - 1][jm - 1] = 
              zqlhs[ik - 1][jm - 1] - zqlhs[jn - 1][jm - 1]*zqlhs[ik - 1][jn - 1];
          }
        }
      }
      
      // Backsubstitution
      //  step 1
      for (jn = 2; jn <= 5; jn += 1) {
        for (jm = 1; jm <= jn - 1; jm += 1) {
          zqxn[jn - 1] = zqxn[jn - 1] - zqlhs[jm - 1][jn - 1]*zqxn[jm - 1];
        }
      }
      //  step 2
      zqxn[5 - 1] = zqxn[5 - 1] / zqlhs[5 - 1][5 - 1];
      for (jn = 5 - 1; jn >= 1; jn += -1) {
        for (jm = jn + 1; jm <= 5; jm += 1) {
          zqxn[jn - 1] = zqxn[jn - 1] - zqlhs[jm - 1][jn - 1]*zqxn[jm - 1];
        }
        zqxn[jn - 1] = zqxn[jn - 1] / zqlhs[jn - 1][jn - 1];
      }
      
      // Ensure no small values (including negatives) remain in cloud variables nor
      // precipitation rates.
      // Evaporate l,i,r,s to water vapour. Latent heating taken into account below
      for (jn = 1; jn <= 5 - 1; jn += 1) {
        if (zqxn[jn - 1] < zepsec) {
          zqxn[5 - 1] = zqxn[5 - 1] + zqxn[jn - 1];
          zqxn[jn - 1] = (double) 0.0;
        }
      }
      
      //--------------------------------
      // variables needed for next level
      //--------------------------------
      for (jm = 1; jm <= 5; jm += 1) {
        zqxnm1[jm - 1] = zqxn[jm - 1];
        zqxn2d[jm - 1][jk - 1] = zqxn[jm - 1];
      }
      
      //------------------------------------------------------------------------
      // 5.3 Precipitation/sedimentation fluxes to next level
      //     diagnostic precipitation fluxes
      //     It is this scaled flux that must be used for source to next layer
      //------------------------------------------------------------------------
      
      for (jm = 1; jm <= 5; jm += 1) {
        zpfplsx[jm - 1][jk + 1 - 1] = zfallsink[jm - 1]*zqxn[jm - 1]*zrdtgdp;
      }
      
      // Ensure precipitation fraction is zero if no precipitation
      zqpretot = zpfplsx[4 - 1][jk + 1 - 1] + zpfplsx[3 - 1][jk + 1 - 1];
      if (zqpretot < zepsec) {
        zcovptot = (double) 0.0;
      }
      
      //######################################################################
      //              6  *** UPDATE TENDANCIES ***
      //######################################################################
      
      //--------------------------------
      // 6.1 Temperature and CLV budgets
      //--------------------------------
      
      for (jm = 1; jm <= 5 - 1; jm += 1) {
        
        // calculate fluxes in and out of box for conservation of TL
        zfluxq[jm - 1] = zpsupsatsrce[jm - 1] + zconvsrce[jm - 1] + zfallsrce[jm - 1] - 
          (zfallsink[jm - 1] + zconvsink[jm - 1])*zqxn[jm - 1];
        
        if (iphase[jm - 1] == 1) {
          tendency_loc_t[ibl - 1][jk - 1][jl - 1] = tendency_loc_t[ibl - 1][jk - 1][jl - 
            1] + ralvdcp*(zqxn[jm - 1] - zqx[jm - 1][jk - 1] - zfluxq[jm - 1])*zqtmst;
        }
        
        if (iphase[jm - 1] == 2) {
          tendency_loc_t[ibl - 1][jk - 1][jl - 1] = tendency_loc_t[ibl - 1][jk - 1][jl - 
            1] + ralsdcp*(zqxn[jm - 1] - zqx[jm - 1][jk - 1] - zfluxq[jm - 1])*zqtmst;
        }
        
        //----------------------------------------------------------------------
        // New prognostic tendencies - ice,liquid rain,snow
        // Note: CLV arrays use PCLV in calculation of tendency while humidity
        //       uses ZQX. This is due to clipping at start of cloudsc which
        //       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
        //----------------------------------------------------------------------
        tendency_loc_cld[ibl - 1][jm - 1][jk - 1][jl - 1] = tendency_loc_cld[ibl - 1][jm 
          - 1][jk - 1][jl - 1] + (zqxn[jm - 1] - zqx0[jm - 1][jk - 1])*zqtmst;
        
      }
      
      //----------------------
      // 6.2 Humidity budget
      //----------------------
      tendency_loc_q[ibl - 1][jk - 1][jl - 1] = tendency_loc_q[ibl - 1][jk - 1][jl - 1] +
         (zqxn[5 - 1] - zqx[5 - 1][jk - 1])*zqtmst;
      
      //-------------------
      // 6.3 cloud cover
      //-----------------------
      tendency_loc_a[ibl - 1][jk - 1][jl - 1] = 
        tendency_loc_a[ibl - 1][jk - 1][jl - 1] + zda*zqtmst;
      
      //--------------------------------------------------
      // Copy precipitation fraction into output variable
      //-------------------------------------------------
      pcovptot[ibl - 1][jk - 1][jl - 1] = zcovptot;
      
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
    for (jk = 1; jk <= klev + 1; jk += 1) {
      pfplsl[ibl - 1][jk - 1][jl - 1] = zpfplsx[3 - 1][jk - 1] + zpfplsx[1 - 1][jk - 1];
      pfplsn[ibl - 1][jk - 1][jl - 1] = zpfplsx[4 - 1][jk - 1] + zpfplsx[2 - 1][jk - 1];
    }
    
    //--------
    // Fluxes:
    //--------
    pfsqlf[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfsqif[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfsqrf[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfsqsf[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfcqlng[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfcqnng[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfcqrng[ibl - 1][1 - 1][jl - 1] = (double) 0.0;      //rain
    pfcqsng[ibl - 1][1 - 1][jl - 1] = (double) 0.0;      //snow
    // fluxes due to turbulence
    pfsqltur[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    pfsqitur[ibl - 1][1 - 1][jl - 1] = (double) 0.0;
    
    for (jk = 1; jk <= klev; jk += 1) {
      
      zgdph_r = 
        -zrg_r*(paph[ibl - 1][jk + 1 - 1][jl - 1] - paph[ibl - 1][jk - 1][jl - 1])*zqtmst
        ;
      pfsqlf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqlf[ibl - 1][jk - 1][jl - 1];
      pfsqif[ibl - 1][jk + 1 - 1][jl - 1] = pfsqif[ibl - 1][jk - 1][jl - 1];
      pfsqrf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqlf[ibl - 1][jk - 1][jl - 1];
      pfsqsf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqif[ibl - 1][jk - 1][jl - 1];
      pfcqlng[ibl - 1][jk + 1 - 1][jl - 1] = pfcqlng[ibl - 1][jk - 1][jl - 1];
      pfcqnng[ibl - 1][jk + 1 - 1][jl - 1] = pfcqnng[ibl - 1][jk - 1][jl - 1];
      pfcqrng[ibl - 1][jk + 1 - 1][jl - 1] = pfcqlng[ibl - 1][jk - 1][jl - 1];
      pfcqsng[ibl - 1][jk + 1 - 1][jl - 1] = pfcqnng[ibl - 1][jk - 1][jl - 1];
      pfsqltur[ibl - 1][jk + 1 - 1][jl - 1] = pfsqltur[ibl - 1][jk - 1][jl - 1];
      pfsqitur[ibl - 1][jk + 1 - 1][jl - 1] = pfsqitur[ibl - 1][jk - 1][jl - 1];
      
      zalfaw = zfoealfa[jk - 1];
      
      // Liquid , LS scheme minus detrainment
      pfsqlf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqlf[ibl - 1][jk + 1 - 1][jl - 1] + 
        (zqxn2d[1 - 1][jk - 1] - zqx0[1 - 1][jk - 1] + pvfl[ibl - 1][jk - 1][jl - 
        1]*ptsphy - zalfaw*plude[ibl - 1][jk - 1][jl - 1])*zgdph_r;
      // liquid, negative numbers
      pfcqlng[ibl - 1][jk + 1 - 1][jl - 1] = 
        pfcqlng[ibl - 1][jk + 1 - 1][jl - 1] + zlneg[1 - 1][jk - 1]*zgdph_r;
      
      // liquid, vertical diffusion
      pfsqltur[ibl - 1][jk + 1 - 1][jl - 1] = pfsqltur[ibl - 1][jk + 1 - 1][jl - 1] + 
        pvfl[ibl - 1][jk - 1][jl - 1]*ptsphy*zgdph_r;
      
      // Rain, LS scheme
      pfsqrf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqrf[ibl - 1][jk + 1 - 1][jl - 1] + 
        (zqxn2d[3 - 1][jk - 1] - zqx0[3 - 1][jk - 1])*zgdph_r;
      // rain, negative numbers
      pfcqrng[ibl - 1][jk + 1 - 1][jl - 1] = 
        pfcqrng[ibl - 1][jk + 1 - 1][jl - 1] + zlneg[3 - 1][jk - 1]*zgdph_r;
      
      // Ice , LS scheme minus detrainment
      pfsqif[ibl - 1][jk + 1 - 1][jl - 1] = pfsqif[ibl - 1][jk + 1 - 1][jl - 1] + 
        (zqxn2d[2 - 1][jk - 1] - zqx0[2 - 1][jk - 1] + pvfi[ibl - 1][jk - 1][jl - 
        1]*ptsphy - ((double) 1.0 - zalfaw)*plude[ibl - 1][jk - 1][jl - 1])*zgdph_r;
      // ice, negative numbers
      pfcqnng[ibl - 1][jk + 1 - 1][jl - 1] = 
        pfcqnng[ibl - 1][jk + 1 - 1][jl - 1] + zlneg[2 - 1][jk - 1]*zgdph_r;
      
      // ice, vertical diffusion
      pfsqitur[ibl - 1][jk + 1 - 1][jl - 1] = pfsqitur[ibl - 1][jk + 1 - 1][jl - 1] + 
        pvfi[ibl - 1][jk - 1][jl - 1]*ptsphy*zgdph_r;
      
      // snow, LS scheme
      pfsqsf[ibl - 1][jk + 1 - 1][jl - 1] = pfsqsf[ibl - 1][jk + 1 - 1][jl - 1] + 
        (zqxn2d[4 - 1][jk - 1] - zqx0[4 - 1][jk - 1])*zgdph_r;
      // snow, negative numbers
      pfcqsng[ibl - 1][jk + 1 - 1][jl - 1] = 
        pfcqsng[ibl - 1][jk + 1 - 1][jl - 1] + zlneg[4 - 1][jk - 1]*zgdph_r;
    }
    
    //-----------------------------------
    // enthalpy flux due to precipitation
    //-----------------------------------
    for (jk = 1; jk <= klev + 1; jk += 1) {
      pfhpsl[ibl - 1][jk - 1][jl - 1] = -rlvtt*pfplsl[ibl - 1][jk - 1][jl - 1];
      pfhpsn[ibl - 1][jk - 1][jl - 1] = -rlstt*pfplsn[ibl - 1][jk - 1][jl - 1];
    }
    
    //===============================================================================
    //IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)
  }
}

