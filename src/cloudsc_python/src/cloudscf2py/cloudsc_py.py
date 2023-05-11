import numpy as np
def cloudsc_py(kidia: np.int32, kfdia: np.int32, klon: np.int32, klev: np.int32, ptsphy: np.float64, pt: np.ndarray, pq: np.ndarray, tendency_tmp_t: np.ndarray, tendency_tmp_q: np.ndarray, tendency_tmp_a: np.ndarray, tendency_tmp_cld: np.ndarray, tendency_loc_t: np.ndarray, 
  tendency_loc_q: np.ndarray, tendency_loc_a: np.ndarray, tendency_loc_cld: np.ndarray, pvfa: np.ndarray, pvfl: np.ndarray, pvfi: np.ndarray, pdyna: np.ndarray, pdynl: np.ndarray, pdyni: np.ndarray, phrsw: np.ndarray, phrlw: np.ndarray, pvervel: np.ndarray, pap: np.ndarray, paph: np.ndarray, 
  plsm: np.ndarray, ldcum: np.ndarray, ktype: np.ndarray, plu: np.ndarray, plude: np.ndarray, psnde: np.ndarray, pmfu: np.ndarray, pmfd: np.ndarray, pa: np.ndarray, pclv: np.ndarray, psupsat: np.ndarray, plcrit_aer: np.ndarray, picrit_aer: np.ndarray, pre_ice: np.ndarray, pccn: np.ndarray, 
  pnice: np.ndarray, pcovptot: np.ndarray, prainfrac_toprfz: np.ndarray, pfsqlf: np.ndarray, pfsqif: np.ndarray, pfcqnng: np.ndarray, pfcqlng: np.ndarray, pfsqrf: np.ndarray, pfsqsf: np.ndarray, pfcqrng: np.ndarray, pfcqsng: np.ndarray, pfsqltur: np.ndarray, pfsqitur: np.ndarray, 
  pfplsl: np.ndarray, pfplsn: np.ndarray, pfhpsl: np.ndarray, pfhpsn: np.ndarray, ydcst, ydthf, yrecldp):
  # USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV
  # USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
  #  & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
  #  & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
  # USE YOECLDP  , ONLY : TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  # USE YOECLDP  , ONLY : TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  
  # USE FCTTRE_MOD, ONLY: FOEDELTA, FOEALFA, FOEEWM, FOEEICE, FOEELIQ, FOELDCP, FOELDCPM, FOEDEM
  # USE FCCLD_MOD, ONLY : FOKOOP
  
  
  
  #-------------------------------------------------------------------------------
  #                 Declare input/output arguments
  #-------------------------------------------------------------------------------
  
  # number of microphysics variables
  nclv = 5
  # liquid cloud water
  ncldql = 1
  # ice cloud water
  ncldqi = 2
  # rain water
  ncldqr = 3
  # snow
  ncldqs = 4
  # vapour
  ncldqv = 5
  
  
  # PLCRIT_AER : critical liquid mmr for rain autoconversion process
  # PICRIT_AER : critical liquid mmr for snow autoconversion process
  # PRE_LIQ : liq Re
  # PRE_ICE : ice Re
  # PCCN    : liquid cloud condensation nuclei
  # PNICE   : ice number concentration (cf. CCN)
  
  
  # Number of grid points
  # Number of levels
  # Physics timestep
  # Land fraction (0-1)
  # Convection active
  # Convection type 0,1,2
  
  
  # Supersat clipped at previous time level in SLTEND
  # Flux diagnostics for DDH budget
  
  # TYPE(TECLDP), INTENT(INOUT) :: YRECLDP
  
  
  #-------------------------------------------------------------------------------
  #                       Declare local variables
  #-------------------------------------------------------------------------------
  
  zlcond1 = np.ndarray(order="F", shape=(klon,))
  zlcond2 = np.ndarray(order="F", shape=(klon,))
  zlevapl = np.ndarray(order="F", shape=(klon,))
  zlevapi = np.ndarray(order="F", shape=(klon,))
  zrainaut = np.ndarray(order="F", shape=(klon,))
  zsnowaut = np.ndarray(order="F", shape=(klon,))
  zliqcld = np.ndarray(order="F", shape=(klon,))
  zicecld = np.ndarray(order="F", shape=(klon,))
  #  condensation and evaporation terms
  # autoconversion terms
  zfokoop = np.ndarray(order="F", shape=(klon,))
  # number concentration of ice nuclei
  zicenuclei = np.ndarray(order="F", shape=(klon,))
  
  zlicld = np.ndarray(order="F", shape=(klon,))
  zlfinalsum = np.ndarray(order="F", shape=(klon,))
  zdqs = np.ndarray(order="F", shape=(klon,))
  ztold = np.ndarray(order="F", shape=(klon,))
  zqold = np.ndarray(order="F", shape=(klon,))
  zdtgdp = np.ndarray(order="F", shape=(klon,))
  zrdtgdp = np.ndarray(order="F", shape=(klon,))
  ztrpaus = np.ndarray(order="F", shape=(klon,))
  zcovpclr = np.ndarray(order="F", shape=(klon,))
  zcovptot = np.ndarray(order="F", shape=(klon,))
  zcovpmax = np.ndarray(order="F", shape=(klon,))
  zqpretot = np.ndarray(order="F", shape=(klon,))
  zldefr = np.ndarray(order="F", shape=(klon,))
  zldifdt = np.ndarray(order="F", shape=(klon,))
  zdtgdpf = np.ndarray(order="F", shape=(klon,))
  zacust = np.ndarray(order="F", shape=(klon,))
  zmf = np.ndarray(order="F", shape=(klon,))
  
  zrho = np.ndarray(order="F", shape=(klon,))
  ztmp1 = np.ndarray(order="F", shape=(klon,))
  ztmp2 = np.ndarray(order="F", shape=(klon,))
  ztmp3 = np.ndarray(order="F", shape=(klon,))
  ztmp4 = np.ndarray(order="F", shape=(klon,))
  ztmp5 = np.ndarray(order="F", shape=(klon,))
  ztmp6 = np.ndarray(order="F", shape=(klon,))
  ztmp7 = np.ndarray(order="F", shape=(klon,))
  zalfawm = np.ndarray(order="F", shape=(klon,))
  
  # Accumulators of A,B,and C factors for cloud equations
  # -ve implicit CC
  zsolab = np.ndarray(order="F", shape=(klon,))
  # linear CC
  zsolac = np.ndarray(order="F", shape=(klon,))
  zanewm1 = np.ndarray(order="F", shape=(klon,))
  
  zgdp = np.ndarray(order="F", shape=(klon,))
  
  #---for flux calculation
  zda = np.ndarray(order="F", shape=(klon,))
  
  llflag = np.ndarray(order="F", shape=(klon,))
  
  
  zdp = np.ndarray(order="F", shape=(klon,))
  zpaphd = np.ndarray(order="F", shape=(klon,))
  
  # & ZALFACU, ZALFALS
  #REAL(KIND=JPRB) :: ZBOTT
  zmin = np.ndarray(order="F", shape=(klon,))
  zsupsat = np.ndarray(order="F", shape=(klon,))
  
  #----------------------------
  # Arrays for new microphysics
  #----------------------------
  # marker for water phase of each species
  iphase = np.ndarray(order="F", shape=(nclv,))
  # 0=vapour, 1=liquid, 2=ice
  
  # marks melting linkage for ice categories
  imelt = np.ndarray(order="F", shape=(nclv,))
  # ice->liquid, snow->rain
  
  # marks falling species
  llfall = np.ndarray(order="F", shape=(nclv,))
  # LLFALL=0, cloud cover must > 0 for zqx > 0
  # LLFALL=1, no cloud needed, zqx can evaporate
  
  
  # Keep the following for possible future total water variance scheme?
  #REAL(KIND=JPRB) :: ZTL(KLON,KLEV)       ! liquid water temperature
  #REAL(KIND=JPRB) :: ZABETA(KLON,KLEV)    ! cloud fraction
  #REAL(KIND=JPRB) :: ZVAR(KLON,KLEV)      ! temporary variance
  #REAL(KIND=JPRB) :: ZQTMIN(KLON,KLEV)
  #REAL(KIND=JPRB) :: ZQTMAX(KLON,KLEV)
  
  zmeltmax = np.ndarray(order="F", shape=(klon,))
  zfrzmax = np.ndarray(order="F", shape=(klon,))
  zicetot = np.ndarray(order="F", shape=(klon,))
  
  
  #REAL(KIND=JPRB) :: ZQSBIN(KLON,KLEV) ! binary switched ice/liq saturation
  
  #REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
  #REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
  #REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice
  
  #REAL(KIND=JPRB) :: ZFOEEICET(KLON,KLEV)
  
  zdqsliqdt = np.ndarray(order="F", shape=(klon,))
  zdqsicedt = np.ndarray(order="F", shape=(klon,))
  zdqsmixdt = np.ndarray(order="F", shape=(klon,))
  zcorqsliq = np.ndarray(order="F", shape=(klon,))
  zcorqsice = np.ndarray(order="F", shape=(klon,))
  #REAL(KIND=JPRB) :: ZCORQSBIN(KLON)
  zcorqsmix = np.ndarray(order="F", shape=(klon,))
  zevaplimliq = np.ndarray(order="F", shape=(klon,))
  zevaplimice = np.ndarray(order="F", shape=(klon,))
  zevaplimmix = np.ndarray(order="F", shape=(klon,))
  
  #-------------------------------------------------------
  # SOURCE/SINK array for implicit and explicit terms
  #-------------------------------------------------------
  # a POSITIVE value entered into the arrays is a...
  #            Source of this variable
  #            |
  #            |   Sink of this variable
  #            |   |
  #            V   V
  # ZSOLQA(JL,IQa,IQb)  = explicit terms
  # ZSOLQB(JL,IQa,IQb)  = implicit terms
  # Thus if ZSOLAB(JL,NCLDQL,IQV)=K where K>0 then this is
  # a source of NCLDQL and a sink of IQV
  # put 'magic' source terms such as PLUDE from
  # detrainment into explicit source/sink array diagnognal
  # ZSOLQA(NCLDQL,NCLDQL)= -PLUDE
  # i.e. A positive value is a sink!????? weird...
  #-------------------------------------------------------
  
  # e.g. microphysical pathways between ice variables.
  # fall speeds of three categories
  zvqx = np.ndarray(order="F", shape=(nclv,))
  
  # for sedimentation source/sink terms
  
  # for convection detrainment source and subsidence source/sink terms
  
  # for supersaturation source term from previous timestep
  
  # Numerical fit to wet bulb temperature
  ztw1 = 1329.31
  ztw2 = 0.0074615
  ztw3 = 0.85E5
  ztw4 = 40.637
  ztw5 = 275.0
  
  # Subsaturation for snow melting term
  # Diff between dry-bulb temperature and
  # temperature when wet-bulb = 0degC
  
  # Variables for deposition term
  # Temperature dependent function for ice PSD
  # PSD correction factor
  # for ice dep
  # Distance from cloud top
  zcldtopdist = np.ndarray(order="F", shape=(klon,))
  # No. of ice nuclei factor for deposition
  
  # Autoconversion/accretion/riming/evaporation
  zrainacc = np.ndarray(order="F", shape=(klon,))
  zraincld = np.ndarray(order="F", shape=(klon,))
  zsnowrime = np.ndarray(order="F", shape=(klon,))
  zsnowcld = np.ndarray(order="F", shape=(klon,))
  
  # Rain freezing
  # True if majority of raindrops are liquid (no ice core)
  llrainliq = np.ndarray(order="F", shape=(klon,))
  
  #----------------------------
  # End: new microphysics
  #----------------------------
  
  #----------------------
  # SCM budget statistics
  #----------------------
  
  
  zrg = np.ndarray(order="F", shape=(klon,))
  
  
  
  psum_solqa = np.ndarray(order="F", shape=(klon,))
  
  # #include "fcttre.func.h"
  # #include "fccld.func.h"
  #*
  #     ------------------------------------------------------------------
  
  #     This COMDECK includes the Thermodynamical functions for the cy39
  #       ECMWF Physics package.
  #       Consistent with YOMCST Basic physics constants, assuming the
  #       partial pressure of water vapour is given by a first order
  #       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
  #       in YOETHF
  #       Two sets of functions are available. In the first set only the
  #       cases water or ice are distinguished by temperature.  This set
  #       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
  #       The second set considers, besides the two cases water and ice
  #       also a mix of both for the temperature range  YDTHF% RTICE < T <  YDTHF% RTWAT.
  #       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
  #       FKOOP modifies the ice saturation mixing ratio for homogeneous
  #       nucleation. FOE_DEWM_DT provides an approximate first derivative
  #       of FOEEWM.
  
  #       Depending on the consideration of mixed phases either the first
  #       set (e.g. surface, post-processing) or the second set
  #       (e.g. clouds, condensation, convection) should be used.
  
  #     ------------------------------------------------------------------
  #     *****************************************************************
  
  #                NO CONSIDERATION OF MIXED PHASES
  
  #     *****************************************************************
  def foedelta(ptare):
    return max(0.0, 1.0*np.sign(ptare - ydcst.rtt))
  
  #                  FOEDELTA = 1    water
  #                  FOEDELTA = 0    ice
  
  #     THERMODYNAMICAL FUNCTIONS .
  
  #     Pressure of water vapour at saturation
  #        INPUT : PTARE = TEMPERATURE
  def foeew(ptare):
    return ydthf.r2es*np.exp((ydthf.r3les*foedelta(ptare) + ydthf.r3ies*(1.0 - foedelta(ptare)))*(ptare - ydcst.rtt) / (ptare - (ydthf.r4les*foedelta(ptare) + ydthf.r4ies*(1.0 - foedelta(ptare)))))
  
  def foede(ptare):
    return (foedelta(ptare)*ydthf.r5alvcp + (1.0 - foedelta(ptare))*ydthf.r5alscp) / (ptare - (ydthf.r4les*foedelta(ptare) + ydthf.r4ies*(1.0 - foedelta(ptare))))**2
  
  def foedesu(ptare):
    return (foedelta(ptare)*ydthf.r5les + (1.0 - foedelta(ptare))*ydthf.r5ies) / (ptare - (ydthf.r4les*foedelta(ptare) + ydthf.r4ies*(1.0 - foedelta(ptare))))**2
  
  def foelh(ptare):
    return foedelta(ptare)*ydcst.rlvtt + (1.0 - foedelta(ptare))*ydcst.rlstt
  
  def foeldcp(ptare):
    return foedelta(ptare)*ydthf.ralvdcp + (1.0 - foedelta(ptare))*ydthf.ralsdcp
  
  #     *****************************************************************
  
  #           CONSIDERATION OF MIXED PHASES
  
  #     *****************************************************************
  
  #     FOEALFA is calculated to distinguish the three cases:
  
  #                       FOEALFA=1            water phase
  #                       FOEALFA=0            ice phase
  #                       0 < FOEALFA < 1      mixed phase
  
  #               INPUT : PTARE = TEMPERATURE
  def foealfa(ptare):
    return min(1.0, ((max(ydthf.rtice, min(ydthf.rtwat, ptare)) - ydthf.rtice)*ydthf.rtwat_rtice_r)**2)
  
  
  #     Pressure of water vapour at saturation
  #        INPUT : PTARE = TEMPERATURE
  def foeewm(ptare):
    return ydthf.r2es*(foealfa(ptare)*np.exp(ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les)) + (1.0 - foealfa(ptare))*np.exp(ydthf.r3ies*(ptare - ydcst.rtt) / (ptare - ydthf.r4ies)))
  
  def foe_dewm_dt(ptare):
    return ydthf.r2es*(ydthf.r3les*foealfa(ptare)*np.exp(ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les))*(ydcst.rtt - ydthf.r4les) / (ptare - ydthf.r4les)**2 + ydthf.r3ies*(1.0 - foealfa(ptare))*np.exp(ydthf.r3ies*(ptare - ydcst.rtt) / (ptare - ydthf.r4ies))*(ydcst.rtt - ydthf.r4ies) / 
      (ptare - ydthf.r4ies)**2)
  
  def foedem(ptare):
    return foealfa(ptare)*ydthf.r5alvcp*(1.0 / (ptare - ydthf.r4les)**2) + (1.0 - foealfa(ptare))*ydthf.r5alscp*(1.0 / (ptare - ydthf.r4ies)**2)
  
  def foeldcpm(ptare):
    return foealfa(ptare)*ydthf.ralvdcp + (1.0 - foealfa(ptare))*ydthf.ralsdcp
  
  def foelhm(ptare):
    return foealfa(ptare)*ydcst.rlvtt + (1.0 - foealfa(ptare))*ydcst.rlstt
  
  
  #     Temperature normalization for humidity background change of variable
  #        INPUT : PTARE = TEMPERATURE
  def foetb(ptare):
    return foealfa(ptare)*ydthf.r3les*(ydcst.rtt - ydthf.r4les)*(1.0 / (ptare - ydthf.r4les)**2) + (1.0 - foealfa(ptare))*ydthf.r3ies*(ydcst.rtt - ydthf.r4ies)*(1.0 / (ptare - ydthf.r4ies)**2)
  
  #     ------------------------------------------------------------------
  #     *****************************************************************
  
  #           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV
  
  #     *****************************************************************
  
  #     FOEALFCU is calculated to distinguish the three cases:
  
  #                       FOEALFCU=1            water phase
  #                       FOEALFCU=0            ice phase
  #                       0 < FOEALFCU < 1      mixed phase
  
  #               INPUT : PTARE = TEMPERATURE
  def foealfcu(ptare):
    return min(1.0, ((max(ydthf.rticecu, min(ydthf.rtwat, ptare)) - ydthf.rticecu)*ydthf.rtwat_rticecu_r)**2)
  
  
  #     Pressure of water vapour at saturation
  #        INPUT : PTARE = TEMPERATURE
  def foeewmcu(ptare):
    return ydthf.r2es*(foealfcu(ptare)*np.exp(ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les)) + (1.0 - foealfcu(ptare))*np.exp(ydthf.r3ies*(ptare - ydcst.rtt) / (ptare - ydthf.r4ies)))
  
  def foedemcu(ptare):
    return foealfcu(ptare)*ydthf.r5alvcp*(1.0 / (ptare - ydthf.r4les)**2) + (1.0 - foealfcu(ptare))*ydthf.r5alscp*(1.0 / (ptare - ydthf.r4ies)**2)
  
  def foeldcpmcu(ptare):
    return foealfcu(ptare)*ydthf.ralvdcp + (1.0 - foealfcu(ptare))*ydthf.ralsdcp
  
  def foelhmcu(ptare):
    return foealfcu(ptare)*ydcst.rlvtt + (1.0 - foealfcu(ptare))*ydcst.rlstt
  #     ------------------------------------------------------------------
  
  #     Pressure of water vapour at saturation
  #     This one is for the WMO definition of saturation, i.e. always
  #     with respect to water.
  #
  #     Duplicate to FOEELIQ and FOEEICE for separate ice variable
  #     FOEELIQ always respect to water
  #     FOEEICE always respect to ice
  #     (could use FOEEW and FOEEWMO, but naming convention unclear)
  #     FOELSON returns e wrt liquid water using D Sonntag (1994, Met. Zeit.)
  #      - now recommended for use with radiosonde data (WMO CIMO guide, 2014)
  #      unlike the FOEE functions does not include 1/( YDCST% RETV+1.0_JPRB) factor
  
  def foeewmo(ptare):
    return ydthf.r2es*np.exp(ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les))
  def foeeliq(ptare):
    return ydthf.r2es*np.exp(ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les))
  def foeeice(ptare):
    return ydthf.r2es*np.exp(ydthf.r3ies*(ptare - ydcst.rtt) / (ptare - ydthf.r4ies))
  def foelson(ptare):
    return np.exp(-6096.9385 / ptare + 21.2409642 - 2.711193E-2*ptare + 1.673952E-5*ptare**2 + 2.433502*log(ptare))
  
  def foeles_v(ptare):
    return ydthf.r3les*(ptare - ydcst.rtt) / (ptare - ydthf.r4les)
  def foeies_v(ptare):
    return ydthf.r3ies*(ptare - ydcst.rtt) / (ptare - ydthf.r4ies)
  def foeewm_v(ptare, exp1, exp2):
    return ydthf.r2es*(foealfa(ptare)*exp1 + (1.0 - foealfa(ptare))*exp2)
  def foeewmcu_v(ptare, exp1, exp2):
    return ydthf.r2es*(foealfcu(ptare)*exp1 + (1.0 - foealfcu(ptare))*exp2)
  # (C) Copyright 1988- ECMWF.
  #
  # This software is licensed under the terms of the Apache Licence Version 2.0
  # which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
  #
  # In applying this licence, ECMWF does not waive the privileges and immunities
  # granted to it by virtue of its status as an intergovernmental organisation
  # nor does it submit to any jurisdiction.
  
  #*
  #     ------------------------------------------------------------------
  #     This COMDECK defines functions to be used in the cloud scheme
  #       other than the standard saturation vapour pressure
  #
  #       FKOOP modifies the ice saturation mixing ratio for homogeneous
  #       nucleation
  #
  #     note: PTARE is temperature and is definited in frttre.h
  #           which MUST be included before this function block
  #
  #     **********************************************
  #     KOOP formula for homogeneous nucleation of ice
  #     **********************************************
  #
  #               INPUT : PTARE = TEMPERATURE
  def fokoop(ptare):
    return min(ydthf.rkoop1 - ydthf.rkoop2*ptare, foeeliq(ptare) / foeeice(ptare))
  #===============================================================================
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)
  zfoealfa = np.ndarray(order="F", shape=(klev + 1, klon,))
  ztp1 = np.ndarray(order="F", shape=(klev, klon,))
  zlcust = np.ndarray(order="F", shape=(nclv, klon,))
  zli = np.ndarray(order="F", shape=(klev, klon,))
  za = np.ndarray(order="F", shape=(klev, klon,))
  zaorig = np.ndarray(order="F", shape=(klev, klon,))
  llindex1 = np.ndarray(order="F", shape=(nclv, klon,))
  llindex3 = np.ndarray(order="F", shape=(nclv, nclv, klon,))
  iorder = np.ndarray(order="F", shape=(nclv, klon,))
  zliqfrac = np.ndarray(order="F", shape=(klev, klon,))
  zicefrac = np.ndarray(order="F", shape=(klev, klon,))
  zqx = np.ndarray(order="F", shape=(nclv, klev, klon,))
  zqx0 = np.ndarray(order="F", shape=(nclv, klev, klon,))
  zqxn = np.ndarray(order="F", shape=(nclv, klon,))
  zqxfg = np.ndarray(order="F", shape=(nclv, klon,))
  zqxnm1 = np.ndarray(order="F", shape=(nclv, klon,))
  zfluxq = np.ndarray(order="F", shape=(nclv, klon,))
  zpfplsx = np.ndarray(order="F", shape=(nclv, klev + 1, klon,))
  zlneg = np.ndarray(order="F", shape=(nclv, klev, klon,))
  zqxn2d = np.ndarray(order="F", shape=(nclv, klev, klon,))
  zqsmix = np.ndarray(order="F", shape=(klev, klon,))
  zqsliq = np.ndarray(order="F", shape=(klev, klon,))
  zqsice = np.ndarray(order="F", shape=(klev, klon,))
  zfoeewmt = np.ndarray(order="F", shape=(klev, klon,))
  zfoeew = np.ndarray(order="F", shape=(klev, klon,))
  zfoeeliqt = np.ndarray(order="F", shape=(klev, klon,))
  zsolqa = np.ndarray(order="F", shape=(nclv, nclv, klon,))
  zsolqb = np.ndarray(order="F", shape=(nclv, nclv, klon,))
  zqlhs = np.ndarray(order="F", shape=(nclv, nclv, klon,))
  zratio = np.ndarray(order="F", shape=(nclv, klon,))
  zsinksum = np.ndarray(order="F", shape=(nclv, klon,))
  zfallsink = np.ndarray(order="F", shape=(nclv, klon,))
  zfallsrce = np.ndarray(order="F", shape=(nclv, klon,))
  zconvsrce = np.ndarray(order="F", shape=(nclv, klon,))
  zconvsink = np.ndarray(order="F", shape=(nclv, klon,))
  zpsupsatsrce = np.ndarray(order="F", shape=(nclv, klon,))
  
  # YDCST,    YDTHF
  
  
  
  
  
  #===============================================================================
  #  0.0     Beginning of timestep book-keeping
  #----------------------------------------------------------------------
  
  
  #######################################################################
  #             0.  *** SET UP CONSTANTS ***
  #######################################################################
  
  # ZEPSILON=100._JPRB*EPSILON(ZEPSILON)
  zepsilon = 1.E-14
  
  # ---------------------------------------------------------------------
  # Set version of warm-rain autoconversion/accretion
  # IWARMRAIN = 1 ! Sundquist
  # IWARMRAIN = 2 ! Khairoutdinov and Kogan (2000)
  # ---------------------------------------------------------------------
  iwarmrain = 2
  # ---------------------------------------------------------------------
  # Set version of rain evaporation
  # IEVAPRAIN = 1 ! Sundquist
  # IEVAPRAIN = 2 ! Abel and Boutle (2013)
  # ---------------------------------------------------------------------
  ievaprain = 2
  # ---------------------------------------------------------------------
  # Set version of snow evaporation
  # IEVAPSNOW = 1 ! Sundquist
  # IEVAPSNOW = 2 ! New
  # ---------------------------------------------------------------------
  ievapsnow = 1
  # ---------------------------------------------------------------------
  # Set version of ice deposition
  # IDEPICE = 1 ! Rotstayn (2001)
  # IDEPICE = 2 ! New
  # ---------------------------------------------------------------------
  idepice = 1
  
  # ---------------------
  # Some simple constants
  # ---------------------
  zqtmst = 1.0 / ptsphy
  zgdcp = ydcst.rg / ydcst.rcpd
  zrdcp = ydcst.rd / ydcst.rcpd
  zcons1a = ydcst.rcpd / (ydcst.rlmlt*ydcst.rg*yrecldp.rtaumel)
  zepsec = 1.E-14
  zrg_r = 1.0 / ydcst.rg
  zrldcp = 1.0 / (ydthf.ralsdcp - ydthf.ralvdcp)
  
  # Note: Defined in module/yoecldp.F90
  # NCLDQL=1    ! liquid cloud water
  # NCLDQI=2    ! ice cloud water
  # NCLDQR=3    ! rain water
  # NCLDQS=4    ! snow
  # NCLDQV=5    ! vapour
  
  # -----------------------------------------------
  # Define species phase, 0=vapour, 1=liquid, 2=ice
  # -----------------------------------------------
  iphase[ncldqv - 1] = 0
  iphase[ncldql - 1] = 1
  iphase[ncldqr - 1] = 1
  iphase[ncldqi - 1] = 2
  iphase[ncldqs - 1] = 2
  
  # ---------------------------------------------------
  # Set up melting/freezing index,
  # if an ice category melts/freezes, where does it go?
  # ---------------------------------------------------
  imelt[ncldqv - 1] = -99
  imelt[ncldql - 1] = ncldqi
  imelt[ncldqr - 1] = ncldqs
  imelt[ncldqi - 1] = ncldqr
  imelt[ncldqs - 1] = ncldqr
  
  # -----------------------------------------------
  # INITIALIZATION OF OUTPUT TENDENCIES
  # -----------------------------------------------
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      tendency_loc_t[jk - 1, jl - 1] = 0.0
      tendency_loc_q[jk - 1, jl - 1] = 0.0
      tendency_loc_a[jk - 1, jl - 1] = 0.0
  for jm in range(1, nclv - 1 + 1):
    for jk in range(1, klev + 1):
      for jl in range(kidia, kfdia + 1):
        tendency_loc_cld[jm - 1, jk - 1, jl - 1] = 0.0
  
  #-- These were uninitialized : meaningful only when we compare error differences
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      pcovptot[jk - 1, jl - 1] = 0.0
      tendency_loc_cld[nclv - 1, jk - 1, jl - 1] = 0.0
  
  # -------------------------
  # set up fall speeds in m/s
  # -------------------------
  zvqx[ncldqv - 1] = 0.0
  zvqx[ncldql - 1] = 0.0
  zvqx[ncldqi - 1] = yrecldp.rvice
  zvqx[ncldqr - 1] = yrecldp.rvrain
  zvqx[ncldqs - 1] = yrecldp.rvsnow
  llfall[:] = False
  for jm in range(1, nclv + 1):
    if zvqx[jm - 1] > 0.0:
      llfall[jm - 1] = True
    # falling species
  # Set LLFALL to false for ice (but ice still sediments!)
  # Need to rationalise this at some point
  llfall[ncldqi - 1] = False
  
  
  #######################################################################
  #             1.  *** INITIAL VALUES FOR VARIABLES ***
  #######################################################################
  
  
  # ----------------------
  # non CLV initialization
  # ----------------------
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      ztp1[jk - 1, jl - 1] = pt[jk - 1, jl - 1] + ptsphy*tendency_tmp_t[jk - 1, jl - 1]
      zqx[ncldqv - 1, jk - 1, jl - 1] = pq[jk - 1, jl - 1] + ptsphy*tendency_tmp_q[jk - 1, jl - 1]
      zqx0[ncldqv - 1, jk - 1, jl - 1] = pq[jk - 1, jl - 1] + ptsphy*tendency_tmp_q[jk - 1, jl - 1]
      za[jk - 1, jl - 1] = pa[jk - 1, jl - 1] + ptsphy*tendency_tmp_a[jk - 1, jl - 1]
      zaorig[jk - 1, jl - 1] = pa[jk - 1, jl - 1] + ptsphy*tendency_tmp_a[jk - 1, jl - 1]
  
  # -------------------------------------
  # initialization for CLV family
  # -------------------------------------
  for jm in range(1, nclv - 1 + 1):
    for jk in range(1, klev + 1):
      for jl in range(kidia, kfdia + 1):
        zqx[jm - 1, jk - 1, jl - 1] = pclv[jm - 1, jk - 1, jl - 1] + ptsphy*tendency_tmp_cld[jm - 1, jk - 1, jl - 1]
        zqx0[jm - 1, jk - 1, jl - 1] = pclv[jm - 1, jk - 1, jl - 1] + ptsphy*tendency_tmp_cld[jm - 1, jk - 1, jl - 1]
  
  #-------------
  # zero arrays
  #-------------
  for jm in range(1, nclv + 1):
    for jk in range(1, klev + 1 + 1):
      for jl in range(kidia, kfdia + 1):
        zpfplsx[jm - 1, jk - 1, jl - 1] = 0.0          # precip fluxes
  
  for jm in range(1, nclv + 1):
    for jk in range(1, klev + 1):
      for jl in range(kidia, kfdia + 1):
        zqxn2d[jm - 1, jk - 1, jl - 1] = 0.0          # end of timestep values in 2D
        zlneg[jm - 1, jk - 1, jl - 1] = 0.0          # negative input check
  
  for jl in range(kidia, kfdia + 1):
    prainfrac_toprfz[jl - 1] = 0.0      # rain fraction at top of refreezing layer
  llrainliq[:] = True    # Assume all raindrops are liquid initially
  
  # ----------------------------------------------------
  # Tidy up very small cloud cover or total cloud water
  # ----------------------------------------------------
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      if zqx[ncldql - 1, jk - 1, jl - 1] + zqx[ncldqi - 1, jk - 1, jl - 1] < yrecldp.rlmin or za[jk - 1, jl - 1] < yrecldp.ramin:
        
        # Evaporate small cloud liquid water amounts
        zlneg[ncldql - 1, jk - 1, jl - 1] = zlneg[ncldql - 1, jk - 1, jl - 1] + zqx[ncldql - 1, jk - 1, jl - 1]
        zqadj = zqx[ncldql - 1, jk - 1, jl - 1]*zqtmst
        tendency_loc_q[jk - 1, jl - 1] = tendency_loc_q[jk - 1, jl - 1] + zqadj
        tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] - ydthf.ralvdcp*zqadj
        zqx[ncldqv - 1, jk - 1, jl - 1] = zqx[ncldqv - 1, jk - 1, jl - 1] + zqx[ncldql - 1, jk - 1, jl - 1]
        zqx[ncldql - 1, jk - 1, jl - 1] = 0.0
        
        # Evaporate small cloud ice water amounts
        zlneg[ncldqi - 1, jk - 1, jl - 1] = zlneg[ncldqi - 1, jk - 1, jl - 1] + zqx[ncldqi - 1, jk - 1, jl - 1]
        zqadj = zqx[ncldqi - 1, jk - 1, jl - 1]*zqtmst
        tendency_loc_q[jk - 1, jl - 1] = tendency_loc_q[jk - 1, jl - 1] + zqadj
        tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] - ydthf.ralsdcp*zqadj
        zqx[ncldqv - 1, jk - 1, jl - 1] = zqx[ncldqv - 1, jk - 1, jl - 1] + zqx[ncldqi - 1, jk - 1, jl - 1]
        zqx[ncldqi - 1, jk - 1, jl - 1] = 0.0
        
        # Set cloud cover to zero
        za[jk - 1, jl - 1] = 0.0
        
  
  # ---------------------------------
  # Tidy up small CLV variables
  # ---------------------------------
  #DIR$ IVDEP
  for jm in range(1, nclv - 1 + 1):
    #DIR$ IVDEP
    for jk in range(1, klev + 1):
      #DIR$ IVDEP
      for jl in range(kidia, kfdia + 1):
        if zqx[jm - 1, jk - 1, jl - 1] < yrecldp.rlmin:
          zlneg[jm - 1, jk - 1, jl - 1] = zlneg[jm - 1, jk - 1, jl - 1] + zqx[jm - 1, jk - 1, jl - 1]
          zqadj = zqx[jm - 1, jk - 1, jl - 1]*zqtmst
          tendency_loc_q[jk - 1, jl - 1] = tendency_loc_q[jk - 1, jl - 1] + zqadj
          if iphase[jm - 1] == 1:
            tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] - ydthf.ralvdcp*zqadj
          if iphase[jm - 1] == 2:
            tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] - ydthf.ralsdcp*zqadj
          zqx[ncldqv - 1, jk - 1, jl - 1] = zqx[ncldqv - 1, jk - 1, jl - 1] + zqx[jm - 1, jk - 1, jl - 1]
          zqx[jm - 1, jk - 1, jl - 1] = 0.0
  
  
  # ------------------------------
  # Define saturation values
  # ------------------------------
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      #----------------------------------------
      # old *diagnostic* mixed phase saturation
      #----------------------------------------
      zfoealfa[jk - 1, jl - 1] = foealfa(ztp1[jk - 1, jl - 1])
      zfoeewmt[jk - 1, jl - 1] = min(foeewm(ztp1[jk - 1, jl - 1]) / pap[jk - 1, jl - 1], 0.5)
      zqsmix[jk - 1, jl - 1] = zfoeewmt[jk - 1, jl - 1]
      zqsmix[jk - 1, jl - 1] = zqsmix[jk - 1, jl - 1] / (1.0 - ydcst.retv*zqsmix[jk - 1, jl - 1])
      
      #---------------------------------------------
      # ice saturation T<273K
      # liquid water saturation for T>273K
      #---------------------------------------------
      zalfa = foedelta(ztp1[jk - 1, jl - 1])
      zfoeew[jk - 1, jl - 1] = min((zalfa*foeeliq(ztp1[jk - 1, jl - 1]) + (1.0 - zalfa)*foeeice(ztp1[jk - 1, jl - 1])) / pap[jk - 1, jl - 1], 0.5)
      zfoeew[jk - 1, jl - 1] = min(0.5, zfoeew[jk - 1, jl - 1])
      zqsice[jk - 1, jl - 1] = zfoeew[jk - 1, jl - 1] / (1.0 - ydcst.retv*zfoeew[jk - 1, jl - 1])
      
      #----------------------------------
      # liquid water saturation
      #----------------------------------
      zfoeeliqt[jk - 1, jl - 1] = min(foeeliq(ztp1[jk - 1, jl - 1]) / pap[jk - 1, jl - 1], 0.5)
      zqsliq[jk - 1, jl - 1] = zfoeeliqt[jk - 1, jl - 1]
      zqsliq[jk - 1, jl - 1] = zqsliq[jk - 1, jl - 1] / (1.0 - ydcst.retv*zqsliq[jk - 1, jl - 1])
      
      #   !----------------------------------
      #   ! ice water saturation
      #   !----------------------------------
      #   ZFOEEICET(JL,JK)=MIN(FOEEICE(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      #   ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
      #   ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-RETV*ZQSICE(JL,JK))
    
  
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      
      
      #------------------------------------------
      # Ensure cloud fraction is between 0 and 1
      #------------------------------------------
      za[jk - 1, jl - 1] = max(0.0, min(1.0, za[jk - 1, jl - 1]))
      
      #-------------------------------------------------------------------
      # Calculate liq/ice fractions (no longer a diagnostic relationship)
      #-------------------------------------------------------------------
      zli[jk - 1, jl - 1] = zqx[ncldql - 1, jk - 1, jl - 1] + zqx[ncldqi - 1, jk - 1, jl - 1]
      if zli[jk - 1, jl - 1] > yrecldp.rlmin:
        zliqfrac[jk - 1, jl - 1] = zqx[ncldql - 1, jk - 1, jl - 1] / zli[jk - 1, jl - 1]
        zicefrac[jk - 1, jl - 1] = 1.0 - zliqfrac[jk - 1, jl - 1]
      else:
        zliqfrac[jk - 1, jl - 1] = 0.0
        zicefrac[jk - 1, jl - 1] = 0.0
      
  
  #######################################################################
  #        2.       *** CONSTANTS AND PARAMETERS ***
  #######################################################################
  #  Calculate L in updrafts of bl-clouds
  #  Specify QS, P/PS for tropopause (for c2)
  #  And initialize variables
  #------------------------------------------
  
  #---------------------------------
  # Find tropopause level (ZTRPAUS)
  #---------------------------------
  for jl in range(kidia, kfdia + 1):
    ztrpaus[jl - 1] = 0.1
    zpaphd[jl - 1] = 1.0 / paph[klev + 1 - 1, jl - 1]
  for jk in range(1, klev - 1 + 1):
    for jl in range(kidia, kfdia + 1):
      zsig = pap[jk - 1, jl - 1]*zpaphd[jl - 1]
      if zsig > 0.1 and zsig < 0.4 and ztp1[jk - 1, jl - 1] > ztp1[jk + 1 - 1, jl - 1]:
        ztrpaus[jl - 1] = zsig
  
  #-----------------------------
  # Reset single level variables
  #-----------------------------
  
  for jl in range(kidia, kfdia + 1):
    zanewm1[jl - 1] = 0.0
    zda[jl - 1] = 0.0
    zcovpclr[jl - 1] = 0.0
    zcovpmax[jl - 1] = 0.0
    zcovptot[jl - 1] = 0.0
    zcldtopdist[jl - 1] = 0.0
  
  #######################################################################
  #           3.       *** PHYSICS ***
  #######################################################################
  
  
  #----------------------------------------------------------------------
  #                       START OF VERTICAL LOOP
  #----------------------------------------------------------------------
  
  for jk in range(yrecldp.ncldtop, klev + 1):
    
    #----------------------------------------------------------------------
    # 3.0 INITIALIZE VARIABLES
    #----------------------------------------------------------------------
    
    #---------------------------------
    # First guess microphysics
    #---------------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zqxfg[jm - 1, jl - 1] = zqx[jm - 1, jk - 1, jl - 1]
    
    #---------------------------------
    # Set KLON arrays to zero
    #---------------------------------
    
    for jl in range(kidia, kfdia + 1):
      zlicld[jl - 1] = 0.0
      zrainaut[jl - 1] = 0.0        # currently needed for diags
      zrainacc[jl - 1] = 0.0        # currently needed for diags
      zsnowaut[jl - 1] = 0.0        # needed
      zldefr[jl - 1] = 0.0
      zacust[jl - 1] = 0.0        # set later when needed
      zqpretot[jl - 1] = 0.0
      zlfinalsum[jl - 1] = 0.0
      
      # Required for first guess call
      zlcond1[jl - 1] = 0.0
      zlcond2[jl - 1] = 0.0
      zsupsat[jl - 1] = 0.0
      zlevapl[jl - 1] = 0.0
      zlevapi[jl - 1] = 0.0
      
      #-------------------------------------
      # solvers for cloud fraction
      #-------------------------------------
      zsolab[jl - 1] = 0.0
      zsolac[jl - 1] = 0.0
      
      zicetot[jl - 1] = 0.0
    
    #------------------------------------------
    # reset matrix so missing pathways are set
    #------------------------------------------
    for jm in range(1, nclv + 1):
      for jn in range(1, nclv + 1):
        for jl in range(kidia, kfdia + 1):
          zsolqb[jm - 1, jn - 1, jl - 1] = 0.0
          zsolqa[jm - 1, jn - 1, jl - 1] = 0.0
    
    #----------------------------------
    # reset new microphysics variables
    #----------------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zfallsrce[jm - 1, jl - 1] = 0.0
        zfallsink[jm - 1, jl - 1] = 0.0
        zconvsrce[jm - 1, jl - 1] = 0.0
        zconvsink[jm - 1, jl - 1] = 0.0
        zpsupsatsrce[jm - 1, jl - 1] = 0.0
        zratio[jm - 1, jl - 1] = 0.0
    
    for jl in range(kidia, kfdia + 1):
      
      #-------------------------
      # derived variables needed
      #-------------------------
      
      zdp[jl - 1] = paph[jk + 1 - 1, jl - 1] - paph[jk - 1, jl - 1]        # dp
      zgdp[jl - 1] = ydcst.rg / zdp[jl - 1]        # g/dp
      zrho[jl - 1] = pap[jk - 1, jl - 1] / (ydcst.rd*ztp1[jk - 1, jl - 1])        # p/RT air density
      
      zdtgdp[jl - 1] = ptsphy*zgdp[jl - 1]        # dt g/dp
      zrdtgdp[jl - 1] = zdp[jl - 1]*(1.0 / (ptsphy*ydcst.rg))        # 1/(dt g/dp)
      
      if jk > 1:
        zdtgdpf[jl - 1] = ptsphy*ydcst.rg / (pap[jk - 1, jl - 1] - pap[jk - 1 - 1, jl - 1])
      
      #------------------------------------
      # Calculate dqs/dT correction factor
      #------------------------------------
      # Reminder: RETV=RV/RD-1
      
      # liquid
      zfacw = ydthf.r5les / ((ztp1[jk - 1, jl - 1] - ydthf.r4les)**2)
      zcor = 1.0 / (1.0 - ydcst.retv*zfoeeliqt[jk - 1, jl - 1])
      zdqsliqdt[jl - 1] = zfacw*zcor*zqsliq[jk - 1, jl - 1]
      zcorqsliq[jl - 1] = 1.0 + ydthf.ralvdcp*zdqsliqdt[jl - 1]
      
      # ice
      zfaci = ydthf.r5ies / ((ztp1[jk - 1, jl - 1] - ydthf.r4ies)**2)
      zcor = 1.0 / (1.0 - ydcst.retv*zfoeew[jk - 1, jl - 1])
      zdqsicedt[jl - 1] = zfaci*zcor*zqsice[jk - 1, jl - 1]
      zcorqsice[jl - 1] = 1.0 + ydthf.ralsdcp*zdqsicedt[jl - 1]
      
      # diagnostic mixed
      zalfaw = zfoealfa[jk - 1, jl - 1]
      zalfawm[jl - 1] = zalfaw
      zfac = zalfaw*zfacw + (1.0 - zalfaw)*zfaci
      zcor = 1.0 / (1.0 - ydcst.retv*zfoeewmt[jk - 1, jl - 1])
      zdqsmixdt[jl - 1] = zfac*zcor*zqsmix[jk - 1, jl - 1]
      zcorqsmix[jl - 1] = 1.0 + foeldcpm(ztp1[jk - 1, jl - 1])*zdqsmixdt[jl - 1]
      
      # evaporation/sublimation limits
      zevaplimmix[jl - 1] = max((zqsmix[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1]) / zcorqsmix[jl - 1], 0.0)
      zevaplimliq[jl - 1] = max((zqsliq[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1]) / zcorqsliq[jl - 1], 0.0)
      zevaplimice[jl - 1] = max((zqsice[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1]) / zcorqsice[jl - 1], 0.0)
      
      #--------------------------------
      # in-cloud consensate amount
      #--------------------------------
      ztmpa = 1.0 / max(za[jk - 1, jl - 1], zepsec)
      zliqcld[jl - 1] = zqx[ncldql - 1, jk - 1, jl - 1]*ztmpa
      zicecld[jl - 1] = zqx[ncldqi - 1, jk - 1, jl - 1]*ztmpa
      zlicld[jl - 1] = zliqcld[jl - 1] + zicecld[jl - 1]
      
    
    #------------------------------------------------
    # Evaporate very small amounts of liquid and ice
    #------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      
      if zqx[ncldql - 1, jk - 1, jl - 1] < yrecldp.rlmin:
        zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zqx[ncldql - 1, jk - 1, jl - 1]
        zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = -zqx[ncldql - 1, jk - 1, jl - 1]
      
      if zqx[ncldqi - 1, jk - 1, jl - 1] < yrecldp.rlmin:
        zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zqx[ncldqi - 1, jk - 1, jl - 1]
        zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = -zqx[ncldqi - 1, jk - 1, jl - 1]
      
    
    #---------------------------------------------------------------------
    #  3.1  ICE SUPERSATURATION ADJUSTMENT
    #---------------------------------------------------------------------
    # Note that the supersaturation adjustment is made with respect to
    # liquid saturation:  when T>0C
    # ice saturation:     when T<0C
    #                     with an adjustment made to allow for ice
    #                     supersaturation in the clear sky
    # Note also that the KOOP factor automatically clips the supersaturation
    # to a maximum set by the liquid water saturation mixing ratio
    # important for temperatures near to but below 0C
    #-----------------------------------------------------------------------
    
    #DIR$ NOFUSION
    for jl in range(kidia, kfdia + 1):
      
      #-----------------------------------
      # 3.1.1 Supersaturation limit (from Koop)
      #-----------------------------------
      # Needs to be set for all temperatures
      zfokoop[jl - 1] = fokoop(ztp1[jk - 1, jl - 1])
    for jl in range(kidia, kfdia + 1):
      
      if ztp1[jk - 1, jl - 1] >= ydcst.rtt or yrecldp.nssopt == 0:
        zfac = 1.0
        zfaci = 1.0
      else:
        zfac = za[jk - 1, jl - 1] + zfokoop[jl - 1]*(1.0 - za[jk - 1, jl - 1])
        zfaci = ptsphy / yrecldp.rkooptau
      
      #-------------------------------------------------------------------
      # 3.1.2 Calculate supersaturation wrt Koop including dqs/dT
      #       correction factor
      # [#Note: QSICE or QSLIQ]
      #-------------------------------------------------------------------
      
      # Calculate supersaturation to add to cloud
      if za[jk - 1, jl - 1] > 1.0 - yrecldp.ramin:
        zsupsat[jl - 1] = max((zqx[ncldqv - 1, jk - 1, jl - 1] - zfac*zqsice[jk - 1, jl - 1]) / zcorqsice[jl - 1], 0.0)
      else:
        # Calculate environmental humidity supersaturation
        zqp1env = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsice[jk - 1, jl - 1]) / max(1.0 - za[jk - 1, jl - 1], zepsilon)
        #& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
        zsupsat[jl - 1] = max((1.0 - za[jk - 1, jl - 1])*(zqp1env - zfac*zqsice[jk - 1, jl - 1]) / zcorqsice[jl - 1], 0.0)
      
      #-------------------------------------------------------------------
      # Here the supersaturation is turned into liquid water
      # However, if the temperature is below the threshold for homogeneous
      # freezing then the supersaturation is turned instantly to ice.
      #--------------------------------------------------------------------
      
      if zsupsat[jl - 1] > zepsec:
        
        if ztp1[jk - 1, jl - 1] > yrecldp.rthomo:
          # Turn supersaturation into liquid water
          zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = zsolqa[ncldqv - 1, ncldql - 1, jl - 1] + zsupsat[jl - 1]
          zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zsolqa[ncldql - 1, ncldqv - 1, jl - 1] - zsupsat[jl - 1]
          # Include liquid in first guess
          zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] + zsupsat[jl - 1]
        else:
          # Turn supersaturation into ice water
          zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] + zsupsat[jl - 1]
          zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] - zsupsat[jl - 1]
          # Add ice to first guess for deposition term
          zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zsupsat[jl - 1]
        
        # Increase cloud amount using RKOOPTAU timescale
        zsolac[jl - 1] = (1.0 - za[jk - 1, jl - 1])*zfaci
        
      
      #-------------------------------------------------------
      # 3.1.3 Include supersaturation from previous timestep
      # (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
      #-------------------------------------------------------
      if psupsat[jk - 1, jl - 1] > zepsec:
        if ztp1[jk - 1, jl - 1] > yrecldp.rthomo:
          # Turn supersaturation into liquid water
          zsolqa[ncldql - 1, ncldql - 1, jl - 1] = zsolqa[ncldql - 1, ncldql - 1, jl - 1] + psupsat[jk - 1, jl - 1]
          zpsupsatsrce[ncldql - 1, jl - 1] = psupsat[jk - 1, jl - 1]
          # Add liquid to first guess for deposition term
          zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] + psupsat[jk - 1, jl - 1]
          # Store cloud budget diagnostics if required
        else:
          # Turn supersaturation into ice water
          zsolqa[ncldqi - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqi - 1, jl - 1] + psupsat[jk - 1, jl - 1]
          zpsupsatsrce[ncldqi - 1, jl - 1] = psupsat[jk - 1, jl - 1]
          # Add ice to first guess for deposition term
          zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + psupsat[jk - 1, jl - 1]
          # Store cloud budget diagnostics if required
        
        # Increase cloud amount using RKOOPTAU timescale
        zsolac[jl - 1] = (1.0 - za[jk - 1, jl - 1])*zfaci
        # Store cloud budget diagnostics if required
      
    # on JL
    
    #---------------------------------------------------------------------
    #  3.2  DETRAINMENT FROM CONVECTION
    #---------------------------------------------------------------------
    # * Diagnostic T-ice/liq split retained for convection
    #    Note: This link is now flexible and a future convection
    #    scheme can detrain explicit seperate budgets of:
    #    cloud water, ice, rain and snow
    # * There is no (1-ZA) multiplier term on the cloud detrainment
    #    term, since is now written in mass-flux terms
    # [#Note: Should use ZFOEALFACU used in convection rather than ZFOEALFA]
    #---------------------------------------------------------------------
    if jk < klev and jk >= yrecldp.ncldtop:
      
      for jl in range(kidia, kfdia + 1):
        
        plude[jk - 1, jl - 1] = plude[jk - 1, jl - 1]*zdtgdp[jl - 1]
        
        if ldcum[jl - 1] and plude[jk - 1, jl - 1] > yrecldp.rlmin and plu[jk + 1 - 1, jl - 1] > zepsec:
          
          zsolac[jl - 1] = zsolac[jl - 1] + plude[jk - 1, jl - 1] / plu[jk + 1 - 1, jl - 1]
          # *diagnostic temperature split*
          zalfaw = zfoealfa[jk - 1, jl - 1]
          zconvsrce[ncldql - 1, jl - 1] = zalfaw*plude[jk - 1, jl - 1]
          zconvsrce[ncldqi - 1, jl - 1] = (1.0 - zalfaw)*plude[jk - 1, jl - 1]
          zsolqa[ncldql - 1, ncldql - 1, jl - 1] = zsolqa[ncldql - 1, ncldql - 1, jl - 1] + zconvsrce[ncldql - 1, jl - 1]
          zsolqa[ncldqi - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqi - 1, jl - 1] + zconvsrce[ncldqi - 1, jl - 1]
          
        else:
          
          plude[jk - 1, jl - 1] = 0.0
          
        # *convective snow detrainment source
        if ldcum[jl - 1]:
          zsolqa[ncldqs - 1, ncldqs - 1, jl - 1] = zsolqa[ncldqs - 1, ncldqs - 1, jl - 1] + psnde[jk - 1, jl - 1]*zdtgdp[jl - 1]
        
      
    # JK<KLEV
    
    #---------------------------------------------------------------------
    #  3.3  SUBSIDENCE COMPENSATING CONVECTIVE UPDRAUGHTS
    #---------------------------------------------------------------------
    # Three terms:
    # * Convective subsidence source of cloud from layer above
    # * Evaporation of cloud within the layer
    # * Subsidence sink of cloud to the layer below (Implicit solution)
    #---------------------------------------------------------------------
    
    #-----------------------------------------------
    # Subsidence source from layer above
    #               and
    # Evaporation of cloud within the layer
    #-----------------------------------------------
    if jk > yrecldp.ncldtop:
      
      for jl in range(kidia, kfdia + 1):
        zmf[jl - 1] = max(0.0, (pmfu[jk - 1, jl - 1] + pmfd[jk - 1, jl - 1])*zdtgdp[jl - 1])
        zacust[jl - 1] = zmf[jl - 1]*zanewm1[jl - 1]
      
      for jm in range(1, nclv + 1):
        if not llfall[jm - 1] and iphase[jm - 1] > 0:
          for jl in range(kidia, kfdia + 1):
            zlcust[jm - 1, jl - 1] = zmf[jl - 1]*zqxnm1[jm - 1, jl - 1]
            # record total flux for enthalpy budget:
            zconvsrce[jm - 1, jl - 1] = zconvsrce[jm - 1, jl - 1] + zlcust[jm - 1, jl - 1]
      
      # Now have to work out how much liquid evaporates at arrival point
      # since there is no prognostic memory for in-cloud humidity, i.e.
      # we always assume cloud is saturated.
      
      for jl in range(kidia, kfdia + 1):
        zdtdp = zrdcp*0.5*(ztp1[jk - 1 - 1, jl - 1] + ztp1[jk - 1, jl - 1]) / paph[jk - 1, jl - 1]
        zdtforc = zdtdp*(pap[jk - 1, jl - 1] - pap[jk - 1 - 1, jl - 1])
        #[#Note: Diagnostic mixed phase should be replaced below]
        zdqs[jl - 1] = zanewm1[jl - 1]*zdtforc*zdqsmixdt[jl - 1]
      
      for jm in range(1, nclv + 1):
        if not llfall[jm - 1] and iphase[jm - 1] > 0:
          for jl in range(kidia, kfdia + 1):
            zlfinal = max(0.0, zlcust[jm - 1, jl - 1] - zdqs[jl - 1])              #lim to zero
            # no supersaturation allowed incloud ---V
            zevap = min((zlcust[jm - 1, jl - 1] - zlfinal), zevaplimmix[jl - 1])
            #          ZEVAP=0.0_JPRB
            zlfinal = zlcust[jm - 1, jl - 1] - zevap
            zlfinalsum[jl - 1] = zlfinalsum[jl - 1] + zlfinal              # sum
            
            zsolqa[jm - 1, jm - 1, jl - 1] = zsolqa[jm - 1, jm - 1, jl - 1] + zlcust[jm - 1, jl - 1]              # whole sum
            zsolqa[jm - 1, ncldqv - 1, jl - 1] = zsolqa[jm - 1, ncldqv - 1, jl - 1] + zevap
            zsolqa[ncldqv - 1, jm - 1, jl - 1] = zsolqa[ncldqv - 1, jm - 1, jl - 1] - zevap
      
      #  Reset the cloud contribution if no cloud water survives to this level:
      for jl in range(kidia, kfdia + 1):
        if zlfinalsum[jl - 1] < zepsec:
          zacust[jl - 1] = 0.0
        zsolac[jl - 1] = zsolac[jl - 1] + zacust[jl - 1]
      
    # on  JK>NCLDTOP
    
    #---------------------------------------------------------------------
    # Subsidence sink of cloud to the layer below
    # (Implicit - re. CFL limit on convective mass flux)
    #---------------------------------------------------------------------
    
    for jl in range(kidia, kfdia + 1):
      
      if jk < klev:
        
        zmfdn = max(0.0, (pmfu[jk + 1 - 1, jl - 1] + pmfd[jk + 1 - 1, jl - 1])*zdtgdp[jl - 1])
        
        zsolab[jl - 1] = zsolab[jl - 1] + zmfdn
        zsolqb[ncldql - 1, ncldql - 1, jl - 1] = zsolqb[ncldql - 1, ncldql - 1, jl - 1] + zmfdn
        zsolqb[ncldqi - 1, ncldqi - 1, jl - 1] = zsolqb[ncldqi - 1, ncldqi - 1, jl - 1] + zmfdn
        
        # Record sink for cloud budget and enthalpy budget diagnostics
        zconvsink[ncldql - 1, jl - 1] = zmfdn
        zconvsink[ncldqi - 1, jl - 1] = zmfdn
        
      
    
    #----------------------------------------------------------------------
    # 3.4  EROSION OF CLOUDS BY TURBULENT MIXING
    #----------------------------------------------------------------------
    # NOTE: In default tiedtke scheme this process decreases the cloud
    #       area but leaves the specific cloud water content
    #       within clouds unchanged
    #----------------------------------------------------------------------
    
    # ------------------------------
    # Define turbulent erosion rate
    # ------------------------------
    for jl in range(kidia, kfdia + 1):
      zldifdt[jl - 1] = yrecldp.rcldiff*ptsphy        #original version
      #Increase by factor of 5 for convective points
      if ktype[jl - 1] > 0 and plude[jk - 1, jl - 1] > zepsec:
        zldifdt[jl - 1] = yrecldp.rcldiff_convi*zldifdt[jl - 1]
    
    # At the moment, works on mixed RH profile and partitioned ice/liq fraction
    # so that it is similar to previous scheme
    # Should apply RHw for liquid cloud and RHi for ice cloud separately
    for jl in range(kidia, kfdia + 1):
      if zli[jk - 1, jl - 1] > zepsec:
        # Calculate environmental humidity
        #      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
        #    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
        #      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
        ze = zldifdt[jl - 1]*max(zqsmix[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1], 0.0)
        zleros = za[jk - 1, jl - 1]*ze
        zleros = min(zleros, zevaplimmix[jl - 1])
        zleros = min(zleros, zli[jk - 1, jl - 1])
        zaeros = zleros / zlicld[jl - 1]          #if linear term
        
        # Erosion is -ve LINEAR in L,A
        zsolac[jl - 1] = zsolac[jl - 1] - zaeros          #linear
        
        zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zsolqa[ncldql - 1, ncldqv - 1, jl - 1] + zliqfrac[jk - 1, jl - 1]*zleros
        zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = zsolqa[ncldqv - 1, ncldql - 1, jl - 1] - zliqfrac[jk - 1, jl - 1]*zleros
        zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] + zicefrac[jk - 1, jl - 1]*zleros
        zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] - zicefrac[jk - 1, jl - 1]*zleros
        
    
    #----------------------------------------------------------------------
    # 3.4  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
    #----------------------------------------------------------------------
    #  calculate dqs/dt
    #  Note: For the separate prognostic Qi and Ql, one would ideally use
    #  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
    #  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
    #  These would then instantaneous freeze if T<-38C or lead to ice growth
    #  by deposition in warmer mixed phase clouds.  However, since we do
    #  not have a separate prognostic equation for in-cloud humidity or a
    #  statistical scheme approach in place, the depositional growth of ice
    #  in the mixed phase can not be modelled and we resort to supersaturation
    #  wrt ice instanteously converting to ice over one timestep
    #  (see Tompkins et al. QJRMS 2007 for details)
    #  Thus for the initial implementation the diagnostic mixed phase is
    #  retained for the moment, and the level of approximation noted.
    #----------------------------------------------------------------------
    
    for jl in range(kidia, kfdia + 1):
      zdtdp = zrdcp*ztp1[jk - 1, jl - 1] / pap[jk - 1, jl - 1]
      zdpmxdt = zdp[jl - 1]*zqtmst
      zmfdn = 0.0
      if jk < klev:
        zmfdn = pmfu[jk + 1 - 1, jl - 1] + pmfd[jk + 1 - 1, jl - 1]
      zwtot = pvervel[jk - 1, jl - 1] + 0.5*ydcst.rg*(pmfu[jk - 1, jl - 1] + pmfd[jk - 1, jl - 1] + zmfdn)
      zwtot = min(zdpmxdt, max(-zdpmxdt, zwtot))
      zzzdt = phrsw[jk - 1, jl - 1] + phrlw[jk - 1, jl - 1]
      zdtdiab = min(zdpmxdt*zdtdp, max(-zdpmxdt*zdtdp, zzzdt))*ptsphy + ydthf.ralfdcp*zldefr[jl - 1]
      # Note: ZLDEFR should be set to the difference between the mixed phase functions
      # in the convection and cloud scheme, but this is not calculated, so is zero and
      # the functions must be the same
      zdtforc = zdtdp*zwtot*ptsphy + zdtdiab
      zqold[jl - 1] = zqsmix[jk - 1, jl - 1]
      ztold[jl - 1] = ztp1[jk - 1, jl - 1]
      ztp1[jk - 1, jl - 1] = ztp1[jk - 1, jl - 1] + zdtforc
      ztp1[jk - 1, jl - 1] = max(ztp1[jk - 1, jl - 1], 160.0)
      llflag[jl - 1] = True
    
    # Formerly a call to CUADJTQ(..., ICALL=5)
    for jl in range(kidia, kfdia + 1):
      zqp = 1.0 / pap[jk - 1, jl - 1]
      zqsat = foeewm(ztp1[jk - 1, jl - 1])*zqp
      zqsat = min(0.5, zqsat)
      zcor = 1.0 / (1.0 - ydcst.retv*zqsat)
      zqsat = zqsat*zcor
      zcond = (zqsmix[jk - 1, jl - 1] - zqsat) / (1.0 + zqsat*zcor*foedem(ztp1[jk - 1, jl - 1]))
      ztp1[jk - 1, jl - 1] = ztp1[jk - 1, jl - 1] + foeldcpm(ztp1[jk - 1, jl - 1])*zcond
      zqsmix[jk - 1, jl - 1] = zqsmix[jk - 1, jl - 1] - zcond
      zqsat = foeewm(ztp1[jk - 1, jl - 1])*zqp
      zqsat = min(0.5, zqsat)
      zcor = 1.0 / (1.0 - ydcst.retv*zqsat)
      zqsat = zqsat*zcor
      zcond1 = (zqsmix[jk - 1, jl - 1] - zqsat) / (1.0 + zqsat*zcor*foedem(ztp1[jk - 1, jl - 1]))
      ztp1[jk - 1, jl - 1] = ztp1[jk - 1, jl - 1] + foeldcpm(ztp1[jk - 1, jl - 1])*zcond1
      zqsmix[jk - 1, jl - 1] = zqsmix[jk - 1, jl - 1] - zcond1
    
    for jl in range(kidia, kfdia + 1):
      zdqs[jl - 1] = zqsmix[jk - 1, jl - 1] - zqold[jl - 1]
      zqsmix[jk - 1, jl - 1] = zqold[jl - 1]
      ztp1[jk - 1, jl - 1] = ztold[jl - 1]
    
    #----------------------------------------------------------------------
    # 3.4a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
    # ----------------------------------------------------------------------
    # Erosion term is LINEAR in L
    # Changed to be uniform distribution in cloud region
    
    for jl in range(kidia, kfdia + 1):
      
      # Previous function based on DELTA DISTRIBUTION in cloud:
      if zdqs[jl - 1] > 0.0:
        #    If subsidence evaporation term is turned off, then need to use updated
        #    liquid and cloud here?
        #    ZLEVAP = MAX(ZA(JL,JK)+ZACUST(JL),1.0_JPRB)*MIN(ZDQS(JL),ZLICLD(JL)+ZLFINALSUM(JL))
        zlevap = za[jk - 1, jl - 1]*min(zdqs[jl - 1], zlicld[jl - 1])
        zlevap = min(zlevap, zevaplimmix[jl - 1])
        zlevap = min(zlevap, max(zqsmix[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1], 0.0))
        
        # For first guess call
        zlevapl[jl - 1] = zliqfrac[jk - 1, jl - 1]*zlevap
        zlevapi[jl - 1] = zicefrac[jk - 1, jl - 1]*zlevap
        
        zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zsolqa[ncldql - 1, ncldqv - 1, jl - 1] + zliqfrac[jk - 1, jl - 1]*zlevap
        zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = zsolqa[ncldqv - 1, ncldql - 1, jl - 1] - zliqfrac[jk - 1, jl - 1]*zlevap
        
        zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] + zicefrac[jk - 1, jl - 1]*zlevap
        zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] - zicefrac[jk - 1, jl - 1]*zlevap
        
      
    
    #----------------------------------------------------------------------
    # 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
    #----------------------------------------------------------------------
    # (1) Increase of cloud water in existing clouds
    for jl in range(kidia, kfdia + 1):
      if za[jk - 1, jl - 1] > zepsec and zdqs[jl - 1] <= -yrecldp.rlmin:
        
        zlcond1[jl - 1] = max(-zdqs[jl - 1], 0.0)          #new limiter
        
        #old limiter (significantly improves upper tropospheric humidity rms)
        if za[jk - 1, jl - 1] > 0.99:
          zcor = 1.0 / (1.0 - ydcst.retv*zqsmix[jk - 1, jl - 1])
          zcdmax = (zqx[ncldqv - 1, jk - 1, jl - 1] - zqsmix[jk - 1, jl - 1]) / (1.0 + zcor*zqsmix[jk - 1, jl - 1]*foedem(ztp1[jk - 1, jl - 1]))
        else:
          zcdmax = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsmix[jk - 1, jl - 1]) / za[jk - 1, jl - 1]
        zlcond1[jl - 1] = max(min(zlcond1[jl - 1], zcdmax), 0.0)
        # end old limiter
        
        zlcond1[jl - 1] = za[jk - 1, jl - 1]*zlcond1[jl - 1]
        if zlcond1[jl - 1] < yrecldp.rlmin:
          zlcond1[jl - 1] = 0.0
        
        #-------------------------------------------------------------------------
        # All increase goes into liquid unless so cold cloud homogeneously freezes
        # Include new liquid formation in first guess value, otherwise liquid
        # remains at cold temperatures until next timestep.
        #-------------------------------------------------------------------------
        if ztp1[jk - 1, jl - 1] > yrecldp.rthomo:
          zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = zsolqa[ncldqv - 1, ncldql - 1, jl - 1] + zlcond1[jl - 1]
          zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zsolqa[ncldql - 1, ncldqv - 1, jl - 1] - zlcond1[jl - 1]
          zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] + zlcond1[jl - 1]
        else:
          zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] + zlcond1[jl - 1]
          zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] - zlcond1[jl - 1]
          zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zlcond1[jl - 1]
    
    # (2) Generation of new clouds (da/dt>0)
    
    for jl in range(kidia, kfdia + 1):
      
      if zdqs[jl - 1] <= -yrecldp.rlmin and za[jk - 1, jl - 1] < 1.0 - zepsec:
        
        #---------------------------
        # Critical relative humidity
        #---------------------------
        zrhc = yrecldp.ramid
        zsigk = pap[jk - 1, jl - 1] / paph[klev + 1 - 1, jl - 1]
        # Increase RHcrit to 1.0 towards the surface (eta>0.8)
        if zsigk > 0.8:
          zrhc = yrecldp.ramid + (1.0 - yrecldp.ramid)*((zsigk - 0.8) / 0.2)**2
        
        # Commented out for CY37R1 to reduce humidity in high trop and strat
        #      ! Increase RHcrit to 1.0 towards the tropopause (trop-0.2) and above
        #      ZBOTT=ZTRPAUS(JL)+0.2_JPRB
        #      IF(ZSIGK < ZBOTT) THEN
        #        ZRHC=RAMID+(1.0_JPRB-RAMID)*MIN(((ZBOTT-ZSIGK)/0.2_JPRB)**2,1.0_JPRB)
        #      ENDIF
        
        #---------------------------
        # Supersaturation options
        #---------------------------
        if yrecldp.nssopt == 0:
          # No scheme
          zqe = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsice[jk - 1, jl - 1]) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
          zqe = max(0.0, zqe)
        elif yrecldp.nssopt == 1:
          # Tompkins
          zqe = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsice[jk - 1, jl - 1]) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
          zqe = max(0.0, zqe)
        elif yrecldp.nssopt == 2:
          # Lohmann and Karcher
          zqe = zqx[ncldqv - 1, jk - 1, jl - 1]
        elif yrecldp.nssopt == 3:
          # Gierens
          zqe = zqx[ncldqv - 1, jk - 1, jl - 1] + zli[jk - 1, jl - 1]
        
        if ztp1[jk - 1, jl - 1] >= ydcst.rtt or yrecldp.nssopt == 0:
          # No ice supersaturation allowed
          zfac = 1.0
        else:
          # Ice supersaturation
          zfac = zfokoop[jl - 1]
        
        if zqe >= zrhc*zqsice[jk - 1, jl - 1]*zfac and zqe < zqsice[jk - 1, jl - 1]*zfac:
          # note: not **2 on 1-a term if ZQE is used.
          # Added correction term ZFAC to numerator 15/03/2010
          zacond = -(1.0 - za[jk - 1, jl - 1])*zfac*zdqs[jl - 1] / max(2.0*(zfac*zqsice[jk - 1, jl - 1] - zqe), zepsec)
          
          zacond = min(zacond, 1.0 - za[jk - 1, jl - 1])            #PUT THE LIMITER BACK
          
          # Linear term:
          # Added correction term ZFAC 15/03/2010
          zlcond2[jl - 1] = -zfac*zdqs[jl - 1]*0.5*zacond            #mine linear
          
          # new limiter formulation
          zzdl = 2.0*(zfac*zqsice[jk - 1, jl - 1] - zqe) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
          # Added correction term ZFAC 15/03/2010
          if zfac*zdqs[jl - 1] < -zzdl:
            # ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
            zlcondlim = (za[jk - 1, jl - 1] - 1.0)*zfac*zdqs[jl - 1] - zfac*zqsice[jk - 1, jl - 1] + zqx[ncldqv - 1, jk - 1, jl - 1]
            zlcond2[jl - 1] = min(zlcond2[jl - 1], zlcondlim)
          zlcond2[jl - 1] = max(zlcond2[jl - 1], 0.0)
          
          if zlcond2[jl - 1] < yrecldp.rlmin or (1.0 - za[jk - 1, jl - 1]) < zepsec:
            zlcond2[jl - 1] = 0.0
            zacond = 0.0
          if zlcond2[jl - 1] == 0.0:
            zacond = 0.0
          
          # Large-scale generation is LINEAR in A and LINEAR in L
          zsolac[jl - 1] = zsolac[jl - 1] + zacond            #linear
          
          #------------------------------------------------------------------------
          # All increase goes into liquid unless so cold cloud homogeneously freezes
          # Include new liquid formation in first guess value, otherwise liquid
          # remains at cold temperatures until next timestep.
          #------------------------------------------------------------------------
          if ztp1[jk - 1, jl - 1] > yrecldp.rthomo:
            zsolqa[ncldqv - 1, ncldql - 1, jl - 1] = zsolqa[ncldqv - 1, ncldql - 1, jl - 1] + zlcond2[jl - 1]
            zsolqa[ncldql - 1, ncldqv - 1, jl - 1] = zsolqa[ncldql - 1, ncldqv - 1, jl - 1] - zlcond2[jl - 1]
            zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] + zlcond2[jl - 1]
          else:
            # homogeneous freezing
            zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqi - 1, jl - 1] + zlcond2[jl - 1]
            zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqi - 1, ncldqv - 1, jl - 1] - zlcond2[jl - 1]
            zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zlcond2[jl - 1]
          
    
    #----------------------------------------------------------------------
    # 3.7 Growth of ice by vapour deposition
    #----------------------------------------------------------------------
    # Following Rotstayn et al. 2001:
    # does not use the ice nuclei number from cloudaer.F90
    # but rather a simple Meyers et al. 1992 form based on the
    # supersaturation and assuming clouds are saturated with
    # respect to liquid water (well mixed), (or Koop adjustment)
    # Growth considered as sink of liquid water if present so
    # Bergeron-Findeisen adjustment in autoconversion term no longer needed
    #----------------------------------------------------------------------
    
    #--------------------------------------------------------
    #-
    #- Ice deposition following Rotstayn et al. (2001)
    #-  (monodisperse ice particle size distribution)
    #-
    #--------------------------------------------------------
    if idepice == 1:
      
      for jl in range(kidia, kfdia + 1):
        
        #--------------------------------------------------------------
        # Calculate distance from cloud top
        # defined by cloudy layer below a layer with cloud frac <0.01
        # ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
        #--------------------------------------------------------------
        
        if za[jk - 1 - 1, jl - 1] < yrecldp.rcldtopcf and za[jk - 1, jl - 1] >= yrecldp.rcldtopcf:
          zcldtopdist[jl - 1] = 0.0
        else:
          zcldtopdist[jl - 1] = zcldtopdist[jl - 1] + zdp[jl - 1] / (zrho[jl - 1]*ydcst.rg)
        
        #--------------------------------------------------------------
        # only treat depositional growth if liquid present. due to fact
        # that can not model ice growth from vapour without additional
        # in-cloud water vapour variable
        #--------------------------------------------------------------
        if ztp1[jk - 1, jl - 1] < ydcst.rtt and zqxfg[ncldql - 1, jl - 1] > yrecldp.rlmin:
          # T<273K
          
          zvpice = foeeice(ztp1[jk - 1, jl - 1])*ydcst.rv / ydcst.rd
          zvpliq = zvpice*zfokoop[jl - 1]
          zicenuclei[jl - 1] = 1000.0*np.exp(12.96*(zvpliq - zvpice) / zvpliq - 0.639)
          
          #------------------------------------------------
          #   2.4e-2 is conductivity of air
          #   8.8 = 700**1/3 = density of ice to the third
          #------------------------------------------------
          zadd = ydcst.rlstt*(ydcst.rlstt / (ydcst.rv*ztp1[jk - 1, jl - 1]) - 1.0) / (2.4E-2*ztp1[jk - 1, jl - 1])
          zbdd = ydcst.rv*ztp1[jk - 1, jl - 1]*pap[jk - 1, jl - 1] / (2.21*zvpice)
          zcvds = 7.8*(zicenuclei[jl - 1] / zrho[jl - 1])**0.666*(zvpliq - zvpice) / (8.87*(zadd + zbdd)*zvpice)
          
          #-----------------------------------------------------
          # RICEINIT=1.E-12_JPRB is initial mass of ice particle
          #-----------------------------------------------------
          zice0 = max(zicecld[jl - 1], zicenuclei[jl - 1]*yrecldp.riceinit / zrho[jl - 1])
          
          #------------------
          # new value of ice:
          #------------------
          zinew = (0.666*zcvds*ptsphy + zice0**0.666)**1.5
          
          #---------------------------
          # grid-mean deposition rate:
          #---------------------------
          zdepos = max(za[jk - 1, jl - 1]*(zinew - zice0), 0.0)
          
          #--------------------------------------------------------------------
          # Limit deposition to liquid water amount
          # If liquid is all frozen, ice would use up reservoir of water
          # vapour in excess of ice saturation mixing ratio - However this
          # can not be represented without a in-cloud humidity variable. Using
          # the grid-mean humidity would imply a large artificial horizontal
          # flux from the clear sky to the cloudy area. We thus rely on the
          # supersaturation check to clean up any remaining supersaturation
          #--------------------------------------------------------------------
          zdepos = min(zdepos, zqxfg[ncldql - 1, jl - 1])            # limit to liquid water amount
          
          #--------------------------------------------------------------------
          # At top of cloud, reduce deposition rate near cloud top to account for
          # small scale turbulent processes, limited ice nucleation and ice fallout
          #--------------------------------------------------------------------
          #      ZDEPOS = ZDEPOS*MIN(RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH,1.0_JPRB)
          # Change to include dependence on ice nuclei concentration
          # to increase deposition rate with decreasing temperatures
          zinfactor = min(zicenuclei[jl - 1] / 15000., 1.0)
          zdepos = zdepos*min(zinfactor + (1.0 - zinfactor)*(yrecldp.rdepliqrefrate + zcldtopdist[jl - 1] / yrecldp.rdepliqrefdepth), 1.0)
          
          #--------------
          # add to matrix
          #--------------
          zsolqa[ncldql - 1, ncldqi - 1, jl - 1] = zsolqa[ncldql - 1, ncldqi - 1, jl - 1] + zdepos
          zsolqa[ncldqi - 1, ncldql - 1, jl - 1] = zsolqa[ncldqi - 1, ncldql - 1, jl - 1] - zdepos
          zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zdepos
          zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] - zdepos
          
      
      #--------------------------------------------------------
      #-
      #- Ice deposition assuming ice PSD
      #-
      #--------------------------------------------------------
    elif idepice == 2:
      
      for jl in range(kidia, kfdia + 1):
        
        #--------------------------------------------------------------
        # Calculate distance from cloud top
        # defined by cloudy layer below a layer with cloud frac <0.01
        # ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
        #--------------------------------------------------------------
        
        if za[jk - 1 - 1, jl - 1] < yrecldp.rcldtopcf and za[jk - 1, jl - 1] >= yrecldp.rcldtopcf:
          zcldtopdist[jl - 1] = 0.0
        else:
          zcldtopdist[jl - 1] = zcldtopdist[jl - 1] + zdp[jl - 1] / (zrho[jl - 1]*ydcst.rg)
        
        #--------------------------------------------------------------
        # only treat depositional growth if liquid present. due to fact
        # that can not model ice growth from vapour without additional
        # in-cloud water vapour variable
        #--------------------------------------------------------------
        if ztp1[jk - 1, jl - 1] < ydcst.rtt and zqxfg[ncldql - 1, jl - 1] > yrecldp.rlmin:
          # T<273K
          
          zvpice = foeeice(ztp1[jk - 1, jl - 1])*ydcst.rv / ydcst.rd
          zvpliq = zvpice*zfokoop[jl - 1]
          zicenuclei[jl - 1] = 1000.0*np.exp(12.96*(zvpliq - zvpice) / zvpliq - 0.639)
          
          #-----------------------------------------------------
          # RICEINIT=1.E-12_JPRB is initial mass of ice particle
          #-----------------------------------------------------
          zice0 = max(zicecld[jl - 1], zicenuclei[jl - 1]*yrecldp.riceinit / zrho[jl - 1])
          
          # Particle size distribution
          ztcg = 1.0
          zfacx1i = 1.0
          
          zaplusb = yrecldp.rcl_apb1*zvpice - yrecldp.rcl_apb2*zvpice*ztp1[jk - 1, jl - 1] + pap[jk - 1, jl - 1]*yrecldp.rcl_apb3*ztp1[jk - 1, jl - 1]**3.
          zcorrfac = (1.0 / zrho[jl - 1])**0.5
          zcorrfac2 = ((ztp1[jk - 1, jl - 1] / 273.0)**1.5)*(393.0 / (ztp1[jk - 1, jl - 1] + 120.0))
          
          zpr02 = zrho[jl - 1]*zice0*yrecldp.rcl_const1i / (ztcg*zfacx1i)
          
          zterm1 = (zvpliq - zvpice)*ztp1[jk - 1, jl - 1]**2.0*zvpice*zcorrfac2*ztcg*yrecldp.rcl_const2i*zfacx1i / (zrho[jl - 1]*zaplusb*zvpice)
          zterm2 = 0.65*yrecldp.rcl_const6i*zpr02**yrecldp.rcl_const4i + yrecldp.rcl_const3i*zcorrfac**0.5*zrho[jl - 1]**0.5*zpr02**yrecldp.rcl_const5i / zcorrfac2**0.5
          
          zdepos = max(za[jk - 1, jl - 1]*zterm1*zterm2*ptsphy, 0.0)
          
          #--------------------------------------------------------------------
          # Limit deposition to liquid water amount
          # If liquid is all frozen, ice would use up reservoir of water
          # vapour in excess of ice saturation mixing ratio - However this
          # can not be represented without a in-cloud humidity variable. Using
          # the grid-mean humidity would imply a large artificial horizontal
          # flux from the clear sky to the cloudy area. We thus rely on the
          # supersaturation check to clean up any remaining supersaturation
          #--------------------------------------------------------------------
          zdepos = min(zdepos, zqxfg[ncldql - 1, jl - 1])            # limit to liquid water amount
          
          #--------------------------------------------------------------------
          # At top of cloud, reduce deposition rate near cloud top to account for
          # small scale turbulent processes, limited ice nucleation and ice fallout
          #--------------------------------------------------------------------
          # Change to include dependence on ice nuclei concentration
          # to increase deposition rate with decreasing temperatures
          zinfactor = min(zicenuclei[jl - 1] / 15000., 1.0)
          zdepos = zdepos*min(zinfactor + (1.0 - zinfactor)*(yrecldp.rdepliqrefrate + zcldtopdist[jl - 1] / yrecldp.rdepliqrefdepth), 1.0)
          
          #--------------
          # add to matrix
          #--------------
          zsolqa[ncldql - 1, ncldqi - 1, jl - 1] = zsolqa[ncldql - 1, ncldqi - 1, jl - 1] + zdepos
          zsolqa[ncldqi - 1, ncldql - 1, jl - 1] = zsolqa[ncldqi - 1, ncldql - 1, jl - 1] - zdepos
          zqxfg[ncldqi - 1, jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zdepos
          zqxfg[ncldql - 1, jl - 1] = zqxfg[ncldql - 1, jl - 1] - zdepos
      
    # on IDEPICE
    
    #######################################################################
    #              4  *** PRECIPITATION PROCESSES ***
    #######################################################################
    
    #----------------------------------
    # revise in-cloud consensate amount
    #----------------------------------
    for jl in range(kidia, kfdia + 1):
      ztmpa = 1.0 / max(za[jk - 1, jl - 1], zepsec)
      zliqcld[jl - 1] = zqxfg[ncldql - 1, jl - 1]*ztmpa
      zicecld[jl - 1] = zqxfg[ncldqi - 1, jl - 1]*ztmpa
      zlicld[jl - 1] = zliqcld[jl - 1] + zicecld[jl - 1]
    
    #----------------------------------------------------------------------
    # 4.2 SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
    #     now that rain, snow, graupel species are prognostic
    #     the precipitation flux can be defined directly level by level
    #     There is no vertical memory required from the flux variable
    #----------------------------------------------------------------------
    
    for jm in range(1, nclv + 1):
      if llfall[jm - 1] or jm == ncldqi:
        for jl in range(kidia, kfdia + 1):
          #------------------------
          # source from layer above
          #------------------------
          if jk > yrecldp.ncldtop:
            zfallsrce[jm - 1, jl - 1] = zpfplsx[jm - 1, jk - 1, jl - 1]*zdtgdp[jl - 1]
            zsolqa[jm - 1, jm - 1, jl - 1] = zsolqa[jm - 1, jm - 1, jl - 1] + zfallsrce[jm - 1, jl - 1]
            zqxfg[jm - 1, jl - 1] = zqxfg[jm - 1, jl - 1] + zfallsrce[jm - 1, jl - 1]
            # use first guess precip----------V
            zqpretot[jl - 1] = zqpretot[jl - 1] + zqxfg[jm - 1, jl - 1]
          #-------------------------------------------------
          # sink to next layer, constant fall speed
          #-------------------------------------------------
          # if aerosol effect then override
          #  note that for T>233K this is the same as above.
          if yrecldp.laericesed and jm == ncldqi:
            zre_ice = pre_ice[jk - 1, jl - 1]
            # The exponent value is from
            # Morrison et al. JAS 2005 Appendix
            zvqx[ncldqi - 1] = 0.002*zre_ice**1.0
          zfall = zvqx[jm - 1]*zrho[jl - 1]
          #-------------------------------------------------
          # modified by Heymsfield and Iaquinta JAS 2000
          #-------------------------------------------------
          # ZFALL = ZFALL*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
          #            &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))
          
          zfallsink[jm - 1, jl - 1] = zdtgdp[jl - 1]*zfall
          # Cloud budget diagnostic stored at end as implicit
        # jl
      # LLFALL
    # jm
    
    #---------------------------------------------------------------
    # Precip cover overlap using MAX-RAN Overlap
    # Since precipitation is now prognostic we must
    #   1) apply an arbitrary minimum coverage (0.3) if precip>0
    #   2) abandon the 2-flux clr/cld treatment
    #   3) Thus, since we have no memory of the clear sky precip
    #      fraction, we mimic the previous method by reducing
    #      ZCOVPTOT(JL), which has the memory, proportionally with
    #      the precip evaporation rate, taking cloud fraction
    #      into account
    #   #3 above leads to much smoother vertical profiles of
    #   precipitation fraction than the Klein-Jakob scheme which
    #   monotonically increases precip fraction and then resets
    #   it to zero in a step function once clear-sky precip reaches
    #   zero.
    #---------------------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      if zqpretot[jl - 1] > zepsec:
        zcovptot[jl - 1] = 1.0 - ((1.0 - zcovptot[jl - 1])*(1.0 - max(za[jk - 1, jl - 1], za[jk - 1 - 1, jl - 1])) / (1.0 - min(za[jk - 1 - 1, jl - 1], 1.0 - 1.E-06)))
        zcovptot[jl - 1] = max(zcovptot[jl - 1], yrecldp.rcovpmin)
        zcovpclr[jl - 1] = max(0.0, zcovptot[jl - 1] - za[jk - 1, jl - 1])          # clear sky proportion
        zraincld[jl - 1] = zqxfg[ncldqr - 1, jl - 1] / zcovptot[jl - 1]
        zsnowcld[jl - 1] = zqxfg[ncldqs - 1, jl - 1] / zcovptot[jl - 1]
        zcovpmax[jl - 1] = max(zcovptot[jl - 1], zcovpmax[jl - 1])
      else:
        zraincld[jl - 1] = 0.0
        zsnowcld[jl - 1] = 0.0
        zcovptot[jl - 1] = 0.0          # no flux - reset cover
        zcovpclr[jl - 1] = 0.0          # reset clear sky proportion
        zcovpmax[jl - 1] = 0.0          # reset max cover for ZZRH calc
    
    #----------------------------------------------------------------------
    # 4.3a AUTOCONVERSION TO SNOW
    #----------------------------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      
      if ztp1[jk - 1, jl - 1] <= ydcst.rtt:
        #-----------------------------------------------------
        #     Snow Autoconversion rate follow Lin et al. 1983
        #-----------------------------------------------------
        if zicecld[jl - 1] > zepsec:
          
          zzco = ptsphy*yrecldp.rsnowlin1*np.exp(yrecldp.rsnowlin2*(ztp1[jk - 1, jl - 1] - ydcst.rtt))
          
          if yrecldp.laericeauto:
            zlcrit = picrit_aer[jk - 1, jl - 1]
            # 0.3 = N**0.333 with N=0.027
            zzco = zzco*(yrecldp.rnice / pnice[jk - 1, jl - 1])**0.333
          else:
            zlcrit = yrecldp.rlcritsnow
          
          zsnowaut[jl - 1] = zzco*(1.0 - np.exp(-(zicecld[jl - 1] / zlcrit)**2))
          zsolqb[ncldqi - 1, ncldqs - 1, jl - 1] = zsolqb[ncldqi - 1, ncldqs - 1, jl - 1] + zsnowaut[jl - 1]
          
      
      #----------------------------------------------------------------------
      # 4.3b AUTOCONVERSION WARM CLOUDS
      #   Collection and accretion will require separate treatment
      #   but for now we keep this simple treatment
      #----------------------------------------------------------------------
      
      if zliqcld[jl - 1] > zepsec:
        
        #--------------------------------------------------------
        #-
        #- Warm-rain process follow Sundqvist (1989)
        #-
        #--------------------------------------------------------
        if iwarmrain == 1:
          
          zzco = yrecldp.rkconv*ptsphy
          
          if yrecldp.laerliqautolsp:
            zlcrit = plcrit_aer[jk - 1, jl - 1]
            # 0.3 = N**0.333 with N=125 cm-3
            zzco = zzco*(yrecldp.rccn / pccn[jk - 1, jl - 1])**0.333
          else:
            # Modify autoconversion threshold dependent on:
            #  land (polluted, high CCN, smaller droplets, higher threshold)
            #  sea  (clean, low CCN, larger droplets, lower threshold)
            if plsm[jl - 1] > 0.5:
              zlcrit = yrecldp.rclcrit_land                # land
            else:
              zlcrit = yrecldp.rclcrit_sea                # ocean
          
          #------------------------------------------------------------------
          # Parameters for cloud collection by rain and snow.
          # Note that with new prognostic variable it is now possible
          # to REPLACE this with an explicit collection parametrization
          #------------------------------------------------------------------
          zprecip = (zpfplsx[ncldqs - 1, jk - 1, jl - 1] + zpfplsx[ncldqr - 1, jk - 1, jl - 1]) / max(zepsec, zcovptot[jl - 1])
          zcfpr = 1.0 + yrecldp.rprc1*np.sqrt(max(zprecip, 0.0))
          #      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
          #       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))
          
          if yrecldp.laerliqcoll:
            # 5.0 = N**0.333 with N=125 cm-3
            zcfpr = zcfpr*(yrecldp.rccn / pccn[jk - 1, jl - 1])**0.333
          
          zzco = zzco*zcfpr
          zlcrit = zlcrit / max(zcfpr, zepsec)
          
          if zliqcld[jl - 1] / zlcrit < 20.0:
            # Security for exp for some compilers
            zrainaut[jl - 1] = zzco*(1.0 - np.exp(-(zliqcld[jl - 1] / zlcrit)**2))
          else:
            zrainaut[jl - 1] = zzco
          
          # rain freezes instantly
          if ztp1[jk - 1, jl - 1] <= ydcst.rtt:
            zsolqb[ncldql - 1, ncldqs - 1, jl - 1] = zsolqb[ncldql - 1, ncldqs - 1, jl - 1] + zrainaut[jl - 1]
          else:
            zsolqb[ncldql - 1, ncldqr - 1, jl - 1] = zsolqb[ncldql - 1, ncldqr - 1, jl - 1] + zrainaut[jl - 1]
          
          #--------------------------------------------------------
          #-
          #- Warm-rain process follow Khairoutdinov and Kogan (2000)
          #-
          #--------------------------------------------------------
        elif iwarmrain == 2:
          
          if plsm[jl - 1] > 0.5:
            # land
            zconst = yrecldp.rcl_kk_cloud_num_land
            zlcrit = yrecldp.rclcrit_land
          else:
            # ocean
            zconst = yrecldp.rcl_kk_cloud_num_sea
            zlcrit = yrecldp.rclcrit_sea
          
          if zliqcld[jl - 1] > zlcrit:
            
            zrainaut[jl - 1] = 1.5*za[jk - 1, jl - 1]*ptsphy*yrecldp.rcl_kkaau*zliqcld[jl - 1]**yrecldp.rcl_kkbauq*zconst**yrecldp.rcl_kkbaun
            
            zrainaut[jl - 1] = min(zrainaut[jl - 1], zqxfg[ncldql - 1, jl - 1])
            if zrainaut[jl - 1] < zepsec:
              zrainaut[jl - 1] = 0.0
            
            zrainacc[jl - 1] = 2.0*za[jk - 1, jl - 1]*ptsphy*yrecldp.rcl_kkaac*(zliqcld[jl - 1]*zraincld[jl - 1])**yrecldp.rcl_kkbac
            
            zrainacc[jl - 1] = min(zrainacc[jl - 1], zqxfg[ncldql - 1, jl - 1])
            if zrainacc[jl - 1] < zepsec:
              zrainacc[jl - 1] = 0.0
            
          else:
            zrainaut[jl - 1] = 0.0
            zrainacc[jl - 1] = 0.0
          
          # If temperature < 0, then autoconversion produces snow rather than rain
          # Explicit
          if ztp1[jk - 1, jl - 1] <= ydcst.rtt:
            zsolqa[ncldql - 1, ncldqs - 1, jl - 1] = zsolqa[ncldql - 1, ncldqs - 1, jl - 1] + zrainaut[jl - 1]
            zsolqa[ncldql - 1, ncldqs - 1, jl - 1] = zsolqa[ncldql - 1, ncldqs - 1, jl - 1] + zrainacc[jl - 1]
            zsolqa[ncldqs - 1, ncldql - 1, jl - 1] = zsolqa[ncldqs - 1, ncldql - 1, jl - 1] - zrainaut[jl - 1]
            zsolqa[ncldqs - 1, ncldql - 1, jl - 1] = zsolqa[ncldqs - 1, ncldql - 1, jl - 1] - zrainacc[jl - 1]
          else:
            zsolqa[ncldql - 1, ncldqr - 1, jl - 1] = zsolqa[ncldql - 1, ncldqr - 1, jl - 1] + zrainaut[jl - 1]
            zsolqa[ncldql - 1, ncldqr - 1, jl - 1] = zsolqa[ncldql - 1, ncldqr - 1, jl - 1] + zrainacc[jl - 1]
            zsolqa[ncldqr - 1, ncldql - 1, jl - 1] = zsolqa[ncldqr - 1, ncldql - 1, jl - 1] - zrainaut[jl - 1]
            zsolqa[ncldqr - 1, ncldql - 1, jl - 1] = zsolqa[ncldqr - 1, ncldql - 1, jl - 1] - zrainacc[jl - 1]
          
        # on IWARMRAIN
        
      # on ZLIQCLD > ZEPSEC
    
    
    #----------------------------------------------------------------------
    # RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
    #      only active if T<0degC and supercooled liquid water is present
    #      AND if not Sundquist autoconversion (as this includes riming)
    #----------------------------------------------------------------------
    if iwarmrain > 1:
      
      for jl in range(kidia, kfdia + 1):
        if ztp1[jk - 1, jl - 1] <= ydcst.rtt and zliqcld[jl - 1] > zepsec:
          
          # Fallspeed air density correction
          zfallcorr = (yrecldp.rdensref / zrho[jl - 1])**0.4
          
          #------------------------------------------------------------------
          # Riming of snow by cloud water - implicit in lwc
          #------------------------------------------------------------------
          if zsnowcld[jl - 1] > zepsec and zcovptot[jl - 1] > 0.01:
            
            # Calculate riming term
            # Factor of liq water taken out because implicit
            zsnowrime[jl - 1] = 0.3*zcovptot[jl - 1]*ptsphy*yrecldp.rcl_const7s*zfallcorr*(zrho[jl - 1]*zsnowcld[jl - 1]*yrecldp.rcl_const1s)**yrecldp.rcl_const8s
            
            # Limit snow riming term
            zsnowrime[jl - 1] = min(zsnowrime[jl - 1], 1.0)
            
            zsolqb[ncldql - 1, ncldqs - 1, jl - 1] = zsolqb[ncldql - 1, ncldqs - 1, jl - 1] + zsnowrime[jl - 1]
            
          
          #------------------------------------------------------------------
          # Riming of ice by cloud water - implicit in lwc
          # NOT YET ACTIVE
          #------------------------------------------------------------------
          #      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
          #
          #        ! Calculate riming term
          #        ! Factor of liq water taken out because implicit
          #        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
          #     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
          #
          #        ! Limit ice riming term
          #        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
          #
          #        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
          #
          #      ENDIF
      
    # on IWARMRAIN > 1
    
    
    #----------------------------------------------------------------------
    # 4.4a  MELTING OF SNOW and ICE
    #       with new implicit solver this also has to treat snow or ice
    #       precipitating from the level above... i.e. local ice AND flux.
    #       in situ ice and snow: could arise from LS advection or warming
    #       falling ice and snow: arrives by precipitation process
    #----------------------------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      
      zicetot[jl - 1] = zqxfg[ncldqi - 1, jl - 1] + zqxfg[ncldqs - 1, jl - 1]
      zmeltmax[jl - 1] = 0.0
      
      # If there are frozen hydrometeors present and dry-bulb temperature > 0degC
      if zicetot[jl - 1] > zepsec and ztp1[jk - 1, jl - 1] > ydcst.rtt:
        
        # Calculate subsaturation
        zsubsat = max(zqsice[jk - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1], 0.0)
        
        # Calculate difference between dry-bulb (ZTP1) and the temperature
        # at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
        # Melting only occurs if the wet-bulb temperature >0
        # i.e. warming of ice particle due to melting > cooling
        # due to evaporation.
        ztdmtw0 = ztp1[jk - 1, jl - 1] - ydcst.rtt - zsubsat*(ztw1 + ztw2*(pap[jk - 1, jl - 1] - ztw3) - ztw4*(ztp1[jk - 1, jl - 1] - ztw5))
        # Not implicit yet...
        # Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
        zcons1 = abs(ptsphy*(1.0 + 0.5*ztdmtw0) / yrecldp.rtaumel)
        zmeltmax[jl - 1] = max(ztdmtw0*zcons1*zrldcp, 0.0)
    
    # Loop over frozen hydrometeors (ice, snow)
    for jm in range(1, nclv + 1):
      if iphase[jm - 1] == 2:
        jn = imelt[jm - 1]
        for jl in range(kidia, kfdia + 1):
          if zmeltmax[jl - 1] > zepsec and zicetot[jl - 1] > zepsec:
            # Apply melting in same proportion as frozen hydrometeor fractions
            zalfa = zqxfg[jm - 1, jl - 1] / zicetot[jl - 1]
            zmelt = min(zqxfg[jm - 1, jl - 1], zalfa*zmeltmax[jl - 1])
            # needed in first guess
            # This implies that zqpretot has to be recalculated below
            # since is not conserved here if ice falls and liquid doesn't
            zqxfg[jm - 1, jl - 1] = zqxfg[jm - 1, jl - 1] - zmelt
            zqxfg[jn - 1, jl - 1] = zqxfg[jn - 1, jl - 1] + zmelt
            zsolqa[jm - 1, jn - 1, jl - 1] = zsolqa[jm - 1, jn - 1, jl - 1] + zmelt
            zsolqa[jn - 1, jm - 1, jl - 1] = zsolqa[jn - 1, jm - 1, jl - 1] - zmelt
    
    #----------------------------------------------------------------------
    # 4.4b  FREEZING of RAIN
    #----------------------------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      
      # If rain present
      if zqx[ncldqr - 1, jk - 1, jl - 1] > zepsec:
        
        if ztp1[jk - 1, jl - 1] <= ydcst.rtt and ztp1[jk - 1 - 1, jl - 1] > ydcst.rtt:
          # Base of melting layer/top of refreezing layer so
          # store rain/snow fraction for precip type diagnosis
          # If mostly rain, then supercooled rain slow to freeze
          # otherwise faster to freeze (snow or ice pellets)
          zqpretot[jl - 1] = max(zqx[ncldqs - 1, jk - 1, jl - 1] + zqx[ncldqr - 1, jk - 1, jl - 1], zepsec)
          prainfrac_toprfz[jl - 1] = zqx[ncldqr - 1, jk - 1, jl - 1] / zqpretot[jl - 1]
          if prainfrac_toprfz[jl - 1] > 0.8:
            llrainliq[jl - 1] = True
          else:
            llrainliq[jl - 1] = False
        
        # If temperature less than zero
        if ztp1[jk - 1, jl - 1] < ydcst.rtt:
          
          if prainfrac_toprfz[jl - 1] > 0.8:
            
            # Majority of raindrops completely melted
            # Refreezing is by slow heterogeneous freezing
            
            # Slope of rain particle size distribution
            zlambda = (yrecldp.rcl_fac1 / (zrho[jl - 1]*zqx[ncldqr - 1, jk - 1, jl - 1]))**yrecldp.rcl_fac2
            
            # Calculate freezing rate based on Bigg(1953) and Wisner(1972)
            ztemp = yrecldp.rcl_fzrab*(ztp1[jk - 1, jl - 1] - ydcst.rtt)
            zfrz = ptsphy*(yrecldp.rcl_const5r / zrho[jl - 1])*(np.exp(ztemp) - 1.)*zlambda**yrecldp.rcl_const6r
            zfrzmax[jl - 1] = max(zfrz, 0.0)
            
          else:
            
            # Majority of raindrops only partially melted
            # Refreeze with a shorter timescale (reverse of melting...for now)
            
            zcons1 = abs(ptsphy*(1.0 + 0.5*(ydcst.rtt - ztp1[jk - 1, jl - 1])) / yrecldp.rtaumel)
            zfrzmax[jl - 1] = max((ydcst.rtt - ztp1[jk - 1, jl - 1])*zcons1*zrldcp, 0.0)
            
          
          if zfrzmax[jl - 1] > zepsec:
            zfrz = min(zqx[ncldqr - 1, jk - 1, jl - 1], zfrzmax[jl - 1])
            zsolqa[ncldqr - 1, ncldqs - 1, jl - 1] = zsolqa[ncldqr - 1, ncldqs - 1, jl - 1] + zfrz
            zsolqa[ncldqs - 1, ncldqr - 1, jl - 1] = zsolqa[ncldqs - 1, ncldqr - 1, jl - 1] - zfrz
        
      
    
    #----------------------------------------------------------------------
    # 4.4c  FREEZING of LIQUID
    #----------------------------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      # not implicit yet...
      zfrzmax[jl - 1] = max((yrecldp.rthomo - ztp1[jk - 1, jl - 1])*zrldcp, 0.0)
    
    jm = ncldql
    jn = imelt[jm - 1]
    for jl in range(kidia, kfdia + 1):
      if zfrzmax[jl - 1] > zepsec and zqxfg[jm - 1, jl - 1] > zepsec:
        zfrz = min(zqxfg[jm - 1, jl - 1], zfrzmax[jl - 1])
        zsolqa[jm - 1, jn - 1, jl - 1] = zsolqa[jm - 1, jn - 1, jl - 1] + zfrz
        zsolqa[jn - 1, jm - 1, jl - 1] = zsolqa[jn - 1, jm - 1, jl - 1] - zfrz
    
    #----------------------------------------------------------------------
    # 4.5   EVAPORATION OF RAIN/SNOW
    #----------------------------------------------------------------------
    
    #----------------------------------------
    # Rain evaporation scheme from Sundquist
    #----------------------------------------
    if ievaprain == 1:
      
      # Rain
      
      for jl in range(kidia, kfdia + 1):
        
        zzrh = yrecldp.rprecrhmax + (1.0 - yrecldp.rprecrhmax)*zcovpmax[jl - 1] / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        zzrh = min(max(zzrh, yrecldp.rprecrhmax), 1.0)
        
        zqe = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsliq[jk - 1, jl - 1]) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        #---------------------------------------------
        # humidity in moistest ZCOVPCLR part of domain
        #---------------------------------------------
        zqe = max(0.0, min(zqe, zqsliq[jk - 1, jl - 1]))
        llo1 = zcovpclr[jl - 1] > zepsec and zqxfg[ncldqr - 1, jl - 1] > zepsec and zqe < zzrh*zqsliq[jk - 1, jl - 1]
        
        if llo1:
          # note: zpreclr is a rain flux
          zpreclr = zqxfg[ncldqr - 1, jl - 1]*zcovpclr[jl - 1] / (max(abs(zcovptot[jl - 1]*zdtgdp[jl - 1]), zepsilon)*np.sign(zcovptot[jl - 1]*zdtgdp[jl - 1]))
          
          #--------------------------------------
          # actual microphysics formula in zbeta
          #--------------------------------------
          
          zbeta1 = np.sqrt(pap[jk - 1, jl - 1] / paph[klev + 1 - 1, jl - 1]) / yrecldp.rvrfactor*zpreclr / max(zcovpclr[jl - 1], zepsec)
          
          zbeta = ydcst.rg*yrecldp.rpecons*0.5*zbeta1**0.5777
          
          zdenom = 1.0 + zbeta*ptsphy*zcorqsliq[jl - 1]
          zdpr = zcovpclr[jl - 1]*zbeta*(zqsliq[jk - 1, jl - 1] - zqe) / zdenom*zdp[jl - 1]*zrg_r
          zdpevap = zdpr*zdtgdp[jl - 1]
          
          #---------------------------------------------------------
          # add evaporation term to explicit sink.
          # this has to be explicit since if treated in the implicit
          # term evaporation can not reduce rain to zero and model
          # produces small amounts of rainfall everywhere.
          #---------------------------------------------------------
          
          # Evaporate rain
          zevap = min(zdpevap, zqxfg[ncldqr - 1, jl - 1])
          
          zsolqa[ncldqr - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqr - 1, ncldqv - 1, jl - 1] + zevap
          zsolqa[ncldqv - 1, ncldqr - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqr - 1, jl - 1] - zevap
          
          #-------------------------------------------------------------
          # Reduce the total precip coverage proportional to evaporation
          # to mimic the previous scheme which had a diagnostic
          # 2-flux treatment, abandoned due to the new prognostic precip
          #-------------------------------------------------------------
          zcovptot[jl - 1] = max(yrecldp.rcovpmin, zcovptot[jl - 1] - max(0.0, (zcovptot[jl - 1] - za[jk - 1, jl - 1])*zevap / zqxfg[ncldqr - 1, jl - 1]))
          
          # Update fg field
          zqxfg[ncldqr - 1, jl - 1] = zqxfg[ncldqr - 1, jl - 1] - zevap
          
      
      
      #---------------------------------------------------------
      # Rain evaporation scheme based on Abel and Boutle (2013)
      #---------------------------------------------------------
    elif ievaprain == 2:
      
      for jl in range(kidia, kfdia + 1):
        
        #-----------------------------------------------------------------------
        # Calculate relative humidity limit for rain evaporation
        # to avoid cloud formation and saturation of the grid box
        #-----------------------------------------------------------------------
        # Limit RH for rain evaporation dependent on precipitation fraction
        zzrh = yrecldp.rprecrhmax + (1.0 - yrecldp.rprecrhmax)*zcovpmax[jl - 1] / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        zzrh = min(max(zzrh, yrecldp.rprecrhmax), 1.0)
        
        # Critical relative humidity
        #ZRHC=RAMID
        #ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
        # Increase RHcrit to 1.0 towards the surface (eta>0.8)
        #IF(ZSIGK > 0.8_JPRB) THEN
        #  ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
        #ENDIF
        #ZZRH = MIN(ZRHC,ZZRH)
        
        # Further limit RH for rain evaporation to 80% (RHcrit in free troposphere)
        zzrh = min(0.8, zzrh)
        
        zqe = max(0.0, min(zqx[ncldqv - 1, jk - 1, jl - 1], zqsliq[jk - 1, jl - 1]))
        
        llo1 = zcovpclr[jl - 1] > zepsec and zqxfg[ncldqr - 1, jl - 1] > zepsec and zqe < zzrh*zqsliq[jk - 1, jl - 1]
        
        if llo1:
          
          #-------------------------------------------
          # Abel and Boutle (2012) evaporation
          #-------------------------------------------
          # Calculate local precipitation (kg/kg)
          zpreclr = zqxfg[ncldqr - 1, jl - 1] / zcovptot[jl - 1]
          
          # Fallspeed air density correction
          zfallcorr = (yrecldp.rdensref / zrho[jl - 1])**0.4
          
          # Saturation vapour pressure with respect to liquid phase
          zesatliq = ydcst.rv / ydcst.rd*foeeliq(ztp1[jk - 1, jl - 1])
          
          # Slope of particle size distribution
          zlambda = (yrecldp.rcl_fac1 / (zrho[jl - 1]*zpreclr))**yrecldp.rcl_fac2            # ZPRECLR=kg/kg
          
          zevap_denom = yrecldp.rcl_cdenom1*zesatliq - yrecldp.rcl_cdenom2*ztp1[jk - 1, jl - 1]*zesatliq + yrecldp.rcl_cdenom3*ztp1[jk - 1, jl - 1]**3.*pap[jk - 1, jl - 1]
          
          # Temperature dependent conductivity
          zcorr2 = (ztp1[jk - 1, jl - 1] / 273.)**1.5*393. / (ztp1[jk - 1, jl - 1] + 120.)
          zka = yrecldp.rcl_ka273*zcorr2
          
          zsubsat = max(zzrh*zqsliq[jk - 1, jl - 1] - zqe, 0.0)
          
          zbeta = (0.5 / zqsliq[jk - 1, jl - 1])*ztp1[jk - 1, jl - 1]**2.*zesatliq*yrecldp.rcl_const1r*(zcorr2 / zevap_denom)*(0.78 / (zlambda**yrecldp.rcl_const4r) + yrecldp.rcl_const2r*(zrho[jl - 1]*zfallcorr)**0.5 / (zcorr2**0.5*zlambda**yrecldp.rcl_const3r))
          
          zdenom = 1.0 + zbeta*ptsphy            #*ZCORQSLIQ(JL)
          zdpevap = zcovpclr[jl - 1]*zbeta*ptsphy*zsubsat / zdenom
          
          #---------------------------------------------------------
          # Add evaporation term to explicit sink.
          # this has to be explicit since if treated in the implicit
          # term evaporation can not reduce rain to zero and model
          # produces small amounts of rainfall everywhere.
          #---------------------------------------------------------
          
          # Limit rain evaporation
          zevap = min(zdpevap, zqxfg[ncldqr - 1, jl - 1])
          
          zsolqa[ncldqr - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqr - 1, ncldqv - 1, jl - 1] + zevap
          zsolqa[ncldqv - 1, ncldqr - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqr - 1, jl - 1] - zevap
          
          #-------------------------------------------------------------
          # Reduce the total precip coverage proportional to evaporation
          # to mimic the previous scheme which had a diagnostic
          # 2-flux treatment, abandoned due to the new prognostic precip
          #-------------------------------------------------------------
          zcovptot[jl - 1] = max(yrecldp.rcovpmin, zcovptot[jl - 1] - max(0.0, (zcovptot[jl - 1] - za[jk - 1, jl - 1])*zevap / zqxfg[ncldqr - 1, jl - 1]))
          
          # Update fg field
          zqxfg[ncldqr - 1, jl - 1] = zqxfg[ncldqr - 1, jl - 1] - zevap
          
      
    # on IEVAPRAIN
    
    #----------------------------------------------------------------------
    # 4.5   EVAPORATION OF SNOW
    #----------------------------------------------------------------------
    # Snow
    if ievapsnow == 1:
      
      for jl in range(kidia, kfdia + 1):
        zzrh = yrecldp.rprecrhmax + (1.0 - yrecldp.rprecrhmax)*zcovpmax[jl - 1] / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        zzrh = min(max(zzrh, yrecldp.rprecrhmax), 1.0)
        zqe = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsice[jk - 1, jl - 1]) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        
        #---------------------------------------------
        # humidity in moistest ZCOVPCLR part of domain
        #---------------------------------------------
        zqe = max(0.0, min(zqe, zqsice[jk - 1, jl - 1]))
        llo1 = zcovpclr[jl - 1] > zepsec and zqxfg[ncldqs - 1, jl - 1] > zepsec and zqe < zzrh*zqsice[jk - 1, jl - 1]
        
        if llo1:
          # note: zpreclr is a rain flux a
          zpreclr = zqxfg[ncldqs - 1, jl - 1]*zcovpclr[jl - 1] / (max(abs(zcovptot[jl - 1]*zdtgdp[jl - 1]), zepsilon)*np.sign(zcovptot[jl - 1]*zdtgdp[jl - 1]))
          
          #--------------------------------------
          # actual microphysics formula in zbeta
          #--------------------------------------
          
          zbeta1 = np.sqrt(pap[jk - 1, jl - 1] / paph[klev + 1 - 1, jl - 1]) / yrecldp.rvrfactor*zpreclr / max(zcovpclr[jl - 1], zepsec)
          
          zbeta = ydcst.rg*yrecldp.rpecons*zbeta1**0.5777
          
          zdenom = 1.0 + zbeta*ptsphy*zcorqsice[jl - 1]
          zdpr = zcovpclr[jl - 1]*zbeta*(zqsice[jk - 1, jl - 1] - zqe) / zdenom*zdp[jl - 1]*zrg_r
          zdpevap = zdpr*zdtgdp[jl - 1]
          
          #---------------------------------------------------------
          # add evaporation term to explicit sink.
          # this has to be explicit since if treated in the implicit
          # term evaporation can not reduce snow to zero and model
          # produces small amounts of snowfall everywhere.
          #---------------------------------------------------------
          
          # Evaporate snow
          zevap = min(zdpevap, zqxfg[ncldqs - 1, jl - 1])
          
          zsolqa[ncldqs - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqs - 1, ncldqv - 1, jl - 1] + zevap
          zsolqa[ncldqv - 1, ncldqs - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqs - 1, jl - 1] - zevap
          
          #-------------------------------------------------------------
          # Reduce the total precip coverage proportional to evaporation
          # to mimic the previous scheme which had a diagnostic
          # 2-flux treatment, abandoned due to the new prognostic precip
          #-------------------------------------------------------------
          zcovptot[jl - 1] = max(yrecldp.rcovpmin, zcovptot[jl - 1] - max(0.0, (zcovptot[jl - 1] - za[jk - 1, jl - 1])*zevap / zqxfg[ncldqs - 1, jl - 1]))
          
          #Update first guess field
          zqxfg[ncldqs - 1, jl - 1] = zqxfg[ncldqs - 1, jl - 1] - zevap
          
      #---------------------------------------------------------
    elif ievapsnow == 2:
      
      
      for jl in range(kidia, kfdia + 1):
        
        #-----------------------------------------------------------------------
        # Calculate relative humidity limit for snow evaporation
        #-----------------------------------------------------------------------
        zzrh = yrecldp.rprecrhmax + (1.0 - yrecldp.rprecrhmax)*zcovpmax[jl - 1] / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        zzrh = min(max(zzrh, yrecldp.rprecrhmax), 1.0)
        zqe = (zqx[ncldqv - 1, jk - 1, jl - 1] - za[jk - 1, jl - 1]*zqsice[jk - 1, jl - 1]) / max(zepsec, 1.0 - za[jk - 1, jl - 1])
        
        #---------------------------------------------
        # humidity in moistest ZCOVPCLR part of domain
        #---------------------------------------------
        zqe = max(0.0, min(zqe, zqsice[jk - 1, jl - 1]))
        llo1 = zcovpclr[jl - 1] > zepsec and zqx[ncldqs - 1, jk - 1, jl - 1] > zepsec and zqe < zzrh*zqsice[jk - 1, jl - 1]
        
        if llo1:
          
          # Calculate local precipitation (kg/kg)
          zpreclr = zqx[ncldqs - 1, jk - 1, jl - 1] / zcovptot[jl - 1]
          zvpice = foeeice(ztp1[jk - 1, jl - 1])*ydcst.rv / ydcst.rd
          
          # Particle size distribution
          # ZTCG increases Ni with colder temperatures - essentially a
          # Fletcher or Meyers scheme?
          ztcg = 1.0            #v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
          # ZFACX1I modification is based on Andrew Barrett's results
          zfacx1s = 1.0            #v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB
          
          zaplusb = yrecldp.rcl_apb1*zvpice - yrecldp.rcl_apb2*zvpice*ztp1[jk - 1, jl - 1] + pap[jk - 1, jl - 1]*yrecldp.rcl_apb3*ztp1[jk - 1, jl - 1]**3
          zcorrfac = (1.0 / zrho[jl - 1])**0.5
          zcorrfac2 = ((ztp1[jk - 1, jl - 1] / 273.0)**1.5)*(393.0 / (ztp1[jk - 1, jl - 1] + 120.0))
          
          zpr02 = zrho[jl - 1]*zpreclr*yrecldp.rcl_const1s / (ztcg*zfacx1s)
          
          zterm1 = (zqsice[jk - 1, jl - 1] - zqe)*ztp1[jk - 1, jl - 1]**2*zvpice*zcorrfac2*ztcg*yrecldp.rcl_const2s*zfacx1s / (zrho[jl - 1]*zaplusb*zqsice[jk - 1, jl - 1])
          zterm2 = 0.65*yrecldp.rcl_const6s*zpr02**yrecldp.rcl_const4s + yrecldp.rcl_const3s*zcorrfac**0.5*zrho[jl - 1]**0.5*zpr02**yrecldp.rcl_const5s / zcorrfac2**0.5
          
          zdpevap = max(zcovpclr[jl - 1]*zterm1*zterm2*ptsphy, 0.0)
          
          #--------------------------------------------------------------------
          # Limit evaporation to snow amount
          #--------------------------------------------------------------------
          zevap = min(zdpevap, zevaplimice[jl - 1])
          zevap = min(zevap, zqx[ncldqs - 1, jk - 1, jl - 1])
          
          
          zsolqa[ncldqs - 1, ncldqv - 1, jl - 1] = zsolqa[ncldqs - 1, ncldqv - 1, jl - 1] + zevap
          zsolqa[ncldqv - 1, ncldqs - 1, jl - 1] = zsolqa[ncldqv - 1, ncldqs - 1, jl - 1] - zevap
          
          #-------------------------------------------------------------
          # Reduce the total precip coverage proportional to evaporation
          # to mimic the previous scheme which had a diagnostic
          # 2-flux treatment, abandoned due to the new prognostic precip
          #-------------------------------------------------------------
          zcovptot[jl - 1] = max(yrecldp.rcovpmin, zcovptot[jl - 1] - max(0.0, (zcovptot[jl - 1] - za[jk - 1, jl - 1])*zevap / zqx[ncldqs - 1, jk - 1, jl - 1]))
          
          #Update first guess field
          zqxfg[ncldqs - 1, jl - 1] = zqxfg[ncldqs - 1, jl - 1] - zevap
          
      
    # on IEVAPSNOW
    
    #--------------------------------------
    # Evaporate small precipitation amounts
    #--------------------------------------
    for jm in range(1, nclv + 1):
      if llfall[jm - 1]:
        for jl in range(kidia, kfdia + 1):
          if zqxfg[jm - 1, jl - 1] < yrecldp.rlmin:
            zsolqa[jm - 1, ncldqv - 1, jl - 1] = zsolqa[jm - 1, ncldqv - 1, jl - 1] + zqxfg[jm - 1, jl - 1]
            zsolqa[ncldqv - 1, jm - 1, jl - 1] = zsolqa[ncldqv - 1, jm - 1, jl - 1] - zqxfg[jm - 1, jl - 1]
    
    #######################################################################
    #            5.0  *** SOLVERS FOR A AND L ***
    # now use an implicit solution rather than exact solution
    # solver is forward in time, upstream difference for advection
    #######################################################################
    
    #---------------------------
    # 5.1 solver for cloud cover
    #---------------------------
    for jl in range(kidia, kfdia + 1):
      zanew = (za[jk - 1, jl - 1] + zsolac[jl - 1]) / (1.0 + zsolab[jl - 1])
      zanew = min(zanew, 1.0)
      if zanew < yrecldp.ramin:
        zanew = 0.0
      zda[jl - 1] = zanew - zaorig[jk - 1, jl - 1]
      #---------------------------------
      # variables needed for next level
      #---------------------------------
      zanewm1[jl - 1] = zanew
    
    #--------------------------------
    # 5.2 solver for the microphysics
    #--------------------------------
    
    #--------------------------------------------------------------
    # Truncate explicit sinks to avoid negatives
    # Note: Species are treated in the order in which they run out
    # since the clipping will alter the balance for the other vars
    #--------------------------------------------------------------
    
    for jm in range(1, nclv + 1):
      for jn in range(1, nclv + 1):
        for jl in range(kidia, kfdia + 1):
          llindex3[jm - 1, jn - 1, jl - 1] = False
      for jl in range(kidia, kfdia + 1):
        zsinksum[jm - 1, jl - 1] = 0.0
    
    #----------------------------
    # collect sink terms and mark
    #----------------------------
    for jm in range(1, nclv + 1):
      for jn in range(1, nclv + 1):
        for jl in range(kidia, kfdia + 1):
          zsinksum[jm - 1, jl - 1] = zsinksum[jm - 1, jl - 1] - zsolqa[jn - 1, jm - 1, jl - 1]            # +ve total is bad
    
    #---------------------------------------
    # calculate overshoot and scaling factor
    #---------------------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zmax = max(zqx[jm - 1, jk - 1, jl - 1], zepsec)
        zrat = max(zsinksum[jm - 1, jl - 1], zmax)
        zratio[jm - 1, jl - 1] = zmax / zrat
    
    #--------------------------------------------
    # scale the sink terms, in the correct order,
    # recalculating the scale factor each time
    #--------------------------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zsinksum[jm - 1, jl - 1] = 0.0
    
    #----------------
    # recalculate sum
    #----------------
    for jm in range(1, nclv + 1):
      psum_solqa[:] = 0.0
      for jn in range(1, nclv + 1):
        for jl in range(kidia, kfdia + 1):
          psum_solqa[jl - 1] = psum_solqa[jl - 1] + zsolqa[jn - 1, jm - 1, jl - 1]
      for jl in range(kidia, kfdia + 1):
        # ZSINKSUM(JL,JM)=ZSINKSUM(JL,JM)-SUM(ZSOLQA(JL,JM,1:NCLV))
        zsinksum[jm - 1, jl - 1] = zsinksum[jm - 1, jl - 1] - psum_solqa[jl - 1]
      #---------------------------
      # recalculate scaling factor
      #---------------------------
      for jl in range(kidia, kfdia + 1):
        zmm = max(zqx[jm - 1, jk - 1, jl - 1], zepsec)
        zrr = max(zsinksum[jm - 1, jl - 1], zmm)
        zratio[jm - 1, jl - 1] = zmm / zrr
      #------
      # scale
      #------
      for jl in range(kidia, kfdia + 1):
        zzratio = zratio[jm - 1, jl - 1]
        #DIR$ IVDEP
        #DIR$ PREFERVECTOR
        for jn in range(1, nclv + 1):
          if zsolqa[jn - 1, jm - 1, jl - 1] < 0.0:
            zsolqa[jn - 1, jm - 1, jl - 1] = zsolqa[jn - 1, jm - 1, jl - 1]*zzratio
            zsolqa[jm - 1, jn - 1, jl - 1] = zsolqa[jm - 1, jn - 1, jl - 1]*zzratio
    
    #--------------------------------------------------------------
    # 5.2.2 Solver
    #------------------------
    
    #------------------------
    # set the LHS of equation
    #------------------------
    for jm in range(1, nclv + 1):
      for jn in range(1, nclv + 1):
        #----------------------------------------------
        # diagonals: microphysical sink terms+transport
        #----------------------------------------------
        if jn == jm:
          for jl in range(kidia, kfdia + 1):
            zqlhs[jm - 1, jn - 1, jl - 1] = 1.0 + zfallsink[jm - 1, jl - 1]
            for jo in range(1, nclv + 1):
              zqlhs[jm - 1, jn - 1, jl - 1] = zqlhs[jm - 1, jn - 1, jl - 1] + zsolqb[jn - 1, jo - 1, jl - 1]
          #------------------------------------------
          # non-diagonals: microphysical source terms
          #------------------------------------------
        else:
          for jl in range(kidia, kfdia + 1):
            zqlhs[jm - 1, jn - 1, jl - 1] = -zsolqb[jm - 1, jn - 1, jl - 1]              # here is the delta T - missing from doc.
    
    #------------------------
    # set the RHS of equation
    #------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        #---------------------------------
        # sum the explicit source and sink
        #---------------------------------
        zexplicit = 0.0
        for jn in range(1, nclv + 1):
          zexplicit = zexplicit + zsolqa[jn - 1, jm - 1, jl - 1]            # sum over middle index
        zqxn[jm - 1, jl - 1] = zqx[jm - 1, jk - 1, jl - 1] + zexplicit
    
    #-----------------------------------
    # *** solve by LU decomposition: ***
    #-----------------------------------
    
    # Note: This fast way of solving NCLVxNCLV system
    #       assumes a good behaviour (i.e. non-zero diagonal
    #       terms with comparable orders) of the matrix stored
    #       in ZQLHS. For the moment this is the case but
    #       be aware to preserve it when doing eventual
    #       modifications.
    
    # Non pivoting recursive factorization
    for jn in range(1, nclv - 1 + 1):
      # number of steps
      for jm in range(jn + 1, nclv + 1):
        # row index
        for jl in range(kidia, kfdia + 1):
          zqlhs[jn - 1, jm - 1, jl - 1] = zqlhs[jn - 1, jm - 1, jl - 1] / zqlhs[jn - 1, jn - 1, jl - 1]
        for ik in range(jn + 1, nclv + 1):
          # column index
          for jl in range(kidia, kfdia + 1):
            zqlhs[ik - 1, jm - 1, jl - 1] = zqlhs[ik - 1, jm - 1, jl - 1] - zqlhs[jn - 1, jm - 1, jl - 1]*zqlhs[ik - 1, jn - 1, jl - 1]
    
    # Backsubstitution
    #  step 1
    for jn in range(2, nclv + 1):
      for jm in range(1, jn - 1 + 1):
        for jl in range(kidia, kfdia + 1):
          zqxn[jn - 1, jl - 1] = zqxn[jn - 1, jl - 1] - zqlhs[jm - 1, jn - 1, jl - 1]*zqxn[jm - 1, jl - 1]
    #  step 2
    for jl in range(kidia, kfdia + 1):
      zqxn[nclv - 1, jl - 1] = zqxn[nclv - 1, jl - 1] / zqlhs[nclv - 1, nclv - 1, jl - 1]
    for jn in range(nclv - 1, 1 + -1, -1):
      for jm in range(jn + 1, nclv + 1):
        for jl in range(kidia, kfdia + 1):
          zqxn[jn - 1, jl - 1] = zqxn[jn - 1, jl - 1] - zqlhs[jm - 1, jn - 1, jl - 1]*zqxn[jm - 1, jl - 1]
      for jl in range(kidia, kfdia + 1):
        zqxn[jn - 1, jl - 1] = zqxn[jn - 1, jl - 1] / zqlhs[jn - 1, jn - 1, jl - 1]
    
    # Ensure no small values (including negatives) remain in cloud variables nor
    # precipitation rates.
    # Evaporate l,i,r,s to water vapour. Latent heating taken into account below
    for jn in range(1, nclv - 1 + 1):
      for jl in range(kidia, kfdia + 1):
        if zqxn[jn - 1, jl - 1] < zepsec:
          zqxn[ncldqv - 1, jl - 1] = zqxn[ncldqv - 1, jl - 1] + zqxn[jn - 1, jl - 1]
          zqxn[jn - 1, jl - 1] = 0.0
    
    #--------------------------------
    # variables needed for next level
    #--------------------------------
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zqxnm1[jm - 1, jl - 1] = zqxn[jm - 1, jl - 1]
        zqxn2d[jm - 1, jk - 1, jl - 1] = zqxn[jm - 1, jl - 1]
    
    #------------------------------------------------------------------------
    # 5.3 Precipitation/sedimentation fluxes to next level
    #     diagnostic precipitation fluxes
    #     It is this scaled flux that must be used for source to next layer
    #------------------------------------------------------------------------
    
    for jm in range(1, nclv + 1):
      for jl in range(kidia, kfdia + 1):
        zpfplsx[jm - 1, jk + 1 - 1, jl - 1] = zfallsink[jm - 1, jl - 1]*zqxn[jm - 1, jl - 1]*zrdtgdp[jl - 1]
    
    # Ensure precipitation fraction is zero if no precipitation
    for jl in range(kidia, kfdia + 1):
      zqpretot[jl - 1] = zpfplsx[ncldqs - 1, jk + 1 - 1, jl - 1] + zpfplsx[ncldqr - 1, jk + 1 - 1, jl - 1]
    for jl in range(kidia, kfdia + 1):
      if zqpretot[jl - 1] < zepsec:
        zcovptot[jl - 1] = 0.0
    
    #######################################################################
    #              6  *** UPDATE TENDANCIES ***
    #######################################################################
    
    #--------------------------------
    # 6.1 Temperature and CLV budgets
    #--------------------------------
    
    for jm in range(1, nclv - 1 + 1):
      for jl in range(kidia, kfdia + 1):
        
        # calculate fluxes in and out of box for conservation of TL
        zfluxq[jm - 1, jl - 1] = zpsupsatsrce[jm - 1, jl - 1] + zconvsrce[jm - 1, jl - 1] + zfallsrce[jm - 1, jl - 1] - (zfallsink[jm - 1, jl - 1] + zconvsink[jm - 1, jl - 1])*zqxn[jm - 1, jl - 1]
      
      if iphase[jm - 1] == 1:
        for jl in range(kidia, kfdia + 1):
          tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] + ydthf.ralvdcp*(zqxn[jm - 1, jl - 1] - zqx[jm - 1, jk - 1, jl - 1] - zfluxq[jm - 1, jl - 1])*zqtmst
      
      if iphase[jm - 1] == 2:
        for jl in range(kidia, kfdia + 1):
          tendency_loc_t[jk - 1, jl - 1] = tendency_loc_t[jk - 1, jl - 1] + ydthf.ralsdcp*(zqxn[jm - 1, jl - 1] - zqx[jm - 1, jk - 1, jl - 1] - zfluxq[jm - 1, jl - 1])*zqtmst
      
      #----------------------------------------------------------------------
      # New prognostic tendencies - ice,liquid rain,snow
      # Note: CLV arrays use PCLV in calculation of tendency while humidity
      #       uses ZQX. This is due to clipping at start of cloudsc which
      #       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
      #----------------------------------------------------------------------
      for jl in range(kidia, kfdia + 1):
        tendency_loc_cld[jm - 1, jk - 1, jl - 1] = tendency_loc_cld[jm - 1, jk - 1, jl - 1] + (zqxn[jm - 1, jl - 1] - zqx0[jm - 1, jk - 1, jl - 1])*zqtmst
      
    
    for jl in range(kidia, kfdia + 1):
      #----------------------
      # 6.2 Humidity budget
      #----------------------
      tendency_loc_q[jk - 1, jl - 1] = tendency_loc_q[jk - 1, jl - 1] + (zqxn[ncldqv - 1, jl - 1] - zqx[ncldqv - 1, jk - 1, jl - 1])*zqtmst
      
      #-------------------
      # 6.3 cloud cover
      #-----------------------
      tendency_loc_a[jk - 1, jl - 1] = tendency_loc_a[jk - 1, jl - 1] + zda[jl - 1]*zqtmst
    
    #--------------------------------------------------
    # Copy precipitation fraction into output variable
    #-------------------------------------------------
    for jl in range(kidia, kfdia + 1):
      pcovptot[jk - 1, jl - 1] = zcovptot[jl - 1]
    
  # on vertical level JK
  #----------------------------------------------------------------------
  #                       END OF VERTICAL LOOP
  #----------------------------------------------------------------------
  
  #######################################################################
  #              8  *** FLUX/DIAGNOSTICS COMPUTATIONS ***
  #######################################################################
  
  #--------------------------------------------------------------------
  # Copy general precip arrays back into PFP arrays for GRIB archiving
  # Add rain and liquid fluxes, ice and snow fluxes
  #--------------------------------------------------------------------
  for jk in range(1, klev + 1 + 1):
    for jl in range(kidia, kfdia + 1):
      pfplsl[jk - 1, jl - 1] = zpfplsx[ncldqr - 1, jk - 1, jl - 1] + zpfplsx[ncldql - 1, jk - 1, jl - 1]
      pfplsn[jk - 1, jl - 1] = zpfplsx[ncldqs - 1, jk - 1, jl - 1] + zpfplsx[ncldqi - 1, jk - 1, jl - 1]
  
  #--------
  # Fluxes:
  #--------
  for jl in range(kidia, kfdia + 1):
    pfsqlf[1 - 1, jl - 1] = 0.0
    pfsqif[1 - 1, jl - 1] = 0.0
    pfsqrf[1 - 1, jl - 1] = 0.0
    pfsqsf[1 - 1, jl - 1] = 0.0
    pfcqlng[1 - 1, jl - 1] = 0.0
    pfcqnng[1 - 1, jl - 1] = 0.0
    pfcqrng[1 - 1, jl - 1] = 0.0      #rain
    pfcqsng[1 - 1, jl - 1] = 0.0      #snow
    # fluxes due to turbulence
    pfsqltur[1 - 1, jl - 1] = 0.0
    pfsqitur[1 - 1, jl - 1] = 0.0
  
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      
      zgdph_r = -zrg_r*(paph[jk + 1 - 1, jl - 1] - paph[jk - 1, jl - 1])*zqtmst
      pfsqlf[jk + 1 - 1, jl - 1] = pfsqlf[jk - 1, jl - 1]
      pfsqif[jk + 1 - 1, jl - 1] = pfsqif[jk - 1, jl - 1]
      pfsqrf[jk + 1 - 1, jl - 1] = pfsqlf[jk - 1, jl - 1]
      pfsqsf[jk + 1 - 1, jl - 1] = pfsqif[jk - 1, jl - 1]
      pfcqlng[jk + 1 - 1, jl - 1] = pfcqlng[jk - 1, jl - 1]
      pfcqnng[jk + 1 - 1, jl - 1] = pfcqnng[jk - 1, jl - 1]
      pfcqrng[jk + 1 - 1, jl - 1] = pfcqlng[jk - 1, jl - 1]
      pfcqsng[jk + 1 - 1, jl - 1] = pfcqnng[jk - 1, jl - 1]
      pfsqltur[jk + 1 - 1, jl - 1] = pfsqltur[jk - 1, jl - 1]
      pfsqitur[jk + 1 - 1, jl - 1] = pfsqitur[jk - 1, jl - 1]
      
      zalfaw = zfoealfa[jk - 1, jl - 1]
      
      # Liquid , LS scheme minus detrainment
      pfsqlf[jk + 1 - 1, jl - 1] = pfsqlf[jk + 1 - 1, jl - 1] + (zqxn2d[ncldql - 1, jk - 1, jl - 1] - zqx0[ncldql - 1, jk - 1, jl - 1] + pvfl[jk - 1, jl - 1]*ptsphy - zalfaw*plude[jk - 1, jl - 1])*zgdph_r
      # liquid, negative numbers
      pfcqlng[jk + 1 - 1, jl - 1] = pfcqlng[jk + 1 - 1, jl - 1] + zlneg[ncldql - 1, jk - 1, jl - 1]*zgdph_r
      
      # liquid, vertical diffusion
      pfsqltur[jk + 1 - 1, jl - 1] = pfsqltur[jk + 1 - 1, jl - 1] + pvfl[jk - 1, jl - 1]*ptsphy*zgdph_r
      
      # Rain, LS scheme
      pfsqrf[jk + 1 - 1, jl - 1] = pfsqrf[jk + 1 - 1, jl - 1] + (zqxn2d[ncldqr - 1, jk - 1, jl - 1] - zqx0[ncldqr - 1, jk - 1, jl - 1])*zgdph_r
      # rain, negative numbers
      pfcqrng[jk + 1 - 1, jl - 1] = pfcqrng[jk + 1 - 1, jl - 1] + zlneg[ncldqr - 1, jk - 1, jl - 1]*zgdph_r
      
      # Ice , LS scheme minus detrainment
      pfsqif[jk + 1 - 1, jl - 1] = pfsqif[jk + 1 - 1, jl - 1] + (zqxn2d[ncldqi - 1, jk - 1, jl - 1] - zqx0[ncldqi - 1, jk - 1, jl - 1] + pvfi[jk - 1, jl - 1]*ptsphy - (1.0 - zalfaw)*plude[jk - 1, jl - 1])*zgdph_r
      # ice, negative numbers
      pfcqnng[jk + 1 - 1, jl - 1] = pfcqnng[jk + 1 - 1, jl - 1] + zlneg[ncldqi - 1, jk - 1, jl - 1]*zgdph_r
      
      # ice, vertical diffusion
      pfsqitur[jk + 1 - 1, jl - 1] = pfsqitur[jk + 1 - 1, jl - 1] + pvfi[jk - 1, jl - 1]*ptsphy*zgdph_r
      
      # snow, LS scheme
      pfsqsf[jk + 1 - 1, jl - 1] = pfsqsf[jk + 1 - 1, jl - 1] + (zqxn2d[ncldqs - 1, jk - 1, jl - 1] - zqx0[ncldqs - 1, jk - 1, jl - 1])*zgdph_r
      # snow, negative numbers
      pfcqsng[jk + 1 - 1, jl - 1] = pfcqsng[jk + 1 - 1, jl - 1] + zlneg[ncldqs - 1, jk - 1, jl - 1]*zgdph_r
  
  #-----------------------------------
  # enthalpy flux due to precipitation
  #-----------------------------------
  for jk in range(1, klev + 1 + 1):
    for jl in range(kidia, kfdia + 1):
      pfhpsl[jk - 1, jl - 1] = -ydcst.rlvtt*pfplsl[jk - 1, jl - 1]
      pfhpsn[jk - 1, jl - 1] = -ydcst.rlstt*pfplsn[jk - 1, jl - 1]
  
  #===============================================================================
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)
  return 
