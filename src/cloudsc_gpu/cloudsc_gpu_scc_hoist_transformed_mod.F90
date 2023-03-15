! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_GPU_SCC_HOIST_MOD

  CONTAINS
  SUBROUTINE CLOUDSC_SCC_HOIST (KIDIA, KFDIA, KLON, KLEV, PTSPHY, PT, PQ, TENDENCY_TMP_T, TENDENCY_TMP_Q, TENDENCY_TMP_A,  &
  & TENDENCY_TMP_CLD, TENDENCY_LOC_T, TENDENCY_LOC_Q, TENDENCY_LOC_A, TENDENCY_LOC_CLD, PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI,  &
  & PHRSW, PHRLW, PVERVEL, PAP, PAPH, PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER,  &
  & PICRIT_AER, PRE_ICE, PCCN, PNICE, PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG,  &
  & PFCQSNG, PFSQLTUR, PFSQITUR, PFPLSL, PFPLSN, PFHPSL, PFHPSN, YDCST, YDTHF, YDECLDP, ZFOEALFA, ZTP1, ZLI, ZA, ZAORIG,  &
  & ZLIQFRAC, ZICEFRAC, ZQX, ZQX0, ZPFPLSX, ZLNEG, ZQXN2D, ZQSMIX, ZQSLIQ, ZQSICE, ZFOEEWMT, ZFOEEW, ZFOEELIQT, JL)
    !---input
    !---prognostic fields
    !-- arrays for aerosol-cloud interactions
    !!! & PQAER,    KAER, &
    !---diagnostic output
    !---resulting fluxes

    !===============================================================================
    !**** *CLOUDSC* -  ROUTINE FOR PARAMATERIZATION OF CLOUD PROCESSES
    !                  FOR PROGNOSTIC CLOUD SCHEME
    !!
    !     M.Tiedtke, C.Jakob, A.Tompkins, R.Forbes     (E.C.M.W.F.)
    !!
    !     PURPOSE
    !     -------
    !          THIS ROUTINE UPDATES THE CONV/STRAT CLOUD FIELDS.
    !          THE FOLLOWING PROCESSES ARE CONSIDERED:
    !        - Detrainment of cloud water from convective updrafts
    !        - Evaporation/condensation of cloud water in connection
    !           with heating/cooling such as by subsidence/ascent
    !        - Erosion of clouds by turbulent mixing of cloud air
    !           with unsaturated environmental air
    !        - Deposition onto ice when liquid water present (Bergeron-Findeison)
    !        - Conversion of cloud water into rain (collision-coalescence)
    !        - Conversion of cloud ice to snow (aggregation)
    !        - Sedimentation of rain, snow and ice
    !        - Evaporation of rain and snow
    !        - Melting of snow and ice
    !        - Freezing of liquid and rain
    !        Note: Turbulent transports of s,q,u,v at cloud tops due to
    !           buoyancy fluxes and lw radiative cooling are treated in
    !           the VDF scheme
    !!
    !     INTERFACE.
    !     ----------
    !          *CLOUDSC* IS CALLED FROM *CALLPAR*
    !     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
    !     T,Q,L,PHI AND DETRAINMENT OF CLOUD WATER FROM THE
    !     CONVECTIVE CLOUDS (MASSFLUX CONVECTION SCHEME), BOUNDARY
    !     LAYER TURBULENT FLUXES OF HEAT AND MOISTURE, RADIATIVE FLUXES,
    !     OMEGA.
    !     IT RETURNS ITS OUTPUT TO:
    !      1.MODIFIED TENDENCIES OF MODEL VARIABLES T AND Q
    !        AS WELL AS CLOUD VARIABLES L AND C
    !      2.GENERATES PRECIPITATION FLUXES FROM STRATIFORM CLOUDS
    !!
    !     EXTERNALS.
    !     ----------
    !          NONE
    !!
    !     MODIFICATIONS.
    !     -------------
    !      M. TIEDTKE    E.C.M.W.F.     8/1988, 2/1990
    !     CH. JAKOB      E.C.M.W.F.     2/1994 IMPLEMENTATION INTO IFS
    !     A.TOMPKINS     E.C.M.W.F.     2002   NEW NUMERICS
    !        01-05-22 : D.Salmond   Safety modifications
    !        02-05-29 : D.Salmond   Optimisation
    !        03-01-13 : J.Hague     MASS Vector Functions  J.Hague
    !        03-10-01 : M.Hamrud    Cleaning
    !        04-12-14 : A.Tompkins  New implicit solver and physics changes
    !        04-12-03 : A.Tompkins & M.Ko"hler  moist PBL
    !     G.Mozdzynski  09-Jan-2006  EXP security fix
    !        19-01-09 : P.Bechtold  Changed increased RCLDIFF value for KTYPE=2
    !        07-07-10 : A.Tompkins/R.Forbes  4-Phase flexible microphysics
    !        01-03-11 : R.Forbes    Mixed phase changes and tidy up
    !        01-10-11 : R.Forbes    Melt ice to rain, allow rain to freeze
    !        01-10-11 : R.Forbes    Limit supersat to avoid excessive values
    !        31-10-11 : M.Ahlgrimm  Add rain, snow and PEXTRA to DDH output
    !        17-02-12 : F.Vana      Simplified/optimized LU factorization
    !        18-05-12 : F.Vana      Cleaning + better support of sequential physics
    !        N.Semane+P.Bechtold     04-10-2012 Add RVRFACTOR factor for small planet
    !        01-02-13 : R.Forbes    New params of autoconv/acc,rain evap,snow riming
    !        15-03-13 : F. Vana     New dataflow + more tendencies from the first call
    !        K. Yessad (July 2014): Move some variables.
    !        F. Vana  05-Mar-2015  Support for single precision
    !        15-01-15 : R.Forbes    Added new options for snow evap & ice deposition
    !        10-01-15 : R.Forbes    New physics for rain freezing
    !        23-10-14 : P. Bechtold remove zeroing of convection arrays
    !
    !     SWITCHES.
    !     --------
    !!
    !     MODEL PARAMETERS
    !     ----------------
    !     RCLDIFF:    PARAMETER FOR EROSION OF CLOUDS
    !     RCLCRIT_SEA:  THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER SEA
    !     RCLCRIT_LAND: THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER LAND
    !     RLCRITSNOW: THRESHOLD VALUE FOR SNOW AUTOCONVERSION
    !     RKCONV:     PARAMETER FOR AUTOCONVERSION OF CLOUDS (KESSLER)
    !     RCLDMAX:    MAXIMUM POSSIBLE CLW CONTENT (MASON,1971)
    !!
    !     REFERENCES.
    !     ----------
    !     TIEDTKE MWR 1993
    !     JAKOB PhD 2000
    !     GREGORY ET AL. QJRMS 2000
    !     TOMPKINS ET AL. QJRMS 2007
    !!
    !===============================================================================

    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMPHYDER, ONLY: state_type
    USE YOECLDP, ONLY: NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
    USE YOECLDP, ONLY: TECLDP
    USE YOEPHLI, ONLY: TEPHLI
    USE YOMCST, ONLY: TOMCST
    USE YOETHF, ONLY: TOETHF





    IMPLICIT NONE

    !-------------------------------------------------------------------------------
    !                 Declare input/output arguments
    !-------------------------------------------------------------------------------

    ! PLCRIT_AER : critical liquid mmr for rain autoconversion process
    ! PICRIT_AER : critical liquid mmr for snow autoconversion process
    ! PRE_LIQ : liq Re
    ! PRE_ICE : ice Re
    ! PCCN    : liquid cloud condensation nuclei
    ! PNICE   : ice number concentration (cf. CCN)

    REAL(KIND=JPRB), INTENT(IN) :: PLCRIT_AER(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PICRIT_AER(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PRE_ICE(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PCCN(KLON, KLEV)    ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), INTENT(IN) :: PNICE(KLON, KLEV)
    ! ice number concentration (cf. CCN)

    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! Number of grid points
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV    ! Number of levels
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
    REAL(KIND=JPRB), INTENT(IN) :: PTSPHY    ! Physics timestep
    REAL(KIND=JPRB), INTENT(IN) :: PT(KLON, KLEV)    ! T at start of callpar
    REAL(KIND=JPRB), INTENT(IN) :: PQ(KLON, KLEV)    ! Q at start of callpar
    REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_T(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_Q(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_A(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_CLD(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_T(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_Q(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_A(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_CLD(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: PVFA(KLON, KLEV)    ! CC from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PVFL(KLON, KLEV)    ! Liq from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PVFI(KLON, KLEV)    ! Ice from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PDYNA(KLON, KLEV)    ! CC from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PDYNL(KLON, KLEV)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PDYNI(KLON, KLEV)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PHRSW(KLON, KLEV)    ! Short-wave heating rate
    REAL(KIND=JPRB), INTENT(IN) :: PHRLW(KLON, KLEV)    ! Long-wave heating rate
    REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(KLON, KLEV)    !Vertical velocity
    REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)    ! Pressure on full levels
    REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)    ! Pressure on half levels
    REAL(KIND=JPRB), INTENT(IN) :: PLSM(KLON)    ! Land fraction (0-1)
    LOGICAL, INTENT(IN) :: LDCUM(KLON)    ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(KLON)    ! Convection type 0,1,2
    REAL(KIND=JPRB), INTENT(IN) :: PLU(KLON, KLEV)    ! Conv. condensate
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(KLON, KLEV)    ! Conv. detrained water
    REAL(KIND=JPRB), INTENT(IN) :: PSNDE(KLON, KLEV)    ! Conv. detrained snow
    REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)    ! Conv. mass flux up
    REAL(KIND=JPRB), INTENT(IN) :: PMFD(KLON, KLEV)    ! Conv. mass flux down
    REAL(KIND=JPRB), INTENT(IN) :: PA(KLON, KLEV)
    ! Original Cloud fraction (t)

    REAL(KIND=JPRB), INTENT(IN) :: PCLV(KLON, KLEV, NCLV)

    ! Supersat clipped at previous time level in SLTEND
    REAL(KIND=JPRB), INTENT(IN) :: PSUPSAT(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: PCOVPTOT(KLON, KLEV)    ! Precip fraction
    REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(KLON)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(KLON, KLEV + 1)    ! Flux of liquid
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(KLON, KLEV + 1)    ! Flux of ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(KLON, KLEV + 1)    ! -ve corr for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(KLON, KLEV + 1)    ! -ve corr for ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(KLON, KLEV + 1)    ! Flux diagnostics
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(KLON, KLEV + 1)    !    for DDH, generic
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(KLON, KLEV + 1)    ! rain
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(KLON, KLEV + 1)    ! snow
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(KLON, KLEV + 1)    ! liquid flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(KLON, KLEV + 1)    ! ice flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(KLON, KLEV + 1)    ! liq+rain sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(KLON, KLEV + 1)    ! ice+snow sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(KLON, KLEV + 1)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(KLON, KLEV + 1)
    ! Enthalp flux for ice

    TYPE(TOMCST), INTENT(IN) :: YDCST
    TYPE(TOETHF), INTENT(IN) :: YDTHF
    TYPE(TECLDP), INTENT(IN) :: YDECLDP

    !-------------------------------------------------------------------------------
    !                       Declare local variables
    !-------------------------------------------------------------------------------

    REAL(KIND=JPRB) :: ZLCOND1, ZLCOND2, ZLEVAP, ZLEROS, ZLEVAPL, ZLEVAPI, ZRAINAUT, ZSNOWAUT, ZLIQCLD, ZICECLD(KLEV)
    !  condensation and evaporation terms
    ! autoconversion terms
    REAL(KIND=JPRB) :: ZFOKOOP(KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZFOEALFA(KLON, KLEV + 1)
    REAL(KIND=JPRB) :: ZICENUCLEI
    ! number concentration of ice nuclei

    REAL(KIND=JPRB) :: ZLICLD(KLEV)
    REAL(KIND=JPRB) :: ZACOND
    REAL(KIND=JPRB) :: ZAEROS
    REAL(KIND=JPRB) :: ZLFINALSUM(KLEV)
    REAL(KIND=JPRB) :: ZDQS
    REAL(KIND=JPRB) :: ZTOLD
    REAL(KIND=JPRB) :: ZQOLD
    REAL(KIND=JPRB) :: ZDTGDP(KLEV)
    REAL(KIND=JPRB) :: ZRDTGDP(KLEV)
    REAL(KIND=JPRB) :: ZTRPAUS
    REAL(KIND=JPRB) :: ZCOVPCLR
    REAL(KIND=JPRB) :: ZPRECLR
    REAL(KIND=JPRB) :: ZCOVPTOT
    REAL(KIND=JPRB) :: ZCOVPMAX
    REAL(KIND=JPRB) :: ZQPRETOT(KLEV)
    REAL(KIND=JPRB) :: ZDPEVAP
    REAL(KIND=JPRB) :: ZDTFORC
    REAL(KIND=JPRB) :: ZDTDIAB
    REAL(KIND=JPRB), INTENT(INOUT) :: ZTP1(KLON, KLEV)
    REAL(KIND=JPRB) :: ZLDEFR(KLEV)
    REAL(KIND=JPRB) :: ZLDIFDT
    REAL(KIND=JPRB) :: ZDTGDPF
    REAL(KIND=JPRB) :: ZLCUST(NCLV)
    REAL(KIND=JPRB) :: ZACUST
    REAL(KIND=JPRB) :: ZMF

    REAL(KIND=JPRB) :: ZRHO(KLEV)
    REAL(KIND=JPRB) :: ZTMP1, ZTMP2, ZTMP3
    REAL(KIND=JPRB) :: ZTMP4, ZTMP5, ZTMP6, ZTMP7
    REAL(KIND=JPRB) :: ZALFAWM

    ! Accumulators of A,B,and C factors for cloud equations
    REAL(KIND=JPRB) :: ZSOLAB(KLEV)    ! -ve implicit CC
    REAL(KIND=JPRB) :: ZSOLAC(KLEV)    ! linear CC
    REAL(KIND=JPRB) :: ZANEW
    REAL(KIND=JPRB) :: ZANEWM1

    REAL(KIND=JPRB) :: ZGDP

    !---for flux calculation
    REAL(KIND=JPRB) :: ZDA(KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZLI(KLON, KLEV), ZA(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZAORIG(KLON, KLEV)
    ! start of scheme value for CC

    LOGICAL :: LLFLAG
    LOGICAL :: LLO1

    INTEGER(KIND=JPIM) :: ICALL, IK, JK, JL, JM, JN, JO, JLEN, IS

    REAL(KIND=JPRB) :: ZDP(KLEV), ZPAPHD

    REAL(KIND=JPRB) :: ZALFA
    ! & ZALFACU, ZALFALS
    REAL(KIND=JPRB) :: ZALFAW
    REAL(KIND=JPRB) :: ZBETA, ZBETA1
    !REAL(KIND=JPRB) :: ZBOTT
    REAL(KIND=JPRB) :: ZCFPR
    REAL(KIND=JPRB) :: ZCOR
    REAL(KIND=JPRB) :: ZCDMAX
    REAL(KIND=JPRB) :: ZMIN
    REAL(KIND=JPRB) :: ZLCONDLIM
    REAL(KIND=JPRB) :: ZDENOM
    REAL(KIND=JPRB) :: ZDPMXDT
    REAL(KIND=JPRB) :: ZDPR
    REAL(KIND=JPRB) :: ZDTDP
    REAL(KIND=JPRB) :: ZE
    REAL(KIND=JPRB) :: ZEPSEC
    REAL(KIND=JPRB) :: ZFAC, ZFACI, ZFACW
    REAL(KIND=JPRB) :: ZGDCP
    REAL(KIND=JPRB) :: ZINEW
    REAL(KIND=JPRB) :: ZLCRIT
    REAL(KIND=JPRB) :: ZMFDN
    REAL(KIND=JPRB) :: ZPRECIP
    REAL(KIND=JPRB) :: ZQE
    REAL(KIND=JPRB) :: ZQSAT, ZQTMST, ZRDCP
    REAL(KIND=JPRB) :: ZRHC, ZSIG, ZSIGK
    REAL(KIND=JPRB) :: ZWTOT
    REAL(KIND=JPRB) :: ZZCO, ZZDL, ZZRH, ZZZDT, ZQADJ
    REAL(KIND=JPRB) :: ZQNEW, ZTNEW
    REAL(KIND=JPRB) :: ZRG_R, ZGDPH_R, ZCONS1, ZCOND, ZCONS1A
    REAL(KIND=JPRB) :: ZLFINAL
    REAL(KIND=JPRB) :: ZMELT
    REAL(KIND=JPRB) :: ZEVAP
    REAL(KIND=JPRB) :: ZFRZ
    REAL(KIND=JPRB) :: ZVPLIQ, ZVPICE
    REAL(KIND=JPRB) :: ZADD, ZBDD, ZCVDS, ZICE0, ZDEPOS
    REAL(KIND=JPRB) :: ZSUPSAT
    REAL(KIND=JPRB) :: ZFALL
    REAL(KIND=JPRB) :: ZRE_ICE
    REAL(KIND=JPRB) :: ZRLDCP
    REAL(KIND=JPRB) :: ZQP1ENV

    !----------------------------
    ! Arrays for new microphysics
    !----------------------------
    INTEGER(KIND=JPIM) :: IPHASE(NCLV)
    ! marker for water phase of each species
    ! 0=vapour, 1=liquid, 2=ice

    INTEGER(KIND=JPIM) :: IMELT(NCLV)
    ! marks melting linkage for ice categories
    ! ice->liquid, snow->rain

    LOGICAL :: LLFALL(NCLV)
    ! marks falling species
    ! LLFALL=0, cloud cover must > 0 for zqx > 0
    ! LLFALL=1, no cloud needed, zqx can evaporate

    LOGICAL :: LLINDEX1(NCLV)    ! index variable
    LOGICAL :: LLINDEX3(NCLV, NCLV)    ! index variable
    REAL(KIND=JPRB) :: ZMAX
    REAL(KIND=JPRB) :: ZRAT
    INTEGER(KIND=JPIM) :: IORDER(NCLV)
    ! array for sorting explicit terms

    REAL(KIND=JPRB), INTENT(INOUT) :: ZLIQFRAC(KLON, KLEV)    ! cloud liquid water fraction: ql/(ql+qi)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZICEFRAC(KLON, KLEV)    ! cloud ice water fraction: qi/(ql+qi)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQX(KLON, KLEV, NCLV)    ! water variables
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQX0(KLON, KLEV, NCLV)    ! water variables at start of scheme
    REAL(KIND=JPRB) :: ZQXN(NCLV, KLEV)    ! new values for zqx at time+1
    REAL(KIND=JPRB) :: ZQXFG(NCLV, KLEV)    ! first guess values including precip
    REAL(KIND=JPRB) :: ZQXNM1(NCLV)    ! new values for zqx at time+1 at level above
    REAL(KIND=JPRB) :: ZFLUXQ(NCLV)
    ! fluxes convergence of species (needed?)
    ! Keep the following for possible future total water variance scheme?
    !REAL(KIND=JPRB) :: ZTL(KLON,KLEV)       ! liquid water temperature
    !REAL(KIND=JPRB) :: ZABETA(KLON,KLEV)    ! cloud fraction
    !REAL(KIND=JPRB) :: ZVAR(KLON,KLEV)      ! temporary variance
    !REAL(KIND=JPRB) :: ZQTMIN(KLON,KLEV)
    !REAL(KIND=JPRB) :: ZQTMAX(KLON,KLEV)

    REAL(KIND=JPRB), INTENT(INOUT) :: ZPFPLSX(KLON, KLEV + 1, NCLV)    ! generalized precipitation flux
    REAL(KIND=JPRB), INTENT(INOUT) :: ZLNEG(KLON, KLEV, NCLV)    ! for negative correction diagnostics
    REAL(KIND=JPRB) :: ZMELTMAX
    REAL(KIND=JPRB) :: ZFRZMAX
    REAL(KIND=JPRB) :: ZICETOT

    REAL(KIND=JPRB), INTENT(INOUT) :: ZQXN2D(KLON, KLEV, NCLV)
    ! water variables store

    REAL(KIND=JPRB), INTENT(INOUT) :: ZQSMIX(KLON, KLEV)
    ! diagnostic mixed phase saturation
    !REAL(KIND=JPRB) :: ZQSBIN(KLON,KLEV) ! binary switched ice/liq saturation
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQSLIQ(KLON, KLEV)    ! liquid water saturation
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQSICE(KLON, KLEV)
    ! ice water saturation

    !REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
    !REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
    !REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice

    REAL(KIND=JPRB), INTENT(INOUT) :: ZFOEEWMT(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZFOEEW(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZFOEELIQT(KLON, KLEV)
    !REAL(KIND=JPRB) :: ZFOEEICET(KLON,KLEV)

    REAL(KIND=JPRB) :: ZDQSLIQDT, ZDQSICEDT, ZDQSMIXDT(KLEV)
    REAL(KIND=JPRB) :: ZCORQSLIQ(KLEV)
    REAL(KIND=JPRB) :: ZCORQSICE(KLEV)
    !REAL(KIND=JPRB) :: ZCORQSBIN(KLON)
    REAL(KIND=JPRB) :: ZCORQSMIX
    REAL(KIND=JPRB) :: ZEVAPLIMLIQ, ZEVAPLIMICE(KLEV), ZEVAPLIMMIX(KLEV)

    !-------------------------------------------------------
    ! SOURCE/SINK array for implicit and explicit terms
    !-------------------------------------------------------
    ! a POSITIVE value entered into the arrays is a...
    !            Source of this variable
    !            |
    !            |   Sink of this variable
    !            |   |
    !            V   V
    ! ZSOLQA(JL,IQa,IQb)  = explicit terms
    ! ZSOLQB(JL,IQa,IQb)  = implicit terms
    ! Thus if ZSOLAB(JL,NCLDQL,IQV)=K where K>0 then this is
    ! a source of NCLDQL and a sink of IQV
    ! put 'magic' source terms such as PLUDE from
    ! detrainment into explicit source/sink array diagnognal
    ! ZSOLQA(NCLDQL,NCLDQL)= -PLUDE
    ! i.e. A positive value is a sink!????? weird...
    !-------------------------------------------------------

    REAL(KIND=JPRB) :: ZSOLQA(NCLV, NCLV, KLEV)    ! explicit sources and sinks
    REAL(KIND=JPRB) :: ZSOLQB(NCLV, NCLV, KLEV)
    ! implicit sources and sinks
    ! e.g. microphysical pathways between ice variables.
    REAL(KIND=JPRB) :: ZQLHS(NCLV, NCLV)    ! n x n matrix storing the LHS of implicit solver
    REAL(KIND=JPRB) :: ZVQX(NCLV)    ! fall speeds of three categories
    REAL(KIND=JPRB) :: ZEXPLICIT, ZRATIO(NCLV), ZSINKSUM(NCLV)

    ! for sedimentation source/sink terms
    REAL(KIND=JPRB) :: ZFALLSINK(NCLV, KLEV)
    REAL(KIND=JPRB) :: ZFALLSRCE(NCLV, KLEV)

    ! for convection detrainment source and subsidence source/sink terms
    REAL(KIND=JPRB) :: ZCONVSRCE(NCLV, KLEV)
    REAL(KIND=JPRB) :: ZCONVSINK(NCLV, KLEV)

    ! for supersaturation source term from previous timestep
    REAL(KIND=JPRB) :: ZPSUPSATSRCE(NCLV, KLEV)

    ! Numerical fit to wet bulb temperature
    REAL(KIND=JPRB), PARAMETER :: ZTW1 = 1329.31_JPRB
    REAL(KIND=JPRB), PARAMETER :: ZTW2 = 0.0074615_JPRB
    REAL(KIND=JPRB), PARAMETER :: ZTW3 = 0.85E5_JPRB
    REAL(KIND=JPRB), PARAMETER :: ZTW4 = 40.637_JPRB
    REAL(KIND=JPRB), PARAMETER :: ZTW5 = 275.0_JPRB

    REAL(KIND=JPRB) :: ZSUBSAT    ! Subsaturation for snow melting term
    REAL(KIND=JPRB) :: ZTDMTW0
    ! Diff between dry-bulb temperature and
    ! temperature when wet-bulb = 0degC

    ! Variables for deposition term
    REAL(KIND=JPRB) :: ZTCG    ! Temperature dependent function for ice PSD
    REAL(KIND=JPRB) :: ZFACX1I, ZFACX1S    ! PSD correction factor
    REAL(KIND=JPRB) :: ZAPLUSB, ZCORRFAC, ZCORRFAC2, ZPR02, ZTERM1, ZTERM2    ! for ice dep
    REAL(KIND=JPRB) :: ZCLDTOPDIST    ! Distance from cloud top
    REAL(KIND=JPRB) :: ZINFACTOR
    ! No. of ice nuclei factor for deposition

    ! Autoconversion/accretion/riming/evaporation
    INTEGER(KIND=JPIM) :: IWARMRAIN
    INTEGER(KIND=JPIM) :: IEVAPRAIN
    INTEGER(KIND=JPIM) :: IEVAPSNOW
    INTEGER(KIND=JPIM) :: IDEPICE
    REAL(KIND=JPRB) :: ZRAINACC
    REAL(KIND=JPRB) :: ZRAINCLD
    REAL(KIND=JPRB) :: ZSNOWRIME
    REAL(KIND=JPRB) :: ZSNOWCLD
    REAL(KIND=JPRB) :: ZESATLIQ
    REAL(KIND=JPRB) :: ZFALLCORR
    REAL(KIND=JPRB) :: ZLAMBDA
    REAL(KIND=JPRB) :: ZEVAP_DENOM
    REAL(KIND=JPRB) :: ZCORR2
    REAL(KIND=JPRB) :: ZKA
    REAL(KIND=JPRB) :: ZCONST
    REAL(KIND=JPRB) :: ZTEMP

    ! Rain freezing
    LOGICAL :: LLRAINLIQ
    ! True if majority of raindrops are liquid (no ice core)

    !----------------------------
    ! End: new microphysics
    !----------------------------

    !----------------------
    ! SCM budget statistics
    !----------------------
    REAL(KIND=JPRB) :: ZRAIN

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    REAL(KIND=JPRB) :: ZTMPL, ZTMPI, ZTMPA

    REAL(KIND=JPRB) :: ZMM, ZRR
    REAL(KIND=JPRB) :: ZRG

    REAL(KIND=JPRB) :: ZZSUM, ZZRATIO
    REAL(KIND=JPRB) :: ZEPSILON

    REAL(KIND=JPRB) :: ZCOND1, ZQP

    REAL(KIND=JPRB) :: PSUM_SOLQA
    !*
    !     ------------------------------------------------------------------

    !     This COMDECK includes the Thermodynamical functions for the cy39
    !       ECMWF Physics package.
    !       Consistent with YOMCST Basic physics constants, assuming the
    !       partial pressure of water vapour is given by a first order
    !       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
    !       in YOETHF
    !       Two sets of functions are available. In the first set only the
    !       cases water or ice are distinguished by temperature.  This set
    !       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
    !       The second set considers, besides the two cases water and ice
    !       also a mix of both for the temperature range  YDTHF% RTICE < T <  YDTHF% RTWAT.
    !       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
    !       FKOOP modifies the ice saturation mixing ratio for homogeneous
    !       nucleation. FOE_DEWM_DT provides an approximate first derivative
    !       of FOEEWM.

    !       Depending on the consideration of mixed phases either the first
    !       set (e.g. surface, post-processing) or the second set
    !       (e.g. clouds, condensation, convection) should be used.

    !     ------------------------------------------------------------------
    !     *****************************************************************

    !                NO CONSIDERATION OF MIXED PHASES

    !     *****************************************************************
    REAL(KIND=JPRB) :: FOEDELTA
    REAL(KIND=JPRB) :: PTARE
    FOEDELTA(PTARE) = MAX(0.0_JPRB, SIGN(1.0_JPRB, PTARE - YDCST%RTT))

    !                  FOEDELTA = 1    water
    !                  FOEDELTA = 0    ice

    !     THERMODYNAMICAL FUNCTIONS .

    !     Pressure of water vapour at saturation
    !        INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOEEW, FOEDE, FOEDESU, FOELH, FOELDCP
    FOEEW(PTARE) = YDTHF%R2ES*EXP((YDTHF%R3LES*FOEDELTA(PTARE) + YDTHF%R3IES*(1.0_JPRB - FOEDELTA(PTARE)))*(PTARE - YDCST%RTT) /  &
    & (PTARE - (YDTHF%R4LES*FOEDELTA(PTARE) + YDTHF%R4IES*(1.0_JPRB - FOEDELTA(PTARE)))))

    FOEDE(PTARE) = (FOEDELTA(PTARE)*YDTHF%R5ALVCP + (1.0_JPRB - FOEDELTA(PTARE))*YDTHF%R5ALSCP) / (PTARE -  &
    & (YDTHF%R4LES*FOEDELTA(PTARE) + YDTHF%R4IES*(1.0_JPRB - FOEDELTA(PTARE))))**2

    FOEDESU(PTARE) = (FOEDELTA(PTARE)*YDTHF%R5LES + (1.0_JPRB - FOEDELTA(PTARE))*YDTHF%R5IES) / (PTARE -  &
    & (YDTHF%R4LES*FOEDELTA(PTARE) + YDTHF%R4IES*(1.0_JPRB - FOEDELTA(PTARE))))**2

    FOELH(PTARE) = FOEDELTA(PTARE)*YDCST%RLVTT + (1.0_JPRB - FOEDELTA(PTARE))*YDCST%RLSTT

    FOELDCP(PTARE) = FOEDELTA(PTARE)*YDTHF%RALVDCP + (1.0_JPRB - FOEDELTA(PTARE))*YDTHF%RALSDCP

    !     *****************************************************************

    !           CONSIDERATION OF MIXED PHASES

    !     *****************************************************************

    !     FOEALFA is calculated to distinguish the three cases:

    !                       FOEALFA=1            water phase
    !                       FOEALFA=0            ice phase
    !                       0 < FOEALFA < 1      mixed phase

    !               INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOEALFA
    FOEALFA(PTARE) = MIN(1.0_JPRB, ((MAX(YDTHF%RTICE, MIN(YDTHF%RTWAT, PTARE)) - YDTHF%RTICE)*YDTHF%RTWAT_RTICE_R)**2)


    !     Pressure of water vapour at saturation
    !        INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOEEWM, FOEDEM, FOELDCPM, FOELHM, FOE_DEWM_DT
    FOEEWM(PTARE) = YDTHF%R2ES*(FOEALFA(PTARE)*EXP(YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES)) + (1.0_JPRB -  &
    & FOEALFA(PTARE))*EXP(YDTHF%R3IES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4IES)))

    FOE_DEWM_DT(PTARE) = YDTHF%R2ES*(YDTHF%R3LES*FOEALFA(PTARE)*EXP(YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES)) &
    & *(YDCST%RTT - YDTHF%R4LES) / (PTARE - YDTHF%R4LES)**2 + YDTHF%R3IES*(1.0 - FOEALFA(PTARE))*EXP(YDTHF%R3IES*(PTARE -  &
    & YDCST%RTT) / (PTARE - YDTHF%R4IES))*(YDCST%RTT - YDTHF%R4IES) / (PTARE - YDTHF%R4IES)**2)

    FOEDEM(PTARE) = FOEALFA(PTARE)*YDTHF%R5ALVCP*(1.0_JPRB / (PTARE - YDTHF%R4LES)**2) + (1.0_JPRB - FOEALFA(PTARE)) &
    & *YDTHF%R5ALSCP*(1.0_JPRB / (PTARE - YDTHF%R4IES)**2)

    FOELDCPM(PTARE) = FOEALFA(PTARE)*YDTHF%RALVDCP + (1.0_JPRB - FOEALFA(PTARE))*YDTHF%RALSDCP

    FOELHM(PTARE) = FOEALFA(PTARE)*YDCST%RLVTT + (1.0_JPRB - FOEALFA(PTARE))*YDCST%RLSTT


    !     Temperature normalization for humidity background change of variable
    !        INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOETB
    FOETB(PTARE) = FOEALFA(PTARE)*YDTHF%R3LES*(YDCST%RTT - YDTHF%R4LES)*(1.0_JPRB / (PTARE - YDTHF%R4LES)**2) + (1.0_JPRB -  &
    & FOEALFA(PTARE))*YDTHF%R3IES*(YDCST%RTT - YDTHF%R4IES)*(1.0_JPRB / (PTARE - YDTHF%R4IES)**2)

    !     ------------------------------------------------------------------
    !     *****************************************************************

    !           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

    !     *****************************************************************

    !     FOEALFCU is calculated to distinguish the three cases:

    !                       FOEALFCU=1            water phase
    !                       FOEALFCU=0            ice phase
    !                       0 < FOEALFCU < 1      mixed phase

    !               INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOEALFCU
    FOEALFCU(PTARE) = MIN(1.0_JPRB, ((MAX(YDTHF%RTICECU, MIN(YDTHF%RTWAT, PTARE)) - YDTHF%RTICECU)*YDTHF%RTWAT_RTICECU_R)**2)


    !     Pressure of water vapour at saturation
    !        INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOEEWMCU, FOEDEMCU, FOELDCPMCU, FOELHMCU
    FOEEWMCU(PTARE) = YDTHF%R2ES*(FOEALFCU(PTARE)*EXP(YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES)) + (1.0_JPRB -  &
    & FOEALFCU(PTARE))*EXP(YDTHF%R3IES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4IES)))

    FOEDEMCU(PTARE) = FOEALFCU(PTARE)*YDTHF%R5ALVCP*(1.0_JPRB / (PTARE - YDTHF%R4LES)**2) + (1.0_JPRB - FOEALFCU(PTARE)) &
    & *YDTHF%R5ALSCP*(1.0_JPRB / (PTARE - YDTHF%R4IES)**2)

    FOELDCPMCU(PTARE) = FOEALFCU(PTARE)*YDTHF%RALVDCP + (1.0_JPRB - FOEALFCU(PTARE))*YDTHF%RALSDCP

    FOELHMCU(PTARE) = FOEALFCU(PTARE)*YDCST%RLVTT + (1.0_JPRB - FOEALFCU(PTARE))*YDCST%RLSTT
    !     ------------------------------------------------------------------

    !     Pressure of water vapour at saturation
    !     This one is for the WMO definition of saturation, i.e. always
    !     with respect to water.
    !
    !     Duplicate to FOEELIQ and FOEEICE for separate ice variable
    !     FOEELIQ always respect to water
    !     FOEEICE always respect to ice
    !     (could use FOEEW and FOEEWMO, but naming convention unclear)
    !     FOELSON returns e wrt liquid water using D Sonntag (1994, Met. Zeit.)
    !      - now recommended for use with radiosonde data (WMO CIMO guide, 2014)
    !      unlike the FOEE functions does not include 1/( YDCST% RETV+1.0_JPRB) factor

    REAL(KIND=JPRB) :: FOEEWMO, FOEELIQ, FOEEICE, FOELSON
    FOEEWMO(PTARE) = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES))
    FOEELIQ(PTARE) = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES))
    FOEEICE(PTARE) = YDTHF%R2ES*EXP(YDTHF%R3IES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4IES))
    FOELSON(PTARE) = EXP(-6096.9385_JPRB / PTARE + 21.2409642_JPRB - 2.711193E-2_JPRB*PTARE + 1.673952E-5_JPRB*PTARE**2 +  &
    & 2.433502_JPRB*LOG(PTARE))

    REAL(KIND=JPRB) :: FOEEWM_V, FOEEWMCU_V, FOELES_V, FOEIES_V
    REAL(KIND=JPRB) :: EXP1, EXP2
    FOELES_V(PTARE) = YDTHF%R3LES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4LES)
    FOEIES_V(PTARE) = YDTHF%R3IES*(PTARE - YDCST%RTT) / (PTARE - YDTHF%R4IES)
    FOEEWM_V(PTARE, EXP1, EXP2) = YDTHF%R2ES*(FOEALFA(PTARE)*EXP1 + (1.0_JPRB - FOEALFA(PTARE))*EXP2)
    FOEEWMCU_V(PTARE, EXP1, EXP2) = YDTHF%R2ES*(FOEALFCU(PTARE)*EXP1 + (1.0_JPRB - FOEALFCU(PTARE))*EXP2)
    ! (C) Copyright 1988- ECMWF.
    !
    ! This software is licensed under the terms of the Apache Licence Version 2.0
    ! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
    !
    ! In applying this licence, ECMWF does not waive the privileges and immunities
    ! granted to it by virtue of its status as an intergovernmental organisation
    ! nor does it submit to any jurisdiction.

    !*
    !     ------------------------------------------------------------------
    !     This COMDECK defines functions to be used in the cloud scheme
    !       other than the standard saturation vapour pressure
    !
    !       FKOOP modifies the ice saturation mixing ratio for homogeneous
    !       nucleation
    !
    !     note: PTARE is temperature and is definited in frttre.h
    !           which MUST be included before this function block
    !
    !     **********************************************
    !     KOOP formula for homogeneous nucleation of ice
    !     **********************************************
    !
    !               INPUT : PTARE = TEMPERATURE
    REAL(KIND=JPRB) :: FOKOOP
    FOKOOP(PTARE) = MIN(YDTHF%RKOOP1 - YDTHF%RKOOP2*PTARE, FOEELIQ(PTARE) / FOEEICE(PTARE))
!$acc routine seq


    !===============================================================================
    !IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)

    !===============================================================================
    !  0.0     Beginning of timestep book-keeping
    !----------------------------------------------------------------------


    !######################################################################
    !             0.  *** SET UP CONSTANTS ***
    !######################################################################

    ZEPSILON = 100._JPRB*EPSILON(ZEPSILON)

    ! ---------------------------------------------------------------------
    ! Set version of warm-rain autoconversion/accretion
    ! IWARMRAIN = 1 ! Sundquist
    ! IWARMRAIN = 2 ! Khairoutdinov and Kogan (2000)
    ! ---------------------------------------------------------------------
    IWARMRAIN = 2
    ! ---------------------------------------------------------------------
    ! Set version of rain evaporation
    ! IEVAPRAIN = 1 ! Sundquist
    ! IEVAPRAIN = 2 ! Abel and Boutle (2013)
    ! ---------------------------------------------------------------------
    IEVAPRAIN = 2
    ! ---------------------------------------------------------------------
    ! Set version of snow evaporation
    ! IEVAPSNOW = 1 ! Sundquist
    ! IEVAPSNOW = 2 ! New
    ! ---------------------------------------------------------------------
    IEVAPSNOW = 1
    ! ---------------------------------------------------------------------
    ! Set version of ice deposition
    ! IDEPICE = 1 ! Rotstayn (2001)
    ! IDEPICE = 2 ! New
    ! ---------------------------------------------------------------------
    IDEPICE = 1

    ! ---------------------
    ! Some simple constants
    ! ---------------------
    ZQTMST = 1.0_JPRB / PTSPHY
    ZGDCP = YDCST%RG / YDCST%RCPD
    ZRDCP = YDCST%RD / YDCST%RCPD
    ZCONS1A = YDCST%RCPD / (YDCST%RLMLT*YDCST%RG*YDECLDP%RTAUMEL)
    ZEPSEC = 1.E-14_JPRB
    ZRG_R = 1.0_JPRB / YDCST%RG
    ZRLDCP = 1.0_JPRB / (YDTHF%RALSDCP - YDTHF%RALVDCP)

    ! Note: Defined in module/yoecldp.F90
    ! NCLDQL=1    ! liquid cloud water
    ! NCLDQI=2    ! ice cloud water
    ! NCLDQR=3    ! rain water
    ! NCLDQS=4    ! snow
    ! NCLDQV=5    ! vapour

    ! -----------------------------------------------
    ! Define species phase, 0=vapour, 1=liquid, 2=ice
    ! -----------------------------------------------
    IPHASE(NCLDQV) = 0
    IPHASE(NCLDQL) = 1
    IPHASE(NCLDQR) = 1
    IPHASE(NCLDQI) = 2
    IPHASE(NCLDQS) = 2

    ! ---------------------------------------------------
    ! Set up melting/freezing index,
    ! if an ice category melts/freezes, where does it go?
    ! ---------------------------------------------------
    IMELT(NCLDQV) = -99
    IMELT(NCLDQL) = NCLDQI
    IMELT(NCLDQR) = NCLDQS
    IMELT(NCLDQI) = NCLDQR
    IMELT(NCLDQS) = NCLDQR

    ! -----------------------------------------------
    ! INITIALIZATION OF OUTPUT TENDENCIES
    ! -----------------------------------------------
!$acc loop seq
    DO JK=1,KLEV
      TENDENCY_LOC_T(JL, JK) = 0.0_JPRB
      TENDENCY_LOC_Q(JL, JK) = 0.0_JPRB
      TENDENCY_LOC_A(JL, JK) = 0.0_JPRB
    END DO
!$acc loop seq
    DO JM=1,NCLV - 1
      DO JK=1,KLEV
        TENDENCY_LOC_CLD(JL, JK, JM) = 0.0_JPRB
      END DO
    END DO

    !-- These were uninitialized : meaningful only when we compare error differences
!$acc loop seq
    DO JK=1,KLEV
      PCOVPTOT(JL, JK) = 0.0_JPRB
      TENDENCY_LOC_CLD(JL, JK, NCLV) = 0.0_JPRB
    END DO

    ! -------------------------
    ! set up fall speeds in m/s
    ! -------------------------
    ZVQX(NCLDQV) = 0.0_JPRB
    ZVQX(NCLDQL) = 0.0_JPRB
    ZVQX(NCLDQI) = YDECLDP%RVICE
    ZVQX(NCLDQR) = YDECLDP%RVRAIN
    ZVQX(NCLDQS) = YDECLDP%RVSNOW
    LLFALL(:) = .false.
!$acc loop seq
    DO JM=1,NCLV
      IF (ZVQX(JM) > 0.0_JPRB)       LLFALL(JM) = .true.
      ! falling species
    END DO
    ! Set LLFALL to false for ice (but ice still sediments!)
    ! Need to rationalise this at some point
    LLFALL(NCLDQI) = .false.

    !######################################################################
    !             1.  *** INITIAL VALUES FOR VARIABLES ***
    !######################################################################


    ! ----------------------
    ! non CLV initialization
    ! ----------------------
!$acc loop seq
    DO JK=1,KLEV
      ZTP1(JL, JK) = PT(JL, JK) + PTSPHY*TENDENCY_TMP_T(JL, JK)
      ZQX(JL, JK, NCLDQV) = PQ(JL, JK) + PTSPHY*TENDENCY_TMP_Q(JL, JK)
      ZQX0(JL, JK, NCLDQV) = PQ(JL, JK) + PTSPHY*TENDENCY_TMP_Q(JL, JK)
      ZA(JL, JK) = PA(JL, JK) + PTSPHY*TENDENCY_TMP_A(JL, JK)
      ZAORIG(JL, JK) = PA(JL, JK) + PTSPHY*TENDENCY_TMP_A(JL, JK)
    END DO

    ! -------------------------------------
    ! initialization for CLV family
    ! -------------------------------------
!$acc loop seq
    DO JM=1,NCLV - 1
      DO JK=1,KLEV
        ZQX(JL, JK, JM) = PCLV(JL, JK, JM) + PTSPHY*TENDENCY_TMP_CLD(JL, JK, JM)
        ZQX0(JL, JK, JM) = PCLV(JL, JK, JM) + PTSPHY*TENDENCY_TMP_CLD(JL, JK, JM)
      END DO
    END DO

    !-------------
    ! zero arrays
    !-------------
!$acc loop seq
    DO JM=1,NCLV
      DO JK=1,KLEV + 1
        ZPFPLSX(JL, JK, JM) = 0.0_JPRB          ! precip fluxes
      END DO
    END DO

!$acc loop seq
    DO JM=1,NCLV
      DO JK=1,KLEV
        ZQXN2D(JL, JK, JM) = 0.0_JPRB          ! end of timestep values in 2D
        ZLNEG(JL, JK, JM) = 0.0_JPRB          ! negative input check
      END DO
    END DO

    PRAINFRAC_TOPRFZ(JL) = 0.0_JPRB      ! rain fraction at top of refreezing layer
    LLRAINLIQ = .true.      ! Assume all raindrops are liquid initially

    ! ----------------------------------------------------
    ! Tidy up very small cloud cover or total cloud water
    ! ----------------------------------------------------
!$acc loop seq
    DO JK=1,KLEV
      IF (ZQX(JL, JK, NCLDQL) + ZQX(JL, JK, NCLDQI) < YDECLDP%RLMIN .or. ZA(JL, JK) < YDECLDP%RAMIN) THEN

        ! Evaporate small cloud liquid water amounts
        ZLNEG(JL, JK, NCLDQL) = ZLNEG(JL, JK, NCLDQL) + ZQX(JL, JK, NCLDQL)
        ZQADJ = ZQX(JL, JK, NCLDQL)*ZQTMST
        TENDENCY_LOC_Q(JL, JK) = TENDENCY_LOC_Q(JL, JK) + ZQADJ
        TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) - YDTHF%RALVDCP*ZQADJ
        ZQX(JL, JK, NCLDQV) = ZQX(JL, JK, NCLDQV) + ZQX(JL, JK, NCLDQL)
        ZQX(JL, JK, NCLDQL) = 0.0_JPRB

        ! Evaporate small cloud ice water amounts
        ZLNEG(JL, JK, NCLDQI) = ZLNEG(JL, JK, NCLDQI) + ZQX(JL, JK, NCLDQI)
        ZQADJ = ZQX(JL, JK, NCLDQI)*ZQTMST
        TENDENCY_LOC_Q(JL, JK) = TENDENCY_LOC_Q(JL, JK) + ZQADJ
        TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) - YDTHF%RALSDCP*ZQADJ
        ZQX(JL, JK, NCLDQV) = ZQX(JL, JK, NCLDQV) + ZQX(JL, JK, NCLDQI)
        ZQX(JL, JK, NCLDQI) = 0.0_JPRB

        ! Set cloud cover to zero
        ZA(JL, JK) = 0.0_JPRB

      END IF
    END DO

    ! ---------------------------------
    ! Tidy up small CLV variables
    ! ---------------------------------
    !DIR$ IVDEP
!$acc loop seq
    DO JM=1,NCLV - 1
      !DIR$ IVDEP
      DO JK=1,KLEV
        !DIR$ IVDEP
        IF (ZQX(JL, JK, JM) < YDECLDP%RLMIN) THEN
          ZLNEG(JL, JK, JM) = ZLNEG(JL, JK, JM) + ZQX(JL, JK, JM)
          ZQADJ = ZQX(JL, JK, JM)*ZQTMST
          TENDENCY_LOC_Q(JL, JK) = TENDENCY_LOC_Q(JL, JK) + ZQADJ
          IF (IPHASE(JM) == 1)           TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) - YDTHF%RALVDCP*ZQADJ
          IF (IPHASE(JM) == 2)           TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) - YDTHF%RALSDCP*ZQADJ
          ZQX(JL, JK, NCLDQV) = ZQX(JL, JK, NCLDQV) + ZQX(JL, JK, JM)
          ZQX(JL, JK, JM) = 0.0_JPRB
        END IF
      END DO
    END DO


    ! ------------------------------
    ! Define saturation values
    ! ------------------------------
!$acc loop seq
    DO JK=1,KLEV
      !----------------------------------------
      ! old *diagnostic* mixed phase saturation
      !----------------------------------------
      ZFOEALFA(JL, JK) = FOEALFA(ZTP1(JL, JK))
      ZFOEEWMT(JL, JK) = MIN(FOEEWM(ZTP1(JL, JK)) / PAP(JL, JK), 0.5_JPRB)
      ZQSMIX(JL, JK) = ZFOEEWMT(JL, JK)
      ZQSMIX(JL, JK) = ZQSMIX(JL, JK) / (1.0_JPRB - YDCST%RETV*ZQSMIX(JL, JK))

      !---------------------------------------------
      ! ice saturation T<273K
      ! liquid water saturation for T>273K
      !---------------------------------------------
      ZALFA = FOEDELTA(ZTP1(JL, JK))
      ZFOEEW(JL, JK) = MIN((ZALFA*FOEELIQ(ZTP1(JL, JK)) + (1.0_JPRB - ZALFA)*FOEEICE(ZTP1(JL, JK))) / PAP(JL, JK), 0.5_JPRB)
      ZFOEEW(JL, JK) = MIN(0.5_JPRB, ZFOEEW(JL, JK))
      ZQSICE(JL, JK) = ZFOEEW(JL, JK) / (1.0_JPRB - YDCST%RETV*ZFOEEW(JL, JK))

      !----------------------------------
      ! liquid water saturation
      !----------------------------------
      ZFOEELIQT(JL, JK) = MIN(FOEELIQ(ZTP1(JL, JK)) / PAP(JL, JK), 0.5_JPRB)
      ZQSLIQ(JL, JK) = ZFOEELIQT(JL, JK)
      ZQSLIQ(JL, JK) = ZQSLIQ(JL, JK) / (1.0_JPRB - YDCST%RETV*ZQSLIQ(JL, JK))

      !   !----------------------------------
      !   ! ice water saturation
      !   !----------------------------------
      !   ZFOEEICET(JL,JK)=MIN(FOEEICE(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      !   ZQSICE(JL,JK)=ZFOEEICET(JL,JK)
      !   ZQSICE(JL,JK)=ZQSICE(JL,JK)/(1.0_JPRB-YDCST%RETV*ZQSICE(JL,JK))

    END DO

!$acc loop seq
    DO JK=1,KLEV


      !------------------------------------------
      ! Ensure cloud fraction is between 0 and 1
      !------------------------------------------
      ZA(JL, JK) = MAX(0.0_JPRB, MIN(1.0_JPRB, ZA(JL, JK)))

      !-------------------------------------------------------------------
      ! Calculate liq/ice fractions (no longer a diagnostic relationship)
      !-------------------------------------------------------------------
      ZLI(JL, JK) = ZQX(JL, JK, NCLDQL) + ZQX(JL, JK, NCLDQI)
      IF (ZLI(JL, JK) > YDECLDP%RLMIN) THEN
        ZLIQFRAC(JL, JK) = ZQX(JL, JK, NCLDQL) / ZLI(JL, JK)
        ZICEFRAC(JL, JK) = 1.0_JPRB - ZLIQFRAC(JL, JK)
      ELSE
        ZLIQFRAC(JL, JK) = 0.0_JPRB
        ZICEFRAC(JL, JK) = 0.0_JPRB
      END IF

    END DO
    !######################################################################
    !        2.       *** CONSTANTS AND PARAMETERS ***
    !######################################################################
    !  Calculate L in updrafts of bl-clouds
    !  Specify QS, P/PS for tropopause (for c2)
    !  And initialize variables
    !------------------------------------------

    !---------------------------------
    ! Find tropopause level (ZTRPAUS)
    !---------------------------------
    ZTRPAUS = 0.1_JPRB
    ZPAPHD = 1.0_JPRB / PAPH(JL, KLEV + 1)

!$acc loop seq
    DO JK=1,KLEV - 1
      ZSIG = PAP(JL, JK)*ZPAPHD
      IF (ZSIG > 0.1_JPRB .and. ZSIG < 0.4_JPRB .and. ZTP1(JL, JK) > ZTP1(JL, JK + 1)) THEN
        ZTRPAUS = ZSIG
      END IF
    END DO
    !-----------------------------
    ! Reset single level variables
    !-----------------------------

    ZANEWM1 = 0.0_JPRB
    ZDA(:) = 0.0_JPRB
    ZCOVPCLR = 0.0_JPRB
    ZCOVPMAX = 0.0_JPRB
    ZCOVPTOT = 0.0_JPRB
    ZCLDTOPDIST = 0.0_JPRB

    !######################################################################
    !           3.       *** PHYSICS ***
    !######################################################################


    !----------------------------------------------------------------------
    !                       START OF VERTICAL LOOP
    !----------------------------------------------------------------------

!$acc loop seq
    DO JK=YDECLDP%NCLDTOP,KLEV

      !----------------------------------------------------------------------
      ! 3.0 INITIALIZE VARIABLES
      !----------------------------------------------------------------------

      !---------------------------------
      ! First guess microphysics
      !---------------------------------
      DO JM=1,NCLV
        ZQXFG(JM, JK) = ZQX(JL, JK, JM)
      END DO

      !---------------------------------
      ! Set KLON arrays to zero
      !---------------------------------
      ZLICLD(JK) = 0.0_JPRB
      ZRAINAUT = 0.0_JPRB        ! currently needed for diags
      ZRAINACC = 0.0_JPRB        ! currently needed for diags
      ZSNOWAUT = 0.0_JPRB        ! needed
      ZLDEFR(JK) = 0.0_JPRB
      ZACUST = 0.0_JPRB        ! set later when needed
      ZQPRETOT(JK) = 0.0_JPRB
      ZLFINALSUM(JK) = 0.0_JPRB

      ! Required for first guess call
      ZLCOND1 = 0.0_JPRB
      ZLCOND2 = 0.0_JPRB
      ZSUPSAT = 0.0_JPRB
      ZLEVAPL = 0.0_JPRB
      ZLEVAPI = 0.0_JPRB

      !-------------------------------------
      ! solvers for cloud fraction
      !-------------------------------------
      ZSOLAB(JK) = 0.0_JPRB
      ZSOLAC(JK) = 0.0_JPRB

      ZICETOT = 0.0_JPRB

      !------------------------------------------
      ! reset matrix so missing pathways are set
      !------------------------------------------
      DO JM=1,NCLV
        DO JN=1,NCLV
          ZSOLQB(JN, JM, JK) = 0.0_JPRB
          ZSOLQA(JN, JM, JK) = 0.0_JPRB
        END DO
      END DO

      !----------------------------------
      ! reset new microphysics variables
      !----------------------------------
      DO JM=1,NCLV
        ZFALLSRCE(JM, JK) = 0.0_JPRB
        ZFALLSINK(JM, JK) = 0.0_JPRB
        ZCONVSRCE(JM, JK) = 0.0_JPRB
        ZCONVSINK(JM, JK) = 0.0_JPRB
        ZPSUPSATSRCE(JM, JK) = 0.0_JPRB
        ZRATIO(JM) = 0.0_JPRB
      END DO


      !-------------------------
      ! derived variables needed
      !-------------------------

      ZDP(JK) = PAPH(JL, JK + 1) - PAPH(JL, JK)        ! dp
      ZGDP = YDCST%RG / ZDP(JK)        ! g/dp
      ZRHO(JK) = PAP(JL, JK) / (YDCST%RD*ZTP1(JL, JK))        ! p/RT air density

      ZDTGDP(JK) = PTSPHY*ZGDP        ! dt g/dp
      ZRDTGDP(JK) = ZDP(JK)*(1.0_JPRB / (PTSPHY*YDCST%RG))        ! 1/(dt g/dp)

      IF (JK > 1)       ZDTGDPF = (PTSPHY*YDCST%RG) / (PAP(JL, JK) - PAP(JL, JK - 1))

      !------------------------------------
      ! Calculate dqs/dT correction factor
      !------------------------------------
      ! Reminder: YDCST%RETV=YDCST%RV/YDCST%RD-1

      ! liquid
      ZFACW = YDTHF%R5LES / ((ZTP1(JL, JK) - YDTHF%R4LES)**2)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZFOEELIQT(JL, JK))
      ZDQSLIQDT = ZFACW*ZCOR*ZQSLIQ(JL, JK)
      ZCORQSLIQ(JK) = 1.0_JPRB + YDTHF%RALVDCP*ZDQSLIQDT

      ! ice
      ZFACI = YDTHF%R5IES / ((ZTP1(JL, JK) - YDTHF%R4IES)**2)
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZFOEEW(JL, JK))
      ZDQSICEDT = ZFACI*ZCOR*ZQSICE(JL, JK)
      ZCORQSICE(JK) = 1.0_JPRB + YDTHF%RALSDCP*ZDQSICEDT

      ! diagnostic mixed
      ZALFAW = ZFOEALFA(JL, JK)
      ZALFAWM = ZALFAW
      ZFAC = ZALFAW*ZFACW + (1.0_JPRB - ZALFAW)*ZFACI
      ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZFOEEWMT(JL, JK))
      ZDQSMIXDT(JK) = ZFAC*ZCOR*ZQSMIX(JL, JK)
      ZCORQSMIX = 1.0_JPRB + FOELDCPM(ZTP1(JL, JK))*ZDQSMIXDT(JK)

      ! evaporation/sublimation limits
      ZEVAPLIMMIX(JK) = MAX((ZQSMIX(JL, JK) - ZQX(JL, JK, NCLDQV)) / ZCORQSMIX, 0.0_JPRB)
      ZEVAPLIMLIQ = MAX((ZQSLIQ(JL, JK) - ZQX(JL, JK, NCLDQV)) / ZCORQSLIQ(JK), 0.0_JPRB)
      ZEVAPLIMICE(JK) = MAX((ZQSICE(JL, JK) - ZQX(JL, JK, NCLDQV)) / ZCORQSICE(JK), 0.0_JPRB)

      !--------------------------------
      ! in-cloud consensate amount
      !--------------------------------
      ZTMPA = 1.0_JPRB / MAX(ZA(JL, JK), ZEPSEC)
      ZLIQCLD = ZQX(JL, JK, NCLDQL)*ZTMPA
      ZICECLD(JK) = ZQX(JL, JK, NCLDQI)*ZTMPA
      ZLICLD(JK) = ZLIQCLD + ZICECLD(JK)


      !------------------------------------------------
      ! Evaporate very small amounts of liquid and ice
      !------------------------------------------------

      IF (ZQX(JL, JK, NCLDQL) < YDECLDP%RLMIN) THEN
        ZSOLQA(NCLDQV, NCLDQL, JK) = ZQX(JL, JK, NCLDQL)
        ZSOLQA(NCLDQL, NCLDQV, JK) = -ZQX(JL, JK, NCLDQL)
      END IF

      IF (ZQX(JL, JK, NCLDQI) < YDECLDP%RLMIN) THEN
        ZSOLQA(NCLDQV, NCLDQI, JK) = ZQX(JL, JK, NCLDQI)
        ZSOLQA(NCLDQI, NCLDQV, JK) = -ZQX(JL, JK, NCLDQI)
      END IF


      !---------------------------------------------------------------------
      !  3.1  ICE SUPERSATURATION ADJUSTMENT
      !---------------------------------------------------------------------
      ! Note that the supersaturation adjustment is made with respect to
      ! liquid saturation:  when T>0C
      ! ice saturation:     when T<0C
      !                     with an adjustment made to allow for ice
      !                     supersaturation in the clear sky
      ! Note also that the KOOP factor automatically clips the supersaturation
      ! to a maximum set by the liquid water saturation mixing ratio
      ! important for temperatures near to but below 0C
      !-----------------------------------------------------------------------

      !DIR$ NOFUSION

      !-----------------------------------
      ! 3.1.1 Supersaturation limit (from Koop)
      !-----------------------------------
      ! Needs to be set for all temperatures
      ZFOKOOP(JK) = FOKOOP(ZTP1(JL, JK))

      IF (ZTP1(JL, JK) >= YDCST%RTT .or. YDECLDP%NSSOPT == 0) THEN
        ZFAC = 1.0_JPRB
        ZFACI = 1.0_JPRB
      ELSE
        ZFAC = ZA(JL, JK) + ZFOKOOP(JK)*(1.0_JPRB - ZA(JL, JK))
        ZFACI = PTSPHY / YDECLDP%RKOOPTAU
      END IF

      !-------------------------------------------------------------------
      ! 3.1.2 Calculate supersaturation wrt Koop including dqs/dT
      !       correction factor
      ! [#Note: QSICE or QSLIQ]
      !-------------------------------------------------------------------

      ! Calculate supersaturation to add to cloud
      IF (ZA(JL, JK) > 1.0_JPRB - YDECLDP%RAMIN) THEN
        ZSUPSAT = MAX((ZQX(JL, JK, NCLDQV) - ZFAC*ZQSICE(JL, JK)) / ZCORQSICE(JK), 0.0_JPRB)
      ELSE
        ! Calculate environmental humidity supersaturation
        ZQP1ENV = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSICE(JL, JK)) / MAX(1.0_JPRB - ZA(JL, JK), ZEPSILON)
        !& SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
        ZSUPSAT = MAX(((1.0_JPRB - ZA(JL, JK))*(ZQP1ENV - ZFAC*ZQSICE(JL, JK))) / ZCORQSICE(JK), 0.0_JPRB)
      END IF

      !-------------------------------------------------------------------
      ! Here the supersaturation is turned into liquid water
      ! However, if the temperature is below the threshold for homogeneous
      ! freezing then the supersaturation is turned instantly to ice.
      !--------------------------------------------------------------------

      IF (ZSUPSAT > ZEPSEC) THEN

        IF (ZTP1(JL, JK) > YDECLDP%RTHOMO) THEN
          ! Turn supersaturation into liquid water
          ZSOLQA(NCLDQL, NCLDQV, JK) = ZSOLQA(NCLDQL, NCLDQV, JK) + ZSUPSAT
          ZSOLQA(NCLDQV, NCLDQL, JK) = ZSOLQA(NCLDQV, NCLDQL, JK) - ZSUPSAT
          ! Include liquid in first guess
          ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) + ZSUPSAT
        ELSE
          ! Turn supersaturation into ice water
          ZSOLQA(NCLDQI, NCLDQV, JK) = ZSOLQA(NCLDQI, NCLDQV, JK) + ZSUPSAT
          ZSOLQA(NCLDQV, NCLDQI, JK) = ZSOLQA(NCLDQV, NCLDQI, JK) - ZSUPSAT
          ! Add ice to first guess for deposition term
          ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + ZSUPSAT
        END IF

        ! Increase cloud amount using RKOOPTAU timescale
        ZSOLAC(JK) = (1.0_JPRB - ZA(JL, JK))*ZFACI

      END IF

      !-------------------------------------------------------
      ! 3.1.3 Include supersaturation from previous timestep
      ! (Calculated in sltENDIF semi-lagrangian LDSLPHY=T)
      !-------------------------------------------------------
      IF (PSUPSAT(JL, JK) > ZEPSEC) THEN
        IF (ZTP1(JL, JK) > YDECLDP%RTHOMO) THEN
          ! Turn supersaturation into liquid water
          ZSOLQA(NCLDQL, NCLDQL, JK) = ZSOLQA(NCLDQL, NCLDQL, JK) + PSUPSAT(JL, JK)
          ZPSUPSATSRCE(NCLDQL, JK) = PSUPSAT(JL, JK)
          ! Add liquid to first guess for deposition term
          ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) + PSUPSAT(JL, JK)
          ! Store cloud budget diagnostics if required
        ELSE
          ! Turn supersaturation into ice water
          ZSOLQA(NCLDQI, NCLDQI, JK) = ZSOLQA(NCLDQI, NCLDQI, JK) + PSUPSAT(JL, JK)
          ZPSUPSATSRCE(NCLDQI, JK) = PSUPSAT(JL, JK)
          ! Add ice to first guess for deposition term
          ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + PSUPSAT(JL, JK)
          ! Store cloud budget diagnostics if required
        END IF

        ! Increase cloud amount using RKOOPTAU timescale
        ZSOLAC(JK) = (1.0_JPRB - ZA(JL, JK))*ZFACI
        ! Store cloud budget diagnostics if required
      END IF

      ! on JL

      !---------------------------------------------------------------------
      !  3.2  DETRAINMENT FROM CONVECTION
      !---------------------------------------------------------------------
      ! * Diagnostic T-ice/liq split retained for convection
      !    Note: This link is now flexible and a future convection
      !    scheme can detrain explicit seperate budgets of:
      !    cloud water, ice, rain and snow
      ! * There is no (1-ZA) multiplier term on the cloud detrainment
      !    term, since is now written in mass-flux terms
      ! [#Note: Should use ZFOEALFACU used in convection rather than ZFOEALFA]
      !---------------------------------------------------------------------
      IF (JK < KLEV .and. JK >= YDECLDP%NCLDTOP) THEN


        PLUDE(JL, JK) = PLUDE(JL, JK)*ZDTGDP(JK)

        IF (LDCUM(JL) .and. PLUDE(JL, JK) > YDECLDP%RLMIN .and. PLU(JL, JK + 1) > ZEPSEC) THEN

          ZSOLAC(JK) = ZSOLAC(JK) + PLUDE(JL, JK) / PLU(JL, JK + 1)
          ! *diagnostic temperature split*
          ZALFAW = ZFOEALFA(JL, JK)
          ZCONVSRCE(NCLDQL, JK) = ZALFAW*PLUDE(JL, JK)
          ZCONVSRCE(NCLDQI, JK) = (1.0_JPRB - ZALFAW)*PLUDE(JL, JK)
          ZSOLQA(NCLDQL, NCLDQL, JK) = ZSOLQA(NCLDQL, NCLDQL, JK) + ZCONVSRCE(NCLDQL, JK)
          ZSOLQA(NCLDQI, NCLDQI, JK) = ZSOLQA(NCLDQI, NCLDQI, JK) + ZCONVSRCE(NCLDQI, JK)

        ELSE

          PLUDE(JL, JK) = 0.0_JPRB

        END IF
        ! *convective snow detrainment source
        IF (LDCUM(JL))         ZSOLQA(NCLDQS, NCLDQS, JK) = ZSOLQA(NCLDQS, NCLDQS, JK) + PSNDE(JL, JK)*ZDTGDP(JK)


      END IF
      ! JK<KLEV

    END DO
    ! Loki - loop-fission
    DO JK=YDECLDP%NCLDTOP,KLEV

      CALL section3p4(KLON, YDTHF, ZICEFRAC, PAPH, ZFOKOOP, YDECLDP, PAP, ZQSICE, PTSPHY, ZQTMST, JK, ZLICLD, ZEVAPLIMMIX, ZLI, ZDP,  &
      & JL, YDCST, PMFU, PVERVEL, ZRDCP, ZLDEFR, PHRSW, PMFD, ZQX, KLEV, ZA, ZEPSEC, PHRLW, ZLIQFRAC, ZQSMIX, ZSOLQA, ZSOLAC,  &
      & ZQXFG, ZTP1, ZLCONDLIM, ZCOR, ZQSAT, ZLEVAPI, ZMFDN, ZQP, ZCDMAX, ZRHC, ZDTFORC, ZDTDIAB, ZLCOND2, ZACOND, ZZDL, ZDQS,  &
      & ZDTDP, ZLEVAPL, ZLEVAP, ZZZDT, ZWTOT, ZCOND1, ZQE, ZDPMXDT, ZLCOND1, ZTOLD, ZSIGK, LLFLAG, ZQOLD, ZFAC, ZCOND)

    END DO
    ! Loki - loop-fission
    DO JK=YDECLDP%NCLDTOP,KLEV

      CALL section3p7(KLON, KLEV, YDTHF, PTSPHY, JK, ZA, ZFOKOOP, ZRHO, YDECLDP, PAP, ZDP, YDCST, JL, IDEPICE, ZICECLD, ZTP1, ZQXFG, ZCLDTOPDIST,  &
      & ZSOLQA, ZCVDS, ZTERM1, ZCORRFAC, ZADD, ZICE0, ZFACX1I, ZVPICE, ZCORRFAC2, ZPR02, ZBDD, ZTCG, ZICENUCLEI, ZTERM2, ZVPLIQ,  &
      & ZAPLUSB, ZINEW, ZDEPOS, ZINFACTOR)

    END DO
    ! Loki - loop-fission
    DO JK=YDECLDP%NCLDTOP,KLEV

!$loki region-to-call name( section3 )
      !---------------------------------------------------------------------
      !  3.3  SUBSIDENCE COMPENSATING CONVECTIVE UPDRAUGHTS
      !---------------------------------------------------------------------
      ! Three terms:
      ! * Convective subsidence source of cloud from layer above
      ! * Evaporation of cloud within the layer
      ! * Subsidence sink of cloud to the layer below (Implicit solution)
      !---------------------------------------------------------------------

      !-----------------------------------------------
      ! Subsidence source from layer above
      !               and
      ! Evaporation of cloud within the layer
      !-----------------------------------------------
      IF (JK > YDECLDP%NCLDTOP) THEN

        ZMF = MAX(0.0_JPRB, (PMFU(JL, JK) + PMFD(JL, JK))*ZDTGDP(JK))
        ZACUST = ZMF*ZANEWM1

        DO JM=1,NCLV
          IF (.not.LLFALL(JM) .and. IPHASE(JM) > 0) THEN
            ZLCUST(JM) = ZMF*ZQXNM1(JM)
            ! record total flux for enthalpy budget:
            ZCONVSRCE(JM, JK) = ZCONVSRCE(JM, JK) + ZLCUST(JM)
          END IF
        END DO

        ! Now have to work out how much liquid evaporates at arrival point
        ! since there is no prognostic memory for in-cloud humidity, i.e.
        ! we always assume cloud is saturated.

        ZDTDP = (ZRDCP*0.5_JPRB*(ZTP1(JL, JK - 1) + ZTP1(JL, JK))) / PAPH(JL, JK)
        ZDTFORC = ZDTDP*(PAP(JL, JK) - PAP(JL, JK - 1))
        ![#Note: Diagnostic mixed phase should be replaced below]
        ZDQS = ZANEWM1*ZDTFORC*ZDQSMIXDT(JK)

        DO JM=1,NCLV
          IF (.not.LLFALL(JM) .and. IPHASE(JM) > 0) THEN
            ZLFINAL = MAX(0.0_JPRB, ZLCUST(JM) - ZDQS)              !lim to zero
            ! no supersaturation allowed incloud ---V
            ZEVAP = MIN((ZLCUST(JM) - ZLFINAL), ZEVAPLIMMIX(JK))
            !          ZEVAP=0.0_JPRB
            ZLFINAL = ZLCUST(JM) - ZEVAP
            ZLFINALSUM(JK) = ZLFINALSUM(JK) + ZLFINAL              ! sum

            ZSOLQA(JM, JM, JK) = ZSOLQA(JM, JM, JK) + ZLCUST(JM)              ! whole sum
            ZSOLQA(NCLDQV, JM, JK) = ZSOLQA(NCLDQV, JM, JK) + ZEVAP
            ZSOLQA(JM, NCLDQV, JK) = ZSOLQA(JM, NCLDQV, JK) - ZEVAP
          END IF
        END DO

        !  Reset the cloud contribution if no cloud water survives to this level:
        IF (ZLFINALSUM(JK) < ZEPSEC)         ZACUST = 0.0_JPRB
        ZSOLAC(JK) = ZSOLAC(JK) + ZACUST

      END IF
      ! on  JK>NCLDTOP

      !---------------------------------------------------------------------
      ! Subsidence sink of cloud to the layer below
      ! (Implicit - re. CFL limit on convective mass flux)
      !---------------------------------------------------------------------


      IF (JK < KLEV) THEN

        ZMFDN = MAX(0.0_JPRB, (PMFU(JL, JK + 1) + PMFD(JL, JK + 1))*ZDTGDP(JK))

        ZSOLAB(JK) = ZSOLAB(JK) + ZMFDN
        ZSOLQB(NCLDQL, NCLDQL, JK) = ZSOLQB(NCLDQL, NCLDQL, JK) + ZMFDN
        ZSOLQB(NCLDQI, NCLDQI, JK) = ZSOLQB(NCLDQI, NCLDQI, JK) + ZMFDN

        ! Record sink for cloud budget and enthalpy budget diagnostics
        ZCONVSINK(NCLDQL, JK) = ZMFDN
        ZCONVSINK(NCLDQI, JK) = ZMFDN

      END IF


      !----------------------------------------------------------------------
      ! 3.4  EROSION OF CLOUDS BY TURBULENT MIXING
      !----------------------------------------------------------------------
      ! NOTE: In default tiedtke scheme this process decreases the cloud
      !       area but leaves the specific cloud water content
      !       within clouds unchanged
      !----------------------------------------------------------------------

      ! ------------------------------
      ! Define turbulent erosion rate
      ! ------------------------------
      ZLDIFDT = YDECLDP%RCLDIFF*PTSPHY        !original version
      !Increase by factor of 5 for convective points
      IF (KTYPE(JL) > 0 .and. PLUDE(JL, JK) > ZEPSEC)       ZLDIFDT = YDECLDP%RCLDIFF_CONVI*ZLDIFDT

      ! At the moment, works on mixed RH profile and partitioned ice/liq fraction
      ! so that it is similar to previous scheme
      ! Should apply RHw for liquid cloud and RHi for ice cloud separately
      IF (ZLI(JL, JK) > ZEPSEC) THEN
        ! Calculate environmental humidity
        !      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
        !    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
        !      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
        ZE = ZLDIFDT*MAX(ZQSMIX(JL, JK) - ZQX(JL, JK, NCLDQV), 0.0_JPRB)
        ZLEROS = ZA(JL, JK)*ZE
        ZLEROS = MIN(ZLEROS, ZEVAPLIMMIX(JK))
        ZLEROS = MIN(ZLEROS, ZLI(JL, JK))
        ZAEROS = ZLEROS / ZLICLD(JK)          !if linear term

        ! Erosion is -ve LINEAR in L,A
        ZSOLAC(JK) = ZSOLAC(JK) - ZAEROS          !linear

        ZSOLQA(NCLDQV, NCLDQL, JK) = ZSOLQA(NCLDQV, NCLDQL, JK) + ZLIQFRAC(JL, JK)*ZLEROS
        ZSOLQA(NCLDQL, NCLDQV, JK) = ZSOLQA(NCLDQL, NCLDQV, JK) - ZLIQFRAC(JL, JK)*ZLEROS
        ZSOLQA(NCLDQV, NCLDQI, JK) = ZSOLQA(NCLDQV, NCLDQI, JK) + ZICEFRAC(JL, JK)*ZLEROS
        ZSOLQA(NCLDQI, NCLDQV, JK) = ZSOLQA(NCLDQI, NCLDQV, JK) - ZICEFRAC(JL, JK)*ZLEROS

      END IF
      ! Loki region-hoist group(section3p4) - region hoisted

      ! Foobar

      ! Loki region-hoist group(section3p7) - region hoisted

      !######################################################################
      !              4  *** PRECIPITATION PROCESSES ***
      !######################################################################

      !----------------------------------
      ! revise in-cloud consensate amount
      !----------------------------------
      ZTMPA = 1.0_JPRB / MAX(ZA(JL, JK), ZEPSEC)
      ZLIQCLD = ZQXFG(NCLDQL, JK)*ZTMPA
      ZICECLD(JK) = ZQXFG(NCLDQI, JK)*ZTMPA
      ZLICLD(JK) = ZLIQCLD + ZICECLD(JK)

      !----------------------------------------------------------------------
      ! 4.2 SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
      !     now that rain, snow, graupel species are prognostic
      !     the precipitation flux can be defined directly level by level
      !     There is no vertical memory required from the flux variable
      !----------------------------------------------------------------------

      DO JM=1,NCLV
        IF (LLFALL(JM) .or. JM == NCLDQI) THEN
          !------------------------
          ! source from layer above
          !------------------------
          IF (JK > YDECLDP%NCLDTOP) THEN
            ZFALLSRCE(JM, JK) = ZPFPLSX(JL, JK, JM)*ZDTGDP(JK)
            ZSOLQA(JM, JM, JK) = ZSOLQA(JM, JM, JK) + ZFALLSRCE(JM, JK)
            ZQXFG(JM, JK) = ZQXFG(JM, JK) + ZFALLSRCE(JM, JK)
            ! use first guess precip----------V
            ZQPRETOT(JK) = ZQPRETOT(JK) + ZQXFG(JM, JK)
          END IF
          !-------------------------------------------------
          ! sink to next layer, constant fall speed
          !-------------------------------------------------
          ! if aerosol effect then override
          !  note that for T>233K this is the same as above.
          IF (YDECLDP%LAERICESED .and. JM == NCLDQI) THEN
            ZRE_ICE = PRE_ICE(JL, JK)
            ! The exponent value is from
            ! Morrison et al. JAS 2005 Appendix
            ZVQX(NCLDQI) = 0.002_JPRB*ZRE_ICE**1.0_JPRB
          END IF
          ZFALL = ZVQX(JM)*ZRHO(JK)
          !-------------------------------------------------
          ! modified by Heymsfield and Iaquinta JAS 2000
          !-------------------------------------------------
          ! ZFALL = ZFALL*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
          !            &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))

          ZFALLSINK(JM, JK) = ZDTGDP(JK)*ZFALL
          ! Cloud budget diagnostic stored at end as implicit
          ! jl
        END IF
        ! LLFALL
      END DO
      ! jm

      !---------------------------------------------------------------
      ! Precip cover overlap using MAX-RAN Overlap
      ! Since precipitation is now prognostic we must
      !   1) apply an arbitrary minimum coverage (0.3) if precip>0
      !   2) abandon the 2-flux clr/cld treatment
      !   3) Thus, since we have no memory of the clear sky precip
      !      fraction, we mimic the previous method by reducing
      !      ZCOVPTOT(JL), which has the memory, proportionally with
      !      the precip evaporation rate, taking cloud fraction
      !      into account
      !   #3 above leads to much smoother vertical profiles of
      !   precipitation fraction than the Klein-Jakob scheme which
      !   monotonically increases precip fraction and then resets
      !   it to zero in a step function once clear-sky precip reaches
      !   zero.
      !---------------------------------------------------------------
      IF (ZQPRETOT(JK) > ZEPSEC) THEN
        ZCOVPTOT = 1.0_JPRB - ((1.0_JPRB - ZCOVPTOT)*(1.0_JPRB - MAX(ZA(JL, JK), ZA(JL, JK - 1)))) / (1.0_JPRB - MIN(ZA(JL, JK -  &
        & 1), 1.0_JPRB - 1.E-06_JPRB))
        ZCOVPTOT = MAX(ZCOVPTOT, YDECLDP%RCOVPMIN)
        ZCOVPCLR = MAX(0.0_JPRB, ZCOVPTOT - ZA(JL, JK))          ! clear sky proportion
        ZRAINCLD = ZQXFG(NCLDQR, JK) / ZCOVPTOT
        ZSNOWCLD = ZQXFG(NCLDQS, JK) / ZCOVPTOT
        ZCOVPMAX = MAX(ZCOVPTOT, ZCOVPMAX)
      ELSE
        ZRAINCLD = 0.0_JPRB
        ZSNOWCLD = 0.0_JPRB
        ZCOVPTOT = 0.0_JPRB          ! no flux - reset cover
        ZCOVPCLR = 0.0_JPRB          ! reset clear sky proportion
        ZCOVPMAX = 0.0_JPRB          ! reset max cover for ZZRH calc
      END IF

      !----------------------------------------------------------------------
      ! 4.3a AUTOCONVERSION TO SNOW
      !----------------------------------------------------------------------

      IF (ZTP1(JL, JK) <= YDCST%RTT) THEN
        !-----------------------------------------------------
        !     Snow Autoconversion rate follow Lin et al. 1983
        !-----------------------------------------------------
        IF (ZICECLD(JK) > ZEPSEC) THEN

          ZZCO = PTSPHY*YDECLDP%RSNOWLIN1*EXP(YDECLDP%RSNOWLIN2*(ZTP1(JL, JK) - YDCST%RTT))

          IF (YDECLDP%LAERICEAUTO) THEN
            ZLCRIT = PICRIT_AER(JL, JK)
            ! 0.3 = N**0.333 with N=0.027
            ZZCO = ZZCO*(YDECLDP%RNICE / PNICE(JL, JK))**0.333_JPRB
          ELSE
            ZLCRIT = YDECLDP%RLCRITSNOW
          END IF

          ZSNOWAUT = ZZCO*(1.0_JPRB - EXP(-(ZICECLD(JK) / ZLCRIT)**2))
          ZSOLQB(NCLDQS, NCLDQI, JK) = ZSOLQB(NCLDQS, NCLDQI, JK) + ZSNOWAUT

        END IF
      END IF

      !----------------------------------------------------------------------
      ! 4.3b AUTOCONVERSION WARM CLOUDS
      !   Collection and accretion will require separate treatment
      !   but for now we keep this simple treatment
      !----------------------------------------------------------------------

      IF (ZLIQCLD > ZEPSEC) THEN

        !--------------------------------------------------------
        !-
        !- Warm-rain process follow Sundqvist (1989)
        !-
        !--------------------------------------------------------
        IF (IWARMRAIN == 1) THEN

          ZZCO = YDECLDP%RKCONV*PTSPHY

          IF (YDECLDP%LAERLIQAUTOLSP) THEN
            ZLCRIT = PLCRIT_AER(JL, JK)
            ! 0.3 = N**0.333 with N=125 cm-3
            ZZCO = ZZCO*(YDECLDP%RCCN / PCCN(JL, JK))**0.333_JPRB
          ELSE
            ! Modify autoconversion threshold dependent on:
            !  land (polluted, high CCN, smaller droplets, higher threshold)
            !  sea  (clean, low CCN, larger droplets, lower threshold)
            IF (PLSM(JL) > 0.5_JPRB) THEN
              ZLCRIT = YDECLDP%RCLCRIT_LAND                ! land
            ELSE
              ZLCRIT = YDECLDP%RCLCRIT_SEA                ! ocean
            END IF
          END IF

          !------------------------------------------------------------------
          ! Parameters for cloud collection by rain and snow.
          ! Note that with new prognostic variable it is now possible
          ! to REPLACE this with an explicit collection parametrization
          !------------------------------------------------------------------
          ZPRECIP = (ZPFPLSX(JL, JK, NCLDQS) + ZPFPLSX(JL, JK, NCLDQR)) / MAX(ZEPSEC, ZCOVPTOT)
          ZCFPR = 1.0_JPRB + YDECLDP%RPRC1*SQRT(MAX(ZPRECIP, 0.0_JPRB))
          !      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))*&
          !       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))

          IF (YDECLDP%LAERLIQCOLL) THEN
            ! 5.0 = N**0.333 with N=125 cm-3
            ZCFPR = ZCFPR*(YDECLDP%RCCN / PCCN(JL, JK))**0.333_JPRB
          END IF

          ZZCO = ZZCO*ZCFPR
          ZLCRIT = ZLCRIT / MAX(ZCFPR, ZEPSEC)

          IF (ZLIQCLD / ZLCRIT < 20.0_JPRB) THEN
            ! Security for exp for some compilers
            ZRAINAUT = ZZCO*(1.0_JPRB - EXP(-(ZLIQCLD / ZLCRIT)**2))
          ELSE
            ZRAINAUT = ZZCO
          END IF

          ! rain freezes instantly
          IF (ZTP1(JL, JK) <= YDCST%RTT) THEN
            ZSOLQB(NCLDQS, NCLDQL, JK) = ZSOLQB(NCLDQS, NCLDQL, JK) + ZRAINAUT
          ELSE
            ZSOLQB(NCLDQR, NCLDQL, JK) = ZSOLQB(NCLDQR, NCLDQL, JK) + ZRAINAUT
          END IF

          !--------------------------------------------------------
          !-
          !- Warm-rain process follow Khairoutdinov and Kogan (2000)
          !-
          !--------------------------------------------------------
        ELSE IF (IWARMRAIN == 2) THEN

          IF (PLSM(JL) > 0.5_JPRB) THEN
            ! land
            ZCONST = YDECLDP%RCL_KK_CLOUD_NUM_LAND
            ZLCRIT = YDECLDP%RCLCRIT_LAND
          ELSE
            ! ocean
            ZCONST = YDECLDP%RCL_KK_CLOUD_NUM_SEA
            ZLCRIT = YDECLDP%RCLCRIT_SEA
          END IF

          IF (ZLIQCLD > ZLCRIT) THEN

            ZRAINAUT = 1.5_JPRB*ZA(JL, JK)*PTSPHY*YDECLDP%RCL_KKAAU*ZLIQCLD**YDECLDP%RCL_KKBAUQ*ZCONST**YDECLDP%RCL_KKBAUN

            ZRAINAUT = MIN(ZRAINAUT, ZQXFG(NCLDQL, JK))
            IF (ZRAINAUT < ZEPSEC)             ZRAINAUT = 0.0_JPRB

            ZRAINACC = 2.0_JPRB*ZA(JL, JK)*PTSPHY*YDECLDP%RCL_KKAAC*(ZLIQCLD*ZRAINCLD)**YDECLDP%RCL_KKBAC

            ZRAINACC = MIN(ZRAINACC, ZQXFG(NCLDQL, JK))
            IF (ZRAINACC < ZEPSEC)             ZRAINACC = 0.0_JPRB

          ELSE
            ZRAINAUT = 0.0_JPRB
            ZRAINACC = 0.0_JPRB
          END IF

          ! If temperature < 0, then autoconversion produces snow rather than rain
          ! Explicit
          IF (ZTP1(JL, JK) <= YDCST%RTT) THEN
            ZSOLQA(NCLDQS, NCLDQL, JK) = ZSOLQA(NCLDQS, NCLDQL, JK) + ZRAINAUT
            ZSOLQA(NCLDQS, NCLDQL, JK) = ZSOLQA(NCLDQS, NCLDQL, JK) + ZRAINACC
            ZSOLQA(NCLDQL, NCLDQS, JK) = ZSOLQA(NCLDQL, NCLDQS, JK) - ZRAINAUT
            ZSOLQA(NCLDQL, NCLDQS, JK) = ZSOLQA(NCLDQL, NCLDQS, JK) - ZRAINACC
          ELSE
            ZSOLQA(NCLDQR, NCLDQL, JK) = ZSOLQA(NCLDQR, NCLDQL, JK) + ZRAINAUT
            ZSOLQA(NCLDQR, NCLDQL, JK) = ZSOLQA(NCLDQR, NCLDQL, JK) + ZRAINACC
            ZSOLQA(NCLDQL, NCLDQR, JK) = ZSOLQA(NCLDQL, NCLDQR, JK) - ZRAINAUT
            ZSOLQA(NCLDQL, NCLDQR, JK) = ZSOLQA(NCLDQL, NCLDQR, JK) - ZRAINACC
          END IF

        END IF
        ! on IWARMRAIN

      END IF
      ! on ZLIQCLD > ZEPSEC


      !----------------------------------------------------------------------
      ! RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
      !      only active if T<0degC and supercooled liquid water is present
      !      AND if not Sundquist autoconversion (as this includes riming)
      !----------------------------------------------------------------------
      IF (IWARMRAIN > 1) THEN

        IF (ZTP1(JL, JK) <= YDCST%RTT .and. ZLIQCLD > ZEPSEC) THEN

          ! Fallspeed air density correction
          ZFALLCORR = (YDECLDP%RDENSREF / ZRHO(JK))**0.4_JPRB

          !------------------------------------------------------------------
          ! Riming of snow by cloud water - implicit in lwc
          !------------------------------------------------------------------
          IF (ZSNOWCLD > ZEPSEC .and. ZCOVPTOT > 0.01_JPRB) THEN

            ! Calculate riming term
            ! Factor of liq water taken out because implicit
            ZSNOWRIME =  &
            & 0.3_JPRB*ZCOVPTOT*PTSPHY*YDECLDP%RCL_CONST7S*ZFALLCORR*(ZRHO(JK)*ZSNOWCLD*YDECLDP%RCL_CONST1S)**YDECLDP%RCL_CONST8S

            ! Limit snow riming term
            ZSNOWRIME = MIN(ZSNOWRIME, 1.0_JPRB)

            ZSOLQB(NCLDQS, NCLDQL, JK) = ZSOLQB(NCLDQS, NCLDQL, JK) + ZSNOWRIME

          END IF

          !------------------------------------------------------------------
          ! Riming of ice by cloud water - implicit in lwc
          ! NOT YET ACTIVE
          !------------------------------------------------------------------
          !      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
          !
          !        ! Calculate riming term
          !        ! Factor of liq water taken out because implicit
          !        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
          !     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
          !
          !        ! Limit ice riming term
          !        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
          !
          !        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
          !
          !      ENDIF
        END IF

      END IF
      ! on IWARMRAIN > 1


      !----------------------------------------------------------------------
      ! 4.4a  MELTING OF SNOW and ICE
      !       with new implicit solver this also has to treat snow or ice
      !       precipitating from the level above... i.e. local ice AND flux.
      !       in situ ice and snow: could arise from LS advection or warming
      !       falling ice and snow: arrives by precipitation process
      !----------------------------------------------------------------------

      ZICETOT = ZQXFG(NCLDQI, JK) + ZQXFG(NCLDQS, JK)
      ZMELTMAX = 0.0_JPRB

      ! If there are frozen hydrometeors present and dry-bulb temperature > 0degC
      IF (ZICETOT > ZEPSEC .and. ZTP1(JL, JK) > YDCST%RTT) THEN

        ! Calculate subsaturation
        ZSUBSAT = MAX(ZQSICE(JL, JK) - ZQX(JL, JK, NCLDQV), 0.0_JPRB)

        ! Calculate difference between dry-bulb (ZTP1) and the temperature
        ! at which the wet-bulb=0degC (YDCST%RTT-ZSUBSAT*....) using an approx.
        ! Melting only occurs if the wet-bulb temperature >0
        ! i.e. warming of ice particle due to melting > cooling
        ! due to evaporation.
        ZTDMTW0 = ZTP1(JL, JK) - YDCST%RTT - ZSUBSAT*(ZTW1 + ZTW2*(PAP(JL, JK) - ZTW3) - ZTW4*(ZTP1(JL, JK) - ZTW5))
        ! Not implicit yet...
        ! Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
        ZCONS1 = ABS((PTSPHY*(1.0_JPRB + 0.5_JPRB*ZTDMTW0)) / YDECLDP%RTAUMEL)
        ZMELTMAX = MAX(ZTDMTW0*ZCONS1*ZRLDCP, 0.0_JPRB)
      END IF

      ! Loop over frozen hydrometeors (ice, snow)
      DO JM=1,NCLV
        IF (IPHASE(JM) == 2) THEN
          JN = IMELT(JM)
          IF (ZMELTMAX > ZEPSEC .and. ZICETOT > ZEPSEC) THEN
            ! Apply melting in same proportion as frozen hydrometeor fractions
            ZALFA = ZQXFG(JM, JK) / ZICETOT
            ZMELT = MIN(ZQXFG(JM, JK), ZALFA*ZMELTMAX)
            ! needed in first guess
            ! This implies that zqpretot has to be recalculated below
            ! since is not conserved here if ice falls and liquid doesn't
            ZQXFG(JM, JK) = ZQXFG(JM, JK) - ZMELT
            ZQXFG(JN, JK) = ZQXFG(JN, JK) + ZMELT
            ZSOLQA(JN, JM, JK) = ZSOLQA(JN, JM, JK) + ZMELT
            ZSOLQA(JM, JN, JK) = ZSOLQA(JM, JN, JK) - ZMELT
          END IF
        END IF
      END DO

      !----------------------------------------------------------------------
      ! 4.4b  FREEZING of RAIN
      !----------------------------------------------------------------------

      ! If rain present
      IF (ZQX(JL, JK, NCLDQR) > ZEPSEC) THEN

        IF (ZTP1(JL, JK) <= YDCST%RTT .and. ZTP1(JL, JK - 1) > YDCST%RTT) THEN
          ! Base of melting layer/top of refreezing layer so
          ! store rain/snow fraction for precip type diagnosis
          ! If mostly rain, then supercooled rain slow to freeze
          ! otherwise faster to freeze (snow or ice pellets)
          ZQPRETOT(JK) = MAX(ZQX(JL, JK, NCLDQS) + ZQX(JL, JK, NCLDQR), ZEPSEC)
          PRAINFRAC_TOPRFZ(JL) = ZQX(JL, JK, NCLDQR) / ZQPRETOT(JK)
          IF (PRAINFRAC_TOPRFZ(JL) > 0.8) THEN
            LLRAINLIQ = .true.
          ELSE
            LLRAINLIQ = .false.
          END IF
        END IF

        ! If temperature less than zero
        IF (ZTP1(JL, JK) < YDCST%RTT) THEN

          IF (PRAINFRAC_TOPRFZ(JL) > 0.8) THEN

            ! Majority of raindrops completely melted
            ! Refreezing is by slow heterogeneous freezing

            ! Slope of rain particle size distribution
            ZLAMBDA = (YDECLDP%RCL_FAC1 / (ZRHO(JK)*ZQX(JL, JK, NCLDQR)))**YDECLDP%RCL_FAC2

            ! Calculate freezing rate based on Bigg(1953) and Wisner(1972)
            ZTEMP = YDECLDP%RCL_FZRAB*(ZTP1(JL, JK) - YDCST%RTT)
            ZFRZ = PTSPHY*(YDECLDP%RCL_CONST5R / ZRHO(JK))*(EXP(ZTEMP) - 1._JPRB)*ZLAMBDA**YDECLDP%RCL_CONST6R
            ZFRZMAX = MAX(ZFRZ, 0.0_JPRB)

          ELSE

            ! Majority of raindrops only partially melted
            ! Refreeze with a shorter timescale (reverse of melting...for now)

            ZCONS1 = ABS((PTSPHY*(1.0_JPRB + 0.5_JPRB*(YDCST%RTT - ZTP1(JL, JK)))) / YDECLDP%RTAUMEL)
            ZFRZMAX = MAX((YDCST%RTT - ZTP1(JL, JK))*ZCONS1*ZRLDCP, 0.0_JPRB)

          END IF

          IF (ZFRZMAX > ZEPSEC) THEN
            ZFRZ = MIN(ZQX(JL, JK, NCLDQR), ZFRZMAX)
            ZSOLQA(NCLDQS, NCLDQR, JK) = ZSOLQA(NCLDQS, NCLDQR, JK) + ZFRZ
            ZSOLQA(NCLDQR, NCLDQS, JK) = ZSOLQA(NCLDQR, NCLDQS, JK) - ZFRZ
          END IF
        END IF

      END IF


      !----------------------------------------------------------------------
      ! 4.4c  FREEZING of LIQUID
      !----------------------------------------------------------------------
      ! not implicit yet...
      ZFRZMAX = MAX((YDECLDP%RTHOMO - ZTP1(JL, JK))*ZRLDCP, 0.0_JPRB)

      JM = NCLDQL
      JN = IMELT(JM)
      IF (ZFRZMAX > ZEPSEC .and. ZQXFG(JM, JK) > ZEPSEC) THEN
        ZFRZ = MIN(ZQXFG(JM, JK), ZFRZMAX)
        ZSOLQA(JN, JM, JK) = ZSOLQA(JN, JM, JK) + ZFRZ
        ZSOLQA(JM, JN, JK) = ZSOLQA(JM, JN, JK) - ZFRZ
      END IF

      !----------------------------------------------------------------------
      ! 4.5   EVAPORATION OF RAIN/SNOW
      !----------------------------------------------------------------------

      !----------------------------------------
      ! Rain evaporation scheme from Sundquist
      !----------------------------------------
      IF (IEVAPRAIN == 1) THEN

        ! Rain


        ZZRH = YDECLDP%RPRECRHMAX + ((1.0_JPRB - YDECLDP%RPRECRHMAX)*ZCOVPMAX) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZZRH = MIN(MAX(ZZRH, YDECLDP%RPRECRHMAX), 1.0_JPRB)

        ZQE = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSLIQ(JL, JK)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        !---------------------------------------------
        ! humidity in moistest ZCOVPCLR part of domain
        !---------------------------------------------
        ZQE = MAX(0.0_JPRB, MIN(ZQE, ZQSLIQ(JL, JK)))
        LLO1 = ZCOVPCLR > ZEPSEC .and. ZQXFG(NCLDQR, JK) > ZEPSEC .and. ZQE < ZZRH*ZQSLIQ(JL, JK)

        IF (LLO1) THEN
          ! note: zpreclr is a rain flux
          ZPRECLR = (ZQXFG(NCLDQR, JK)*ZCOVPCLR) / SIGN(MAX(ABS(ZCOVPTOT*ZDTGDP(JK)), ZEPSILON), ZCOVPTOT*ZDTGDP(JK))

          !--------------------------------------
          ! actual microphysics formula in zbeta
          !--------------------------------------

          ZBETA1 = ((SQRT(PAP(JL, JK) / PAPH(JL, KLEV + 1)) / YDECLDP%RVRFACTOR)*ZPRECLR) / MAX(ZCOVPCLR, ZEPSEC)

          ZBETA = YDCST%RG*YDECLDP%RPECONS*0.5_JPRB*ZBETA1**0.5777_JPRB

          ZDENOM = 1.0_JPRB + ZBETA*PTSPHY*ZCORQSLIQ(JK)
          ZDPR = ((ZCOVPCLR*ZBETA*(ZQSLIQ(JL, JK) - ZQE)) / ZDENOM)*ZDP(JK)*ZRG_R
          ZDPEVAP = ZDPR*ZDTGDP(JK)

          !---------------------------------------------------------
          ! add evaporation term to explicit sink.
          ! this has to be explicit since if treated in the implicit
          ! term evaporation can not reduce rain to zero and model
          ! produces small amounts of rainfall everywhere.
          !---------------------------------------------------------

          ! Evaporate rain
          ZEVAP = MIN(ZDPEVAP, ZQXFG(NCLDQR, JK))

          ZSOLQA(NCLDQV, NCLDQR, JK) = ZSOLQA(NCLDQV, NCLDQR, JK) + ZEVAP
          ZSOLQA(NCLDQR, NCLDQV, JK) = ZSOLQA(NCLDQR, NCLDQV, JK) - ZEVAP

          !-------------------------------------------------------------
          ! Reduce the total precip coverage proportional to evaporation
          ! to mimic the previous scheme which had a diagnostic
          ! 2-flux treatment, abandoned due to the new prognostic precip
          !-------------------------------------------------------------
          ZCOVPTOT = MAX(YDECLDP%RCOVPMIN, ZCOVPTOT - MAX(0.0_JPRB, ((ZCOVPTOT - ZA(JL, JK))*ZEVAP) / ZQXFG(NCLDQR, JK)))

          ! Update fg field
          ZQXFG(NCLDQR, JK) = ZQXFG(NCLDQR, JK) - ZEVAP

        END IF


        !---------------------------------------------------------
        ! Rain evaporation scheme based on Abel and Boutle (2013)
        !---------------------------------------------------------
      ELSE IF (IEVAPRAIN == 2) THEN


        !-----------------------------------------------------------------------
        ! Calculate relative humidity limit for rain evaporation
        ! to avoid cloud formation and saturation of the grid box
        !-----------------------------------------------------------------------
        ! Limit RH for rain evaporation dependent on precipitation fraction
        ZZRH = YDECLDP%RPRECRHMAX + ((1.0_JPRB - YDECLDP%RPRECRHMAX)*ZCOVPMAX) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZZRH = MIN(MAX(ZZRH, YDECLDP%RPRECRHMAX), 1.0_JPRB)

        ! Critical relative humidity
        !ZRHC=RAMID
        !ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
        ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
        !IF(ZSIGK > 0.8_JPRB) THEN
        !  ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
        !ENDIF
        !ZZRH = MIN(ZRHC,ZZRH)

        ! Further limit RH for rain evaporation to 80% (RHcrit in free troposphere)
        ZZRH = MIN(0.8_JPRB, ZZRH)

        ZQE = MAX(0.0_JPRB, MIN(ZQX(JL, JK, NCLDQV), ZQSLIQ(JL, JK)))

        LLO1 = ZCOVPCLR > ZEPSEC .and. ZQXFG(NCLDQR, JK) > ZEPSEC .and. ZQE < ZZRH*ZQSLIQ(JL, JK)

        IF (LLO1) THEN

          !-------------------------------------------
          ! Abel and Boutle (2012) evaporation
          !-------------------------------------------
          ! Calculate local precipitation (kg/kg)
          ZPRECLR = ZQXFG(NCLDQR, JK) / ZCOVPTOT

          ! Fallspeed air density correction
          ZFALLCORR = (YDECLDP%RDENSREF / ZRHO(JK))**0.4

          ! Saturation vapour pressure with respect to liquid phase
          ZESATLIQ = (YDCST%RV / YDCST%RD)*FOEELIQ(ZTP1(JL, JK))

          ! Slope of particle size distribution
          ZLAMBDA = (YDECLDP%RCL_FAC1 / (ZRHO(JK)*ZPRECLR))**YDECLDP%RCL_FAC2            ! ZPRECLR=kg/kg

          ZEVAP_DENOM = YDECLDP%RCL_CDENOM1*ZESATLIQ - YDECLDP%RCL_CDENOM2*ZTP1(JL, JK)*ZESATLIQ + YDECLDP%RCL_CDENOM3*ZTP1(JL,  &
          & JK)**3._JPRB*PAP(JL, JK)

          ! Temperature dependent conductivity
          ZCORR2 = ((ZTP1(JL, JK) / 273._JPRB)**1.5_JPRB*393._JPRB) / (ZTP1(JL, JK) + 120._JPRB)
          ZKA = YDECLDP%RCL_KA273*ZCORR2

          ZSUBSAT = MAX(ZZRH*ZQSLIQ(JL, JK) - ZQE, 0.0_JPRB)

          ZBETA = (0.5_JPRB / ZQSLIQ(JL, JK))*ZTP1(JL, JK)**2._JPRB*ZESATLIQ*YDECLDP%RCL_CONST1R*(ZCORR2 / ZEVAP_DENOM) &
          & *(0.78_JPRB / (ZLAMBDA**YDECLDP%RCL_CONST4R) + (YDECLDP%RCL_CONST2R*(ZRHO(JK)*ZFALLCORR)**0.5_JPRB) /  &
          & (ZCORR2**0.5_JPRB*ZLAMBDA**YDECLDP%RCL_CONST3R))

          ZDENOM = 1.0_JPRB + ZBETA*PTSPHY            !*ZCORQSLIQ(JL)
          ZDPEVAP = (ZCOVPCLR*ZBETA*PTSPHY*ZSUBSAT) / ZDENOM

          !---------------------------------------------------------
          ! Add evaporation term to explicit sink.
          ! this has to be explicit since if treated in the implicit
          ! term evaporation can not reduce rain to zero and model
          ! produces small amounts of rainfall everywhere.
          !---------------------------------------------------------

          ! Limit rain evaporation
          ZEVAP = MIN(ZDPEVAP, ZQXFG(NCLDQR, JK))

          ZSOLQA(NCLDQV, NCLDQR, JK) = ZSOLQA(NCLDQV, NCLDQR, JK) + ZEVAP
          ZSOLQA(NCLDQR, NCLDQV, JK) = ZSOLQA(NCLDQR, NCLDQV, JK) - ZEVAP

          !-------------------------------------------------------------
          ! Reduce the total precip coverage proportional to evaporation
          ! to mimic the previous scheme which had a diagnostic
          ! 2-flux treatment, abandoned due to the new prognostic precip
          !-------------------------------------------------------------
          ZCOVPTOT = MAX(YDECLDP%RCOVPMIN, ZCOVPTOT - MAX(0.0_JPRB, ((ZCOVPTOT - ZA(JL, JK))*ZEVAP) / ZQXFG(NCLDQR, JK)))

          ! Update fg field
          ZQXFG(NCLDQR, JK) = ZQXFG(NCLDQR, JK) - ZEVAP

        END IF

      END IF
      ! on IEVAPRAIN

      !----------------------------------------------------------------------
      ! 4.5   EVAPORATION OF SNOW
      !----------------------------------------------------------------------
      ! Snow
      IF (IEVAPSNOW == 1) THEN

        ZZRH = YDECLDP%RPRECRHMAX + ((1.0_JPRB - YDECLDP%RPRECRHMAX)*ZCOVPMAX) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZZRH = MIN(MAX(ZZRH, YDECLDP%RPRECRHMAX), 1.0_JPRB)
        ZQE = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSICE(JL, JK)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))

        !---------------------------------------------
        ! humidity in moistest ZCOVPCLR part of domain
        !---------------------------------------------
        ZQE = MAX(0.0_JPRB, MIN(ZQE, ZQSICE(JL, JK)))
        LLO1 = ZCOVPCLR > ZEPSEC .and. ZQXFG(NCLDQS, JK) > ZEPSEC .and. ZQE < ZZRH*ZQSICE(JL, JK)

        IF (LLO1) THEN
          ! note: zpreclr is a rain flux a
          ZPRECLR = (ZQXFG(NCLDQS, JK)*ZCOVPCLR) / SIGN(MAX(ABS(ZCOVPTOT*ZDTGDP(JK)), ZEPSILON), ZCOVPTOT*ZDTGDP(JK))

          !--------------------------------------
          ! actual microphysics formula in zbeta
          !--------------------------------------

          ZBETA1 = ((SQRT(PAP(JL, JK) / PAPH(JL, KLEV + 1)) / YDECLDP%RVRFACTOR)*ZPRECLR) / MAX(ZCOVPCLR, ZEPSEC)

          ZBETA = YDCST%RG*YDECLDP%RPECONS*ZBETA1**0.5777_JPRB

          ZDENOM = 1.0_JPRB + ZBETA*PTSPHY*ZCORQSICE(JK)
          ZDPR = ((ZCOVPCLR*ZBETA*(ZQSICE(JL, JK) - ZQE)) / ZDENOM)*ZDP(JK)*ZRG_R
          ZDPEVAP = ZDPR*ZDTGDP(JK)

          !---------------------------------------------------------
          ! add evaporation term to explicit sink.
          ! this has to be explicit since if treated in the implicit
          ! term evaporation can not reduce snow to zero and model
          ! produces small amounts of snowfall everywhere.
          !---------------------------------------------------------

          ! Evaporate snow
          ZEVAP = MIN(ZDPEVAP, ZQXFG(NCLDQS, JK))

          ZSOLQA(NCLDQV, NCLDQS, JK) = ZSOLQA(NCLDQV, NCLDQS, JK) + ZEVAP
          ZSOLQA(NCLDQS, NCLDQV, JK) = ZSOLQA(NCLDQS, NCLDQV, JK) - ZEVAP

          !-------------------------------------------------------------
          ! Reduce the total precip coverage proportional to evaporation
          ! to mimic the previous scheme which had a diagnostic
          ! 2-flux treatment, abandoned due to the new prognostic precip
          !-------------------------------------------------------------
          ZCOVPTOT = MAX(YDECLDP%RCOVPMIN, ZCOVPTOT - MAX(0.0_JPRB, ((ZCOVPTOT - ZA(JL, JK))*ZEVAP) / ZQXFG(NCLDQS, JK)))

          !Update first guess field
          ZQXFG(NCLDQS, JK) = ZQXFG(NCLDQS, JK) - ZEVAP

        END IF
        !---------------------------------------------------------
      ELSE IF (IEVAPSNOW == 2) THEN



        !-----------------------------------------------------------------------
        ! Calculate relative humidity limit for snow evaporation
        !-----------------------------------------------------------------------
        ZZRH = YDECLDP%RPRECRHMAX + ((1.0_JPRB - YDECLDP%RPRECRHMAX)*ZCOVPMAX) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZZRH = MIN(MAX(ZZRH, YDECLDP%RPRECRHMAX), 1.0_JPRB)
        ZQE = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSICE(JL, JK)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))

        !---------------------------------------------
        ! humidity in moistest ZCOVPCLR part of domain
        !---------------------------------------------
        ZQE = MAX(0.0_JPRB, MIN(ZQE, ZQSICE(JL, JK)))
        LLO1 = ZCOVPCLR > ZEPSEC .and. ZQX(JL, JK, NCLDQS) > ZEPSEC .and. ZQE < ZZRH*ZQSICE(JL, JK)

        IF (LLO1) THEN

          ! Calculate local precipitation (kg/kg)
          ZPRECLR = ZQX(JL, JK, NCLDQS) / ZCOVPTOT
          ZVPICE = (FOEEICE(ZTP1(JL, JK))*YDCST%RV) / YDCST%RD

          ! Particle size distribution
          ! ZTCG increases Ni with colder temperatures - essentially a
          ! Fletcher or Meyers scheme?
          ZTCG = 1.0_JPRB            !v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
          ! ZFACX1I modification is based on Andrew Barrett's results
          ZFACX1S = 1.0_JPRB            !v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB

          ZAPLUSB = YDECLDP%RCL_APB1*ZVPICE - YDECLDP%RCL_APB2*ZVPICE*ZTP1(JL, JK) + PAP(JL, JK)*YDECLDP%RCL_APB3*ZTP1(JL, JK)**3
          ZCORRFAC = (1.0 / ZRHO(JK))**0.5
          ZCORRFAC2 = ((ZTP1(JL, JK) / 273.0)**1.5)*(393.0 / (ZTP1(JL, JK) + 120.0))

          ZPR02 = (ZRHO(JK)*ZPRECLR*YDECLDP%RCL_CONST1S) / (ZTCG*ZFACX1S)

          ZTERM1 = ((ZQSICE(JL, JK) - ZQE)*ZTP1(JL, JK)**2*ZVPICE*ZCORRFAC2*ZTCG*YDECLDP%RCL_CONST2S*ZFACX1S) / (ZRHO(JK) &
          & *ZAPLUSB*ZQSICE(JL, JK))
          ZTERM2 = 0.65*YDECLDP%RCL_CONST6S*ZPR02**YDECLDP%RCL_CONST4S + (YDECLDP%RCL_CONST3S*ZCORRFAC**0.5*ZRHO(JK) &
          & **0.5*ZPR02**YDECLDP%RCL_CONST5S) / ZCORRFAC2**0.5

          ZDPEVAP = MAX(ZCOVPCLR*ZTERM1*ZTERM2*PTSPHY, 0.0_JPRB)

          !--------------------------------------------------------------------
          ! Limit evaporation to snow amount
          !--------------------------------------------------------------------
          ZEVAP = MIN(ZDPEVAP, ZEVAPLIMICE(JK))
          ZEVAP = MIN(ZEVAP, ZQX(JL, JK, NCLDQS))


          ZSOLQA(NCLDQV, NCLDQS, JK) = ZSOLQA(NCLDQV, NCLDQS, JK) + ZEVAP
          ZSOLQA(NCLDQS, NCLDQV, JK) = ZSOLQA(NCLDQS, NCLDQV, JK) - ZEVAP

          !-------------------------------------------------------------
          ! Reduce the total precip coverage proportional to evaporation
          ! to mimic the previous scheme which had a diagnostic
          ! 2-flux treatment, abandoned due to the new prognostic precip
          !-------------------------------------------------------------
          ZCOVPTOT = MAX(YDECLDP%RCOVPMIN, ZCOVPTOT - MAX(0.0_JPRB, ((ZCOVPTOT - ZA(JL, JK))*ZEVAP) / ZQX(JL, JK, NCLDQS)))

          !Update first guess field
          ZQXFG(NCLDQS, JK) = ZQXFG(NCLDQS, JK) - ZEVAP

        END IF

      END IF
      ! on IEVAPSNOW

      !--------------------------------------
      ! Evaporate small precipitation amounts
      !--------------------------------------
      DO JM=1,NCLV
        IF (LLFALL(JM)) THEN
          IF (ZQXFG(JM, JK) < YDECLDP%RLMIN) THEN
            ZSOLQA(NCLDQV, JM, JK) = ZSOLQA(NCLDQV, JM, JK) + ZQXFG(JM, JK)
            ZSOLQA(JM, NCLDQV, JK) = ZSOLQA(JM, NCLDQV, JK) - ZQXFG(JM, JK)
          END IF
        END IF
      END DO

      !######################################################################
      !            5.0  *** SOLVERS FOR A AND L ***
      ! now use an implicit solution rather than exact solution
      ! solver is forward in time, upstream difference for advection
      !######################################################################

      !---------------------------
      ! 5.1 solver for cloud cover
      !---------------------------
      ZANEW = (ZA(JL, JK) + ZSOLAC(JK)) / (1.0_JPRB + ZSOLAB(JK))
      ZANEW = MIN(ZANEW, 1.0_JPRB)
      IF (ZANEW < YDECLDP%RAMIN)       ZANEW = 0.0_JPRB
      ZDA(JK) = ZANEW - ZAORIG(JL, JK)
      !---------------------------------
      ! variables needed for next level
      !---------------------------------
      ZANEWM1 = ZANEW

      !--------------------------------
      ! 5.2 solver for the microphysics
      !--------------------------------

      !--------------------------------------------------------------
      ! Truncate explicit sinks to avoid negatives
      ! Note: Species are treated in the order in which they run out
      ! since the clipping will alter the balance for the other vars
      !--------------------------------------------------------------

      DO JM=1,NCLV
!$claw nodep
        DO JN=1,NCLV
          LLINDEX3(JN, JM) = .false.
        END DO
        ZSINKSUM(JM) = 0.0_JPRB
      END DO

      !----------------------------
      ! collect sink terms and mark
      !----------------------------
      DO JM=1,NCLV
        DO JN=1,NCLV
          ZSINKSUM(JM) = ZSINKSUM(JM) - ZSOLQA(JM, JN, JK)            ! +ve total is bad
        END DO
      END DO

      !---------------------------------------
      ! calculate overshoot and scaling factor
      !---------------------------------------
      DO JM=1,NCLV
        ZMAX = MAX(ZQX(JL, JK, JM), ZEPSEC)
        ZRAT = MAX(ZSINKSUM(JM), ZMAX)
        ZRATIO(JM) = ZMAX / ZRAT
      END DO

      !--------------------------------------------
      ! scale the sink terms, in the correct order,
      ! recalculating the scale factor each time
      !--------------------------------------------
      DO JM=1,NCLV
        ZSINKSUM(JM) = 0.0_JPRB
      END DO

      !----------------
      ! recalculate sum
      !----------------
      DO JM=1,NCLV
        PSUM_SOLQA = 0.0
        DO JN=1,NCLV
          PSUM_SOLQA = PSUM_SOLQA + ZSOLQA(JM, JN, JK)
        END DO
        ! ZSINKSUM(JL,JM)=ZSINKSUM(JL,JM)-SUM(ZSOLQA(JL,JM,1:NCLV))
        ZSINKSUM(JM) = ZSINKSUM(JM) - PSUM_SOLQA
        !---------------------------
        ! recalculate scaling factor
        !---------------------------
        ZMM = MAX(ZQX(JL, JK, JM), ZEPSEC)
        ZRR = MAX(ZSINKSUM(JM), ZMM)
        ZRATIO(JM) = ZMM / ZRR
        !------
        ! scale
        !------
        ZZRATIO = ZRATIO(JM)
        !DIR$ IVDEP
        !DIR$ PREFERVECTOR
        DO JN=1,NCLV
          IF (ZSOLQA(JM, JN, JK) < 0.0_JPRB) THEN
            ZSOLQA(JM, JN, JK) = ZSOLQA(JM, JN, JK)*ZZRATIO
            ZSOLQA(JN, JM, JK) = ZSOLQA(JN, JM, JK)*ZZRATIO
          END IF
        END DO
      END DO

      !--------------------------------------------------------------
      ! 5.2.2 Solver
      !------------------------

      !------------------------
      ! set the LHS of equation
      !------------------------
      DO JM=1,NCLV
        DO JN=1,NCLV
          !----------------------------------------------
          ! diagonals: microphysical sink terms+transport
          !----------------------------------------------
          IF (JN == JM) THEN
            ZQLHS(JN, JM) = 1.0_JPRB + ZFALLSINK(JM, JK)
            DO JO=1,NCLV
              ZQLHS(JN, JM) = ZQLHS(JN, JM) + ZSOLQB(JO, JN, JK)
            END DO
            !------------------------------------------
            ! non-diagonals: microphysical source terms
            !------------------------------------------
          ELSE
            ZQLHS(JN, JM) = -ZSOLQB(JN, JM, JK)              ! here is the delta T - missing from doc.
          END IF
        END DO
      END DO

      !------------------------
      ! set the RHS of equation
      !------------------------
      DO JM=1,NCLV
        !---------------------------------
        ! sum the explicit source and sink
        !---------------------------------
        ZEXPLICIT = 0.0_JPRB
        DO JN=1,NCLV
          ZEXPLICIT = ZEXPLICIT + ZSOLQA(JM, JN, JK)            ! sum over middle index
        END DO
        ZQXN(JM, JK) = ZQX(JL, JK, JM) + ZEXPLICIT
      END DO

      !-----------------------------------
      ! *** solve by LU decomposition: ***
      !-----------------------------------

      ! Note: This fast way of solving NCLVxNCLV system
      !       assumes a good behaviour (i.e. non-zero diagonal
      !       terms with comparable orders) of the matrix stored
      !       in ZQLHS. For the moment this is the case but
      !       be aware to preserve it when doing eventual
      !       modifications.

      ! Non pivoting recursive factorization
      DO JN=1,NCLV - 1
        ! number of steps
        DO JM=JN + 1,NCLV
          ! row index
          ZQLHS(JM, JN) = ZQLHS(JM, JN) / ZQLHS(JN, JN)
          DO IK=JN + 1,NCLV
            ! column index
            ZQLHS(JM, IK) = ZQLHS(JM, IK) - ZQLHS(JM, JN)*ZQLHS(JN, IK)
          END DO
        END DO
      END DO

      ! Backsubstitution
      !  step 1
      DO JN=2,NCLV
        DO JM=1,JN - 1
          ZQXN(JN, JK) = ZQXN(JN, JK) - ZQLHS(JN, JM)*ZQXN(JM, JK)
        END DO
      END DO
      !  step 2
      ZQXN(NCLV, JK) = ZQXN(NCLV, JK) / ZQLHS(NCLV, NCLV)
      DO JN=NCLV - 1,1,-1
        DO JM=JN + 1,NCLV
          ZQXN(JN, JK) = ZQXN(JN, JK) - ZQLHS(JN, JM)*ZQXN(JM, JK)
        END DO
        ZQXN(JN, JK) = ZQXN(JN, JK) / ZQLHS(JN, JN)
      END DO

      ! Ensure no small values (including negatives) remain in cloud variables nor
      ! precipitation rates.
      ! Evaporate l,i,r,s to water vapour. Latent heating taken into account below
      DO JN=1,NCLV - 1
        IF (ZQXN(JN, JK) < ZEPSEC) THEN
          ZQXN(NCLDQV, JK) = ZQXN(NCLDQV, JK) + ZQXN(JN, JK)
          ZQXN(JN, JK) = 0.0_JPRB
        END IF
      END DO

      !--------------------------------
      ! variables needed for next level
      !--------------------------------
      DO JM=1,NCLV
        ZQXNM1(JM) = ZQXN(JM, JK)
        ZQXN2D(JL, JK, JM) = ZQXN(JM, JK)
      END DO

      !------------------------------------------------------------------------
      ! 5.3 Precipitation/sedimentation fluxes to next level
      !     diagnostic precipitation fluxes
      !     It is this scaled flux that must be used for source to next layer
      !------------------------------------------------------------------------

      DO JM=1,NCLV
        ZPFPLSX(JL, JK + 1, JM) = ZFALLSINK(JM, JK)*ZQXN(JM, JK)*ZRDTGDP(JK)
      END DO

      ! Ensure precipitation fraction is zero if no precipitation
      ZQPRETOT(JK) = ZPFPLSX(JL, JK + 1, NCLDQS) + ZPFPLSX(JL, JK + 1, NCLDQR)
      IF (ZQPRETOT(JK) < ZEPSEC) THEN
        ZCOVPTOT = 0.0_JPRB
      END IF

      ! Loki region-hoist group(pcovptot)
      !--------------------------------------------------
      ! Copy precipitation fraction into output variable
      !-------------------------------------------------
      PCOVPTOT(JL, JK) = ZCOVPTOT
      ! Loki end region-hoist

    END DO
    ! Loki - loop-fission
    DO JK=YDECLDP%NCLDTOP,KLEV

      CALL section6(KLON, KLEV, YDTHF, ZQTMST, JK, ZQXN, ZPSUPSATSRCE, ZQX0, ZCONVSINK, ZDA, JL, IPHASE, ZQX, ZCONVSRCE, ZFALLSRCE,  &
      & ZFALLSINK, TENDENCY_LOC_T, TENDENCY_LOC_a, TENDENCY_LOC_q, TENDENCY_LOC_CLD, ZFLUXQ)
    END DO

    CALL section8(KLON, ZFOEALFA, ZRG_R, PAPH, PTSPHY, ZQTMST, PVFL, ZQX0, ZLNEG, PVFI, ZQXN2D, YDCST, JL, ZPFPLSX, PLUDE, KLEV,  &
    & PFCQSNG, PFCQLNG, PFPLSN, PFCQRNG, PFCQNNG, PFSQIF, PFSQLTUR, ZGDPH_R, PFSQSF, PFSQRF, PFHPSN, PFSQITUR, PFPLSL, PFSQLF,  &
    & PFHPSL, ZALFAW)

    !===============================================================================
    !IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)
  END SUBROUTINE CLOUDSC_SCC_HOIST
  SUBROUTINE section3p4 (KLON, YDTHF, ZICEFRAC, PAPH, ZFOKOOP, YDECLDP, PAP, ZQSICE, PTSPHY, ZQTMST, JK, ZLICLD, ZEVAPLIMMIX, ZLI,  &
  & ZDP, JL, YDCST, PMFU, PVERVEL, ZRDCP, ZLDEFR, PHRSW, PMFD, ZQX, KLEV, ZA, ZEPSEC, PHRLW, ZLIQFRAC, ZQSMIX, ZSOLQA, ZSOLAC,  &
  & ZQXFG, ZTP1, ZLCONDLIM, ZCOR, ZQSAT, ZLEVAPI, ZMFDN, ZQP, ZCDMAX, ZRHC, ZDTFORC, ZDTDIAB, ZLCOND2, ZACOND, ZZDL, ZDQS,  &
  & ZDTDP, ZLEVAPL, ZLEVAP, ZZZDT, ZWTOT, ZCOND1, ZQE, ZDPMXDT, ZLCOND1, ZTOLD, ZSIGK, LLFLAG, ZQOLD, ZFAC, ZCOND)
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMPHYDER, ONLY: state_type
    USE YOECLDP, ONLY: NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
    USE YOECLDP, ONLY: TECLDP
    USE YOEPHLI, ONLY: TEPHLI
    USE YOMCST, ONLY: TOMCST
    USE YOETHF, ONLY: TOETHF
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! Number of grid points
    REAL(KIND=JPRB), INTENT(INOUT) :: ZSOLAC(KLEV)
    TYPE(TOETHF), INTENT(IN) :: YDTHF
    REAL(KIND=JPRB), INTENT(OUT) :: ZLCONDLIM
    REAL(KIND=JPRB), INTENT(OUT) :: ZCOR
    REAL(KIND=JPRB), INTENT(OUT) :: ZQSAT
    REAL(KIND=JPRB), INTENT(IN) :: ZLIQFRAC(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZQX(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZCDMAX
    REAL(KIND=JPRB), INTENT(OUT) :: ZLEVAPI
    REAL(KIND=JPRB), INTENT(OUT) :: ZQP
    REAL(KIND=JPRB), INTENT(OUT) :: ZMFDN
    TYPE(TECLDP), INTENT(IN) :: YDECLDP
    REAL(KIND=JPRB), INTENT(OUT) :: ZRHC
    REAL(KIND=JPRB), INTENT(IN) :: PMFU(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZDP(KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
    REAL(KIND=JPRB), INTENT(IN) :: PHRSW(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZDTFORC
    REAL(KIND=JPRB), INTENT(IN) :: ZQTMST
    REAL(KIND=JPRB), INTENT(OUT) :: ZDTDIAB
    INTEGER(KIND=JPIM), INTENT(IN) :: JK
    REAL(KIND=JPRB), INTENT(IN) :: PHRLW(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZLCOND2
    REAL(KIND=JPRB), INTENT(OUT) :: ZACOND
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQXFG(NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZZDL
    INTEGER(KIND=JPIM), INTENT(IN) :: JL
    TYPE(TOMCST), INTENT(IN) :: YDCST
    REAL(KIND=JPRB), INTENT(OUT) :: ZDQS
    REAL(KIND=JPRB), INTENT(IN) :: ZRDCP
    REAL(KIND=JPRB), INTENT(OUT) :: ZDTDP
    REAL(KIND=JPRB), INTENT(OUT) :: ZLEVAPL
    REAL(KIND=JPRB), INTENT(IN) :: ZA(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZLEVAP
    REAL(KIND=JPRB), INTENT(IN) :: ZLICLD(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZZZDT
    REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZWTOT
    REAL(KIND=JPRB), INTENT(IN) :: ZQSICE(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZFOKOOP(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZCOND1
    REAL(KIND=JPRB), INTENT(IN) :: ZEVAPLIMMIX(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZQE
    REAL(KIND=JPRB), INTENT(OUT) :: ZDPMXDT
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQSMIX(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZLCOND1
    REAL(KIND=JPRB), INTENT(OUT) :: ZTOLD
    REAL(KIND=JPRB), INTENT(IN) :: ZEPSEC
    REAL(KIND=JPRB), INTENT(IN) :: ZICEFRAC(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZSIGK
    REAL(KIND=JPRB), INTENT(IN) :: ZLI(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: ZSOLQA(NCLV, NCLV, KLEV)
    LOGICAL, INTENT(OUT) :: LLFLAG
    REAL(KIND=JPRB), INTENT(OUT) :: ZQOLD
    REAL(KIND=JPRB), INTENT(OUT) :: ZFAC
    REAL(KIND=JPRB), INTENT(INOUT) :: ZTP1(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZLDEFR(KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PMFD(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZCOND
#include "fcttre.ycst.h"
#include "fccld.ydthf.h"
!$acc routine seq
    ! Loki region-hoist group(section3p4)

    !----------------------------------------------------------------------
    ! 3.4  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
    !----------------------------------------------------------------------
    !  calculate dqs/dt
    !  Note: For the separate prognostic Qi and Ql, one would ideally use
    !  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
    !  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
    !  These would then instantaneous freeze if T<-38C or lead to ice growth
    !  by deposition in warmer mixed phase clouds.  However, since we do
    !  not have a separate prognostic equation for in-cloud humidity or a
    !  statistical scheme approach in place, the depositional growth of ice
    !  in the mixed phase can not be modelled and we resort to supersaturation
    !  wrt ice instanteously converting to ice over one timestep
    !  (see Tompkins et al. QJRMS 2007 for details)
    !  Thus for the initial implementation the diagnostic mixed phase is
    !  retained for the moment, and the level of approximation noted.
    !----------------------------------------------------------------------

    ZDTDP = (ZRDCP*ZTP1(JL, JK)) / PAP(JL, JK)
    ZDPMXDT = ZDP(JK)*ZQTMST
    ZMFDN = 0.0_JPRB
    IF (JK < KLEV)     ZMFDN = PMFU(JL, JK + 1) + PMFD(JL, JK + 1)
    ZWTOT = PVERVEL(JL, JK) + 0.5_JPRB*YDCST%RG*(PMFU(JL, JK) + PMFD(JL, JK) + ZMFDN)
    ZWTOT = MIN(ZDPMXDT, MAX(-ZDPMXDT, ZWTOT))
    ZZZDT = PHRSW(JL, JK) + PHRLW(JL, JK)
    ZDTDIAB = MIN(ZDPMXDT*ZDTDP, MAX(-ZDPMXDT*ZDTDP, ZZZDT))*PTSPHY + YDTHF%RALFDCP*ZLDEFR(JK)
    ! Note: ZLDEFR should be set to the difference between the mixed phase functions
    ! in the convection and cloud scheme, but this is not calculated, so is zero and
    ! the functions must be the same
    ZDTFORC = ZDTDP*ZWTOT*PTSPHY + ZDTDIAB
    ZQOLD = ZQSMIX(JL, JK)
    ZTOLD = ZTP1(JL, JK)
    ZTP1(JL, JK) = ZTP1(JL, JK) + ZDTFORC
    ZTP1(JL, JK) = MAX(ZTP1(JL, JK), 160.0_JPRB)
    LLFLAG = .true.

    ! Formerly a call to CUADJTQ(..., ICALL=5)
    ZQP = 1.0_JPRB / PAP(JL, JK)
    ZQSAT = FOEEWM(ZTP1(JL, JK))*ZQP
    ZQSAT = MIN(0.5_JPRB, ZQSAT)
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    ZCOND = (ZQSMIX(JL, JK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(ZTP1(JL, JK)))
    ZTP1(JL, JK) = ZTP1(JL, JK) + FOELDCPM(ZTP1(JL, JK))*ZCOND
    ZQSMIX(JL, JK) = ZQSMIX(JL, JK) - ZCOND
    ZQSAT = FOEEWM(ZTP1(JL, JK))*ZQP
    ZQSAT = MIN(0.5_JPRB, ZQSAT)
    ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSAT)
    ZQSAT = ZQSAT*ZCOR
    ZCOND1 = (ZQSMIX(JL, JK) - ZQSAT) / (1.0_JPRB + ZQSAT*ZCOR*FOEDEM(ZTP1(JL, JK)))
    ZTP1(JL, JK) = ZTP1(JL, JK) + FOELDCPM(ZTP1(JL, JK))*ZCOND1
    ZQSMIX(JL, JK) = ZQSMIX(JL, JK) - ZCOND1

    ZDQS = ZQSMIX(JL, JK) - ZQOLD
    ZQSMIX(JL, JK) = ZQOLD
    ZTP1(JL, JK) = ZTOLD

    !----------------------------------------------------------------------
    ! 3.4a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
    ! ----------------------------------------------------------------------
    ! Erosion term is LINEAR in L
    ! Changed to be uniform distribution in cloud region


    ! Previous function based on DELTA DISTRIBUTION in cloud:
    IF (ZDQS > 0.0_JPRB) THEN
      !    If subsidence evaporation term is turned off, then need to use updated
      !    liquid and cloud here?
      !    ZLEVAP = MAX(ZA(JL,JK)+ZACUST(JL),1.0_JPRB)*MIN(ZDQS(JL),ZLICLD(JL)+ZLFINALSUM(JL))
      ZLEVAP = ZA(JL, JK)*MIN(ZDQS, ZLICLD(JK))
      ZLEVAP = MIN(ZLEVAP, ZEVAPLIMMIX(JK))
      ZLEVAP = MIN(ZLEVAP, MAX(ZQSMIX(JL, JK) - ZQX(JL, JK, NCLDQV), 0.0_JPRB))

      ! For first guess call
      ZLEVAPL = ZLIQFRAC(JL, JK)*ZLEVAP
      ZLEVAPI = ZICEFRAC(JL, JK)*ZLEVAP

      ZSOLQA(NCLDQV, NCLDQL, JK) = ZSOLQA(NCLDQV, NCLDQL, JK) + ZLIQFRAC(JL, JK)*ZLEVAP
      ZSOLQA(NCLDQL, NCLDQV, JK) = ZSOLQA(NCLDQL, NCLDQV, JK) - ZLIQFRAC(JL, JK)*ZLEVAP

      ZSOLQA(NCLDQV, NCLDQI, JK) = ZSOLQA(NCLDQV, NCLDQI, JK) + ZICEFRAC(JL, JK)*ZLEVAP
      ZSOLQA(NCLDQI, NCLDQV, JK) = ZSOLQA(NCLDQI, NCLDQV, JK) - ZICEFRAC(JL, JK)*ZLEVAP

    END IF


    !----------------------------------------------------------------------
    ! 3.4b ZDQS(JL) < 0: FORMATION OF CLOUDS
    !----------------------------------------------------------------------
    ! (1) Increase of cloud water in existing clouds
    IF (ZA(JL, JK) > ZEPSEC .and. ZDQS <= -YDECLDP%RLMIN) THEN

      ZLCOND1 = MAX(-ZDQS, 0.0_JPRB)        !new limiter

      !old limiter (significantly improves upper tropospheric humidity rms)
      IF (ZA(JL, JK) > 0.99_JPRB) THEN
        ZCOR = 1.0_JPRB / (1.0_JPRB - YDCST%RETV*ZQSMIX(JL, JK))
        ZCDMAX = (ZQX(JL, JK, NCLDQV) - ZQSMIX(JL, JK)) / (1.0_JPRB + ZCOR*ZQSMIX(JL, JK)*FOEDEM(ZTP1(JL, JK)))
      ELSE
        ZCDMAX = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSMIX(JL, JK)) / ZA(JL, JK)
      END IF
      ZLCOND1 = MAX(MIN(ZLCOND1, ZCDMAX), 0.0_JPRB)
      ! end old limiter

      ZLCOND1 = ZA(JL, JK)*ZLCOND1
      IF (ZLCOND1 < YDECLDP%RLMIN)       ZLCOND1 = 0.0_JPRB

      !-------------------------------------------------------------------------
      ! All increase goes into liquid unless so cold cloud homogeneously freezes
      ! Include new liquid formation in first guess value, otherwise liquid
      ! remains at cold temperatures until next timestep.
      !-------------------------------------------------------------------------
      IF (ZTP1(JL, JK) > YDECLDP%RTHOMO) THEN
        ZSOLQA(NCLDQL, NCLDQV, JK) = ZSOLQA(NCLDQL, NCLDQV, JK) + ZLCOND1
        ZSOLQA(NCLDQV, NCLDQL, JK) = ZSOLQA(NCLDQV, NCLDQL, JK) - ZLCOND1
        ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) + ZLCOND1
      ELSE
        ZSOLQA(NCLDQI, NCLDQV, JK) = ZSOLQA(NCLDQI, NCLDQV, JK) + ZLCOND1
        ZSOLQA(NCLDQV, NCLDQI, JK) = ZSOLQA(NCLDQV, NCLDQI, JK) - ZLCOND1
        ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + ZLCOND1
      END IF
    END IF

    ! (2) Generation of new clouds (da/dt>0)


    IF (ZDQS <= -YDECLDP%RLMIN .and. ZA(JL, JK) < 1.0_JPRB - ZEPSEC) THEN

      !---------------------------
      ! Critical relative humidity
      !---------------------------
      ZRHC = YDECLDP%RAMID
      ZSIGK = PAP(JL, JK) / PAPH(JL, KLEV + 1)
      ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
      IF (ZSIGK > 0.8_JPRB) THEN
        ZRHC = YDECLDP%RAMID + (1.0_JPRB - YDECLDP%RAMID)*((ZSIGK - 0.8_JPRB) / 0.2_JPRB)**2
      END IF

      ! Commented out for CY37R1 to reduce humidity in high trop and strat
      !      ! Increase RHcrit to 1.0 towards the tropopause (trop-0.2) and above
      !      ZBOTT=ZTRPAUS(JL)+0.2_JPRB
      !      IF(ZSIGK < ZBOTT) THEN
      !        ZRHC=RAMID+(1.0_JPRB-RAMID)*MIN(((ZBOTT-ZSIGK)/0.2_JPRB)**2,1.0_JPRB)
      !      ENDIF

      !---------------------------
      ! Supersaturation options
      !---------------------------
      IF (YDECLDP%NSSOPT == 0) THEN
        ! No scheme
        ZQE = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSICE(JL, JK)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZQE = MAX(0.0_JPRB, ZQE)
      ELSE IF (YDECLDP%NSSOPT == 1) THEN
        ! Tompkins
        ZQE = (ZQX(JL, JK, NCLDQV) - ZA(JL, JK)*ZQSICE(JL, JK)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ZQE = MAX(0.0_JPRB, ZQE)
      ELSE IF (YDECLDP%NSSOPT == 2) THEN
        ! Lohmann and Karcher
        ZQE = ZQX(JL, JK, NCLDQV)
      ELSE IF (YDECLDP%NSSOPT == 3) THEN
        ! Gierens
        ZQE = ZQX(JL, JK, NCLDQV) + ZLI(JL, JK)
      END IF

      IF (ZTP1(JL, JK) >= YDCST%RTT .or. YDECLDP%NSSOPT == 0) THEN
        ! No ice supersaturation allowed
        ZFAC = 1.0_JPRB
      ELSE
        ! Ice supersaturation
        ZFAC = ZFOKOOP(JK)
      END IF

      IF (ZQE >= ZRHC*ZQSICE(JL, JK)*ZFAC .and. ZQE < ZQSICE(JL, JK)*ZFAC) THEN
        ! note: not **2 on 1-a term if ZQE is used.
        ! Added correction term ZFAC to numerator 15/03/2010
        ZACOND = -((1.0_JPRB - ZA(JL, JK))*ZFAC*ZDQS) / MAX(2.0_JPRB*(ZFAC*ZQSICE(JL, JK) - ZQE), ZEPSEC)

        ZACOND = MIN(ZACOND, 1.0_JPRB - ZA(JL, JK))          !PUT THE LIMITER BACK

        ! Linear term:
        ! Added correction term ZFAC 15/03/2010
        ZLCOND2 = -ZFAC*ZDQS*0.5_JPRB*ZACOND          !mine linear

        ! new limiter formulation
        ZZDL = (2.0_JPRB*(ZFAC*ZQSICE(JL, JK) - ZQE)) / MAX(ZEPSEC, 1.0_JPRB - ZA(JL, JK))
        ! Added correction term ZFAC 15/03/2010
        IF (ZFAC*ZDQS < -ZZDL) THEN
          ! ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
          ZLCONDLIM = (ZA(JL, JK) - 1.0_JPRB)*ZFAC*ZDQS - ZFAC*ZQSICE(JL, JK) + ZQX(JL, JK, NCLDQV)
          ZLCOND2 = MIN(ZLCOND2, ZLCONDLIM)
        END IF
        ZLCOND2 = MAX(ZLCOND2, 0.0_JPRB)

        IF (ZLCOND2 < YDECLDP%RLMIN .or. (1.0_JPRB - ZA(JL, JK)) < ZEPSEC) THEN
          ZLCOND2 = 0.0_JPRB
          ZACOND = 0.0_JPRB
        END IF
        IF (ZLCOND2 == 0.0_JPRB)         ZACOND = 0.0_JPRB

        ! Large-scale generation is LINEAR in A and LINEAR in L
        ZSOLAC(JK) = ZSOLAC(JK) + ZACOND          !linear

        !------------------------------------------------------------------------
        ! All increase goes into liquid unless so cold cloud homogeneously freezes
        ! Include new liquid formation in first guess value, otherwise liquid
        ! remains at cold temperatures until next timestep.
        !------------------------------------------------------------------------
        IF (ZTP1(JL, JK) > YDECLDP%RTHOMO) THEN
          ZSOLQA(NCLDQL, NCLDQV, JK) = ZSOLQA(NCLDQL, NCLDQV, JK) + ZLCOND2
          ZSOLQA(NCLDQV, NCLDQL, JK) = ZSOLQA(NCLDQV, NCLDQL, JK) - ZLCOND2
          ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) + ZLCOND2
        ELSE
          ! homogeneous freezing
          ZSOLQA(NCLDQI, NCLDQV, JK) = ZSOLQA(NCLDQI, NCLDQV, JK) + ZLCOND2
          ZSOLQA(NCLDQV, NCLDQI, JK) = ZSOLQA(NCLDQV, NCLDQI, JK) - ZLCOND2
          ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + ZLCOND2
        END IF

      END IF
    END IF

    ! Loki end region-hoist group(section3p4)
  END SUBROUTINE section3p4
  SUBROUTINE section3p7 (KLON, KLEV, YDTHF, PTSPHY, JK, ZA, ZFOKOOP, ZRHO, YDECLDP, PAP, ZDP, YDCST, JL, IDEPICE, ZICECLD, ZTP1, ZQXFG,  &
  & ZCLDTOPDIST, ZSOLQA, ZCVDS, ZTERM1, ZCORRFAC, ZADD, ZICE0, ZFACX1I, ZVPICE, ZCORRFAC2, ZPR02, ZBDD, ZTCG, ZICENUCLEI,  &
  & ZTERM2, ZVPLIQ, ZAPLUSB, ZINEW, ZDEPOS, ZINFACTOR)
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMPHYDER, ONLY: state_type
    USE YOECLDP, ONLY: NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
    USE YOECLDP, ONLY: TECLDP
    USE YOEPHLI, ONLY: TEPHLI
    USE YOMCST, ONLY: TOMCST
    USE YOETHF, ONLY: TOETHF
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! Number of grid points
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV    ! Number of levels
    TYPE(TOETHF), INTENT(IN) :: YDTHF
    REAL(KIND=JPRB), INTENT(OUT) :: ZPR02
    TYPE(TECLDP), INTENT(IN) :: YDECLDP
    REAL(KIND=JPRB), INTENT(OUT) :: ZBDD
    REAL(KIND=JPRB), INTENT(IN) :: PAP(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZDP(KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
    REAL(KIND=JPRB), INTENT(OUT) :: ZADD
    INTEGER(KIND=JPIM), INTENT(IN) :: JK
    REAL(KIND=JPRB), INTENT(OUT) :: ZVPICE
    REAL(KIND=JPRB), INTENT(INOUT) :: ZCLDTOPDIST
    REAL(KIND=JPRB), INTENT(INOUT) :: ZQXFG(NCLV, KLEV)
    INTEGER(KIND=JPIM), INTENT(IN) :: JL
    TYPE(TOMCST), INTENT(IN) :: YDCST
    REAL(KIND=JPRB), INTENT(OUT) :: ZICENUCLEI
    REAL(KIND=JPRB), INTENT(OUT) :: ZTERM2
    REAL(KIND=JPRB), INTENT(OUT) :: ZDEPOS
    REAL(KIND=JPRB), INTENT(OUT) :: ZTERM1
    REAL(KIND=JPRB), INTENT(IN) :: ZA(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZICE0
    REAL(KIND=JPRB), INTENT(OUT) :: ZFACX1I
    INTEGER(KIND=JPIM), INTENT(IN) :: IDEPICE
    REAL(KIND=JPRB), INTENT(IN) :: ZFOKOOP(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZCVDS
    REAL(KIND=JPRB), INTENT(OUT) :: ZINEW
    REAL(KIND=JPRB), INTENT(OUT) :: ZAPLUSB
    REAL(KIND=JPRB), INTENT(IN) :: ZRHO(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZINFACTOR
    REAL(KIND=JPRB), INTENT(IN) :: ZICECLD(KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZCORRFAC
    REAL(KIND=JPRB), INTENT(OUT) :: ZCORRFAC2
    REAL(KIND=JPRB), INTENT(INOUT) :: ZSOLQA(NCLV, NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZTCG
    REAL(KIND=JPRB), INTENT(IN) :: ZTP1(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZVPLIQ
#include "fcttre.ycst.h"
!$acc routine seq
    ! Loki region-hoist group(section3p7)

    !----------------------------------------------------------------------
    ! 3.7 Growth of ice by vapour deposition
    !----------------------------------------------------------------------
    ! Following Rotstayn et al. 2001:
    ! does not use the ice nuclei number from cloudaer.F90
    ! but rather a simple Meyers et al. 1992 form based on the
    ! supersaturation and assuming clouds are saturated with
    ! respect to liquid water (well mixed), (or Koop adjustment)
    ! Growth considered as sink of liquid water if present so
    ! Bergeron-Findeisen adjustment in autoconversion term no longer needed
    !----------------------------------------------------------------------

    !--------------------------------------------------------
    !-
    !- Ice deposition following Rotstayn et al. (2001)
    !-  (monodisperse ice particle size distribution)
    !-
    !--------------------------------------------------------
    IF (IDEPICE == 1) THEN


      !--------------------------------------------------------------
      ! Calculate distance from cloud top
      ! defined by cloudy layer below a layer with cloud frac <0.01
      ! ZDZ = ZDP(JL)/(ZRHO(JL)*YDCST%RG)
      !--------------------------------------------------------------

      IF (ZA(JL, JK - 1) < YDECLDP%RCLDTOPCF .and. ZA(JL, JK) >= YDECLDP%RCLDTOPCF) THEN
        ZCLDTOPDIST = 0.0_JPRB
      ELSE
        ZCLDTOPDIST = ZCLDTOPDIST + ZDP(JK) / (ZRHO(JK)*YDCST%RG)
      END IF

      !--------------------------------------------------------------
      ! only treat depositional growth if liquid present. due to fact
      ! that can not model ice growth from vapour without additional
      ! in-cloud water vapour variable
      !--------------------------------------------------------------
      IF (ZTP1(JL, JK) < YDCST%RTT .and. ZQXFG(NCLDQL, JK) > YDECLDP%RLMIN) THEN
        ! T<273K

        ZVPICE = (FOEEICE(ZTP1(JL, JK))*YDCST%RV) / YDCST%RD
        ZVPLIQ = ZVPICE*ZFOKOOP(JK)
        ZICENUCLEI = 1000.0_JPRB*EXP((12.96_JPRB*(ZVPLIQ - ZVPICE)) / ZVPLIQ - 0.639_JPRB)

        !------------------------------------------------
        !   2.4e-2 is conductivity of air
        !   8.8 = 700**1/3 = density of ice to the third
        !------------------------------------------------
        ZADD = (YDCST%RLSTT*(YDCST%RLSTT / (YDCST%RV*ZTP1(JL, JK)) - 1.0_JPRB)) / (2.4E-2_JPRB*ZTP1(JL, JK))
        ZBDD = (YDCST%RV*ZTP1(JL, JK)*PAP(JL, JK)) / (2.21_JPRB*ZVPICE)
        ZCVDS = (7.8_JPRB*(ZICENUCLEI / ZRHO(JK))**0.666_JPRB*(ZVPLIQ - ZVPICE)) / (8.87_JPRB*(ZADD + ZBDD)*ZVPICE)

        !-----------------------------------------------------
        ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
        !-----------------------------------------------------
        ZICE0 = MAX(ZICECLD(JK), (ZICENUCLEI*YDECLDP%RICEINIT) / ZRHO(JK))

        !------------------
        ! new value of ice:
        !------------------
        ZINEW = (0.666_JPRB*ZCVDS*PTSPHY + ZICE0**0.666_JPRB)**1.5_JPRB

        !---------------------------
        ! grid-mean deposition rate:
        !---------------------------
        ZDEPOS = MAX(ZA(JL, JK)*(ZINEW - ZICE0), 0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit deposition to liquid water amount
        ! If liquid is all frozen, ice would use up reservoir of water
        ! vapour in excess of ice saturation mixing ratio - However this
        ! can not be represented without a in-cloud humidity variable. Using
        ! the grid-mean humidity would imply a large artificial horizontal
        ! flux from the clear sky to the cloudy area. We thus rely on the
        ! supersaturation check to clean up any remaining supersaturation
        !--------------------------------------------------------------------
        ZDEPOS = MIN(ZDEPOS, ZQXFG(NCLDQL, JK))          ! limit to liquid water amount

        !--------------------------------------------------------------------
        ! At top of cloud, reduce deposition rate near cloud top to account for
        ! small scale turbulent processes, limited ice nucleation and ice fallout
        !--------------------------------------------------------------------
        !      ZDEPOS = ZDEPOS*MIN(RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH,1.0_JPRB)
        ! Change to include dependence on ice nuclei concentration
        ! to increase deposition rate with decreasing temperatures
        ZINFACTOR = MIN(ZICENUCLEI / 15000._JPRB, 1.0_JPRB)
        ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB - ZINFACTOR)*(YDECLDP%RDEPLIQREFRATE + ZCLDTOPDIST / YDECLDP%RDEPLIQREFDEPTH),  &
        & 1.0_JPRB)

        !--------------
        ! add to matrix
        !--------------
        ZSOLQA(NCLDQI, NCLDQL, JK) = ZSOLQA(NCLDQI, NCLDQL, JK) + ZDEPOS
        ZSOLQA(NCLDQL, NCLDQI, JK) = ZSOLQA(NCLDQL, NCLDQI, JK) - ZDEPOS
        ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + ZDEPOS
        ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) - ZDEPOS

      END IF

      !--------------------------------------------------------
      !-
      !- Ice deposition assuming ice PSD
      !-
      !--------------------------------------------------------
    ELSE IF (IDEPICE == 2) THEN


      !--------------------------------------------------------------
      ! Calculate distance from cloud top
      ! defined by cloudy layer below a layer with cloud frac <0.01
      ! ZDZ = ZDP(JL)/(ZRHO(JL)*YDCST%RG)
      !--------------------------------------------------------------

      IF (ZA(JL, JK - 1) < YDECLDP%RCLDTOPCF .and. ZA(JL, JK) >= YDECLDP%RCLDTOPCF) THEN
        ZCLDTOPDIST = 0.0_JPRB
      ELSE
        ZCLDTOPDIST = ZCLDTOPDIST + ZDP(JK) / (ZRHO(JK)*YDCST%RG)
      END IF

      !--------------------------------------------------------------
      ! only treat depositional growth if liquid present. due to fact
      ! that can not model ice growth from vapour without additional
      ! in-cloud water vapour variable
      !--------------------------------------------------------------
      IF (ZTP1(JL, JK) < YDCST%RTT .and. ZQXFG(NCLDQL, JK) > YDECLDP%RLMIN) THEN
        ! T<273K

        ZVPICE = (FOEEICE(ZTP1(JL, JK))*YDCST%RV) / YDCST%RD
        ZVPLIQ = ZVPICE*ZFOKOOP(JK)
        ZICENUCLEI = 1000.0_JPRB*EXP((12.96_JPRB*(ZVPLIQ - ZVPICE)) / ZVPLIQ - 0.639_JPRB)

        !-----------------------------------------------------
        ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
        !-----------------------------------------------------
        ZICE0 = MAX(ZICECLD(JK), (ZICENUCLEI*YDECLDP%RICEINIT) / ZRHO(JK))

        ! Particle size distribution
        ZTCG = 1.0_JPRB
        ZFACX1I = 1.0_JPRB

        ZAPLUSB =  &
        & YDECLDP%RCL_APB1*ZVPICE - YDECLDP%RCL_APB2*ZVPICE*ZTP1(JL, JK) + PAP(JL, JK)*YDECLDP%RCL_APB3*ZTP1(JL, JK)**3._JPRB
        ZCORRFAC = (1.0_JPRB / ZRHO(JK))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL, JK) / 273.0_JPRB)**1.5_JPRB)*(393.0_JPRB / (ZTP1(JL, JK) + 120.0_JPRB))

        ZPR02 = (ZRHO(JK)*ZICE0*YDECLDP%RCL_CONST1I) / (ZTCG*ZFACX1I)

        ZTERM1 = ((ZVPLIQ - ZVPICE)*ZTP1(JL, JK)**2.0_JPRB*ZVPICE*ZCORRFAC2*ZTCG*YDECLDP%RCL_CONST2I*ZFACX1I) / (ZRHO(JK) &
        & *ZAPLUSB*ZVPICE)
        ZTERM2 = 0.65_JPRB*YDECLDP%RCL_CONST6I*ZPR02**YDECLDP%RCL_CONST4I + (YDECLDP%RCL_CONST3I*ZCORRFAC**0.5_JPRB*ZRHO(JK) &
        & **0.5_JPRB*ZPR02**YDECLDP%RCL_CONST5I) / ZCORRFAC2**0.5_JPRB

        ZDEPOS = MAX(ZA(JL, JK)*ZTERM1*ZTERM2*PTSPHY, 0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit deposition to liquid water amount
        ! If liquid is all frozen, ice would use up reservoir of water
        ! vapour in excess of ice saturation mixing ratio - However this
        ! can not be represented without a in-cloud humidity variable. Using
        ! the grid-mean humidity would imply a large artificial horizontal
        ! flux from the clear sky to the cloudy area. We thus rely on the
        ! supersaturation check to clean up any remaining supersaturation
        !--------------------------------------------------------------------
        ZDEPOS = MIN(ZDEPOS, ZQXFG(NCLDQL, JK))          ! limit to liquid water amount

        !--------------------------------------------------------------------
        ! At top of cloud, reduce deposition rate near cloud top to account for
        ! small scale turbulent processes, limited ice nucleation and ice fallout
        !--------------------------------------------------------------------
        ! Change to include dependence on ice nuclei concentration
        ! to increase deposition rate with decreasing temperatures
        ZINFACTOR = MIN(ZICENUCLEI / 15000._JPRB, 1.0_JPRB)
        ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB - ZINFACTOR)*(YDECLDP%RDEPLIQREFRATE + ZCLDTOPDIST / YDECLDP%RDEPLIQREFDEPTH),  &
        & 1.0_JPRB)

        !--------------
        ! add to matrix
        !--------------
        ZSOLQA(NCLDQI, NCLDQL, JK) = ZSOLQA(NCLDQI, NCLDQL, JK) + ZDEPOS
        ZSOLQA(NCLDQL, NCLDQI, JK) = ZSOLQA(NCLDQL, NCLDQI, JK) - ZDEPOS
        ZQXFG(NCLDQI, JK) = ZQXFG(NCLDQI, JK) + ZDEPOS
        ZQXFG(NCLDQL, JK) = ZQXFG(NCLDQL, JK) - ZDEPOS
      END IF

    END IF
    ! on IDEPICE

    ! Loki end region-hoist group(section3p7)
  END SUBROUTINE section3p7
  SUBROUTINE section6 (KLON, KLEV, YDTHF, ZQTMST, JK, ZQXN, ZPSUPSATSRCE, ZQX0, ZCONVSINK, ZDA, JL, IPHASE, ZQX, ZCONVSRCE, ZFALLSRCE,  &
  & ZFALLSINK, TENDENCY_LOC_T, TENDENCY_LOC_a, TENDENCY_LOC_q, TENDENCY_LOC_CLD, ZFLUXQ)
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMPHYDER, ONLY: state_type
    USE YOECLDP, ONLY: NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
    USE YOECLDP, ONLY: TECLDP
    USE YOEPHLI, ONLY: TEPHLI
    USE YOMCST, ONLY: TOMCST
    USE YOETHF, ONLY: TOETHF
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! Number of grid points
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV    ! Number of levels
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_a(KLON, KLEV)
    TYPE(TOETHF), INTENT(IN) :: YDTHF
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_T(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZQX(KLON, KLEV, NCLV)
    INTEGER(KIND=JPIM), INTENT(IN) :: IPHASE(NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: ZCONVSINK(NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: ZFLUXQ(NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: ZQTMST
    INTEGER(KIND=JPIM), INTENT(IN) :: JK
    REAL(KIND=JPRB), INTENT(IN) :: ZPSUPSATSRCE(NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZQXN(NCLV, KLEV)
    INTEGER(KIND=JPIM), INTENT(IN) :: JL
    REAL(KIND=JPRB), INTENT(IN) :: ZQX0(KLON, KLEV, NCLV)
    INTEGER(KIND=JPIM) :: JM
    REAL(KIND=JPRB), INTENT(IN) :: ZFALLSINK(NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_q(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_CLD(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: ZCONVSRCE(NCLV, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZDA(KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: ZFALLSRCE(NCLV, KLEV)
#include "fccld.ydthf.h"
!$acc routine seq
    !######################################################################
    !              6  *** UPDATE TENDANCIES ***
    !######################################################################

    !--------------------------------
    ! 6.1 Temperature and CLV budgets
    !--------------------------------

    DO JM=1,NCLV - 1

      ! calculate fluxes in and out of box for conservation of TL
      ZFLUXQ(JM) =  &
      & ZPSUPSATSRCE(JM, JK) + ZCONVSRCE(JM, JK) + ZFALLSRCE(JM, JK) - (ZFALLSINK(JM, JK) + ZCONVSINK(JM, JK))*ZQXN(JM, JK)

      IF (IPHASE(JM) == 1) THEN
        TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) + YDTHF%RALVDCP*(ZQXN(JM, JK) - ZQX(JL, JK, JM) - ZFLUXQ(JM))*ZQTMST
      END IF

      IF (IPHASE(JM) == 2) THEN
        TENDENCY_LOC_T(JL, JK) = TENDENCY_LOC_T(JL, JK) + YDTHF%RALSDCP*(ZQXN(JM, JK) - ZQX(JL, JK, JM) - ZFLUXQ(JM))*ZQTMST
      END IF

      !----------------------------------------------------------------------
      ! New prognostic tendencies - ice,liquid rain,snow
      ! Note: CLV arrays use PCLV in calculation of tendency while humidity
      !       uses ZQX. This is due to clipping at start of cloudsc which
      !       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
      !----------------------------------------------------------------------
      TENDENCY_LOC_CLD(JL, JK, JM) = TENDENCY_LOC_CLD(JL, JK, JM) + (ZQXN(JM, JK) - ZQX0(JL, JK, JM))*ZQTMST

    END DO

    !----------------------
    ! 6.2 Humidity budget
    !----------------------
    TENDENCY_LOC_q(JL, JK) = TENDENCY_LOC_Q(JL, JK) + (ZQXN(NCLDQV, JK) - ZQX(JL, JK, NCLDQV))*ZQTMST

    !-------------------
    ! 6.3 cloud cover
    !-----------------------
    TENDENCY_LOC_a(JL, JK) = TENDENCY_LOC_A(JL, JK) + ZDA(JK)*ZQTMST

    ! Loki region-hoist group(pcovptot) - region hoisted
    ! on vertical level JK
    !----------------------------------------------------------------------
    !                       END OF VERTICAL LOOP
    !----------------------------------------------------------------------

  END SUBROUTINE section6
  SUBROUTINE section8 (KLON, ZFOEALFA, ZRG_R, PAPH, PTSPHY, ZQTMST, PVFL, ZQX0, ZLNEG, PVFI, ZQXN2D, YDCST, JL, ZPFPLSX, PLUDE, KLEV,  &
  & PFCQSNG, PFCQLNG, PFPLSN, PFCQRNG, PFCQNNG, PFSQIF, PFSQLTUR, ZGDPH_R, PFSQSF, PFSQRF, PFHPSN, PFSQITUR, PFPLSL, PFSQLF,  &
  & PFHPSL, ZALFAW)
    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMPHYDER, ONLY: state_type
    USE YOECLDP, ONLY: NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
    USE YOECLDP, ONLY: TECLDP
    USE YOEPHLI, ONLY: TEPHLI
    USE YOMCST, ONLY: TOMCST
    USE YOETHF, ONLY: TOETHF
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! Number of grid points
    REAL(KIND=JPRB), INTENT(IN) :: ZRG_R
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: ZGDPH_R
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: ZALFAW
    REAL(KIND=JPRB), INTENT(IN) :: PVFL(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(IN) :: PTSPHY
    REAL(KIND=JPRB), INTENT(IN) :: ZQTMST
    INTEGER(KIND=JPIM) :: JK
    REAL(KIND=JPRB), INTENT(IN) :: PLUDE(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(KLON, KLEV + 1)
    INTEGER(KIND=JPIM), INTENT(IN) :: JL
    TYPE(TOMCST), INTENT(IN) :: YDCST
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: ZQXN2D(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: ZQX0(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: ZLNEG(KLON, KLEV, NCLV)
    REAL(KIND=JPRB), INTENT(IN) :: ZPFPLSX(KLON, KLEV + 1, NCLV)
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(KLON, KLEV + 1)
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: ZFOEALFA(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: PVFI(KLON, KLEV)
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(KLON, KLEV + 1)
    REAL(KIND=JPRB), INTENT(IN) :: PAPH(KLON, KLEV + 1)
!$acc routine seq

    !######################################################################
    !              8  *** FLUX/DIAGNOSTICS COMPUTATIONS ***
    !######################################################################

    !--------------------------------------------------------------------
    ! Copy general precip arrays back into PFP arrays for GRIB archiving
    ! Add rain and liquid fluxes, ice and snow fluxes
    !--------------------------------------------------------------------
!$acc loop seq
    DO JK=1,KLEV + 1
      PFPLSL(JL, JK) = ZPFPLSX(JL, JK, NCLDQR) + ZPFPLSX(JL, JK, NCLDQL)
      PFPLSN(JL, JK) = ZPFPLSX(JL, JK, NCLDQS) + ZPFPLSX(JL, JK, NCLDQI)
    END DO

    !--------
    ! Fluxes:
    !--------
    PFSQLF(JL, 1) = 0.0_JPRB
    PFSQIF(JL, 1) = 0.0_JPRB
    PFSQRF(JL, 1) = 0.0_JPRB
    PFSQSF(JL, 1) = 0.0_JPRB
    PFCQLNG(JL, 1) = 0.0_JPRB
    PFCQNNG(JL, 1) = 0.0_JPRB
    PFCQRNG(JL, 1) = 0.0_JPRB      !rain
    PFCQSNG(JL, 1) = 0.0_JPRB      !snow
    ! fluxes due to turbulence
    PFSQLTUR(JL, 1) = 0.0_JPRB
    PFSQITUR(JL, 1) = 0.0_JPRB

!$acc loop seq
    DO JK=1,KLEV

      ZGDPH_R = -ZRG_R*(PAPH(JL, JK + 1) - PAPH(JL, JK))*ZQTMST
      PFSQLF(JL, JK + 1) = PFSQLF(JL, JK)
      PFSQIF(JL, JK + 1) = PFSQIF(JL, JK)
      PFSQRF(JL, JK + 1) = PFSQLF(JL, JK)
      PFSQSF(JL, JK + 1) = PFSQIF(JL, JK)
      PFCQLNG(JL, JK + 1) = PFCQLNG(JL, JK)
      PFCQNNG(JL, JK + 1) = PFCQNNG(JL, JK)
      PFCQRNG(JL, JK + 1) = PFCQLNG(JL, JK)
      PFCQSNG(JL, JK + 1) = PFCQNNG(JL, JK)
      PFSQLTUR(JL, JK + 1) = PFSQLTUR(JL, JK)
      PFSQITUR(JL, JK + 1) = PFSQITUR(JL, JK)

      ZALFAW = ZFOEALFA(JL, JK)

      ! Liquid , LS scheme minus detrainment
      PFSQLF(JL, JK + 1) =  &
      & PFSQLF(JL, JK + 1) + (ZQXN2D(JL, JK, NCLDQL) - ZQX0(JL, JK, NCLDQL) + PVFL(JL, JK)*PTSPHY - ZALFAW*PLUDE(JL, JK))*ZGDPH_R
      ! liquid, negative numbers
      PFCQLNG(JL, JK + 1) = PFCQLNG(JL, JK + 1) + ZLNEG(JL, JK, NCLDQL)*ZGDPH_R

      ! liquid, vertical diffusion
      PFSQLTUR(JL, JK + 1) = PFSQLTUR(JL, JK + 1) + PVFL(JL, JK)*PTSPHY*ZGDPH_R

      ! Rain, LS scheme
      PFSQRF(JL, JK + 1) = PFSQRF(JL, JK + 1) + (ZQXN2D(JL, JK, NCLDQR) - ZQX0(JL, JK, NCLDQR))*ZGDPH_R
      ! rain, negative numbers
      PFCQRNG(JL, JK + 1) = PFCQRNG(JL, JK + 1) + ZLNEG(JL, JK, NCLDQR)*ZGDPH_R

      ! Ice , LS scheme minus detrainment
      PFSQIF(JL, JK + 1) = PFSQIF(JL, JK + 1) + (ZQXN2D(JL, JK, NCLDQI) - ZQX0(JL, JK, NCLDQI) + PVFI(JL, JK)*PTSPHY - (1.0_JPRB  &
      & - ZALFAW)*PLUDE(JL, JK))*ZGDPH_R
      ! ice, negative numbers
      PFCQNNG(JL, JK + 1) = PFCQNNG(JL, JK + 1) + ZLNEG(JL, JK, NCLDQI)*ZGDPH_R

      ! ice, vertical diffusion
      PFSQITUR(JL, JK + 1) = PFSQITUR(JL, JK + 1) + PVFI(JL, JK)*PTSPHY*ZGDPH_R

      ! snow, LS scheme
      PFSQSF(JL, JK + 1) = PFSQSF(JL, JK + 1) + (ZQXN2D(JL, JK, NCLDQS) - ZQX0(JL, JK, NCLDQS))*ZGDPH_R
      ! snow, negative numbers
      PFCQSNG(JL, JK + 1) = PFCQSNG(JL, JK + 1) + ZLNEG(JL, JK, NCLDQS)*ZGDPH_R
    END DO

    !-----------------------------------
    ! enthalpy flux due to precipitation
    !-----------------------------------
!$acc loop seq
    DO JK=1,KLEV + 1
      PFHPSL(JL, JK) = -YDCST%RLVTT*PFPLSL(JL, JK)
      PFHPSN(JL, JK) = -YDCST%RLSTT*PFPLSN(JL, JK)
    END DO

  END SUBROUTINE section8
END MODULE CLOUDSC_GPU_SCC_HOIST_MOD
