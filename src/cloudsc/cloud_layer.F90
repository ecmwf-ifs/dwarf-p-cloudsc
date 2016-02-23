SUBROUTINE CLOUD_LAYER( &
 ! Input quantities
 &  KDIM, LDSLPHY, LDMAINCALL, PAUX, state, tendency_cml, tendency_tmp, tendency_dyn, tendency_vdf, PRAD, &
 &  PSAVTENDSAT, PSURF, LLKEYS, &
 ! Input/Output quantities
 &  AUXL, FLUX, PDIAG, & 
 ! Output tendencies
 &  tendency_loc)

!**** *CLOUD_LAYER* - Layer routine calling cloud scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : key to activate SL physics
! LDMAINCALL:  key to distinguish between first simplified call or the full one
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state
! tendency_cml : D. V. for model tendencies  from processes before 
! tendency_tmp : D. V. for model tendencies (entering cloud) from processes before 
! tendency_dyn : D. V. for model tendencies from explicit dynamics
! tendency_vdf : D. V. for model tendencies from turbulence scheme
! PRAD     : D. V. for radiative quantities
! PSAVTENDSAT : Supersat clipped at previous time level in SLTEND
! PSURF    : D.V. for surface quantities
! LLKEYS   : D.V. for local switches

!     ==== Input/output ====
! AUXL         : Derived variable for local quantites
! FLUX         : Derived variable for fluxes
! PDIAG        : Derived variable for diagnostics

!    ==== Output tendencies from convection ====
! tendency_loc :  Output process tendencies


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-28  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      2015-01-10  R. Forbes  Precip type diags now here so called every timestp

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
!USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : RTT

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &  FLUX_TYPE, AUX_DIAG_LOCAL_TYPE, AUX_RAD_TYPE, KEYS_LOCAL_TYPE, &
  &  SURF_AND_MORE_TYPE, AUX_DIAG_TYPE
USE YOMPHY2  , ONLY : TSPHY
USE YOECLDP  , ONLY : YRECLDP, NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE SURFACE_FIELDS_MIX , ONLY : YRSURF
USE YOMMP0 , ONLY : MYPROC
USE YOMCT3 , ONLY : NSTEP

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
LOGICAL                        , INTENT (IN)   :: LDSLPHY
LOGICAL                        , INTENT (IN)   :: LDMAINCALL
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: state
TYPE (STATE_TYPE)              , INTENT (IN)   :: tendency_cml
TYPE (STATE_TYPE)              , INTENT (IN)   :: tendency_tmp
TYPE (STATE_TYPE)              , INTENT (IN)   :: tendency_dyn
TYPE (STATE_TYPE)              , INTENT (IN)   :: tendency_vdf
TYPE (AUX_RAD_TYPE)            , INTENT (IN)   :: PRAD
REAL(KIND=JPRB)                , INTENT (IN)   :: PSAVTENDSAT(KDIM%KLON,KDIM%KLEV)
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (KEYS_LOCAL_TYPE)         , INTENT (IN)   :: LLKEYS
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE)              , INTENT (OUT)  :: tendency_loc

!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB)    :: ZRAIN, ZSNOW, ZTOT, ZRAINFRAC_SFC

! Rain fraction at top of refreezing layer
REAL(KIND=JPRB)    :: ZRAINFRAC_TOPRFZ(KDIM%KLON) 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cloudsc.intfb.h"

!     ------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('CLOUD_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VD=>YRSURF%YSD_VD, YSD_VF=>YRSURF%YSD_VF)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CLOUDSC


if (LDMAINCALL .and. MYPROC==1 .and. NSTEP==3 .and. KDIM%KBL == 1) then

CALL CLOUDSC_OUT &
  & (KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV,&
  & TSPHY,&
  & state%T, state%q , tendency_cml, tendency_tmp,tendency_loc, &
  & tendency_vdf%a, tendency_vdf%cld(:,:,NCLDQL),tendency_vdf%cld(:,:,NCLDQI), tendency_dyn%a, &
  &     tendency_dyn%cld(:,:,NCLDQL), tendency_dyn%cld(:,:,NCLDQI),&
  & PRAD%PHRSW,    PRAD%PHRLW,&
  & PAUX%PVERVEL,  PAUX%PRSF1,    PAUX%PRS1,&
  & PSURF%PSD_VF(:,YSD_VF%YLSM%MP),     LLKEYS%LLCUM,    PDIAG%ITYPE,&
  & PDIAG%ZLU,     PDIAG%ZLUDE,    PDIAG%ZSNDE,   PDIAG%PMFU,     PDIAG%PMFD,&
  & LDSLPHY,  LDMAINCALL,&
  & state%a, &
  & state%cld, &
  & PSAVTENDSAT,&
!-- arrays for aerosol-cloud interactions
!!-- aerosol climatology     & ZQAER, 6, &
  & AUXL%ZLCRIT_AER,AUXL%ZICRIT_AER, &     
  & AUXL%ZRE_ICE, &     
  & AUXL%ZCCN,     AUXL%ZNICE,&
!--
  & PDIAG%PCOVPTOT,ZRAINFRAC_TOPRFZ,&
  & FLUX%PFCSQL,   FLUX%PFCSQN ,  FLUX%PFCQNNG,  FLUX%PFCQLNG,&
  & FLUX%PFSQRF,   FLUX%PFSQSF ,  FLUX%PFCQRNG,  FLUX%PFCQSNG,&
  & FLUX%PFSQLTUR, FLUX%PFSQITUR , &
  & FLUX%PFPLSL,   FLUX%PFPLSN,   FLUX%PFHPSL,   FLUX%PFHPSN,&
  & PSURF%PSD_XA, KDIM%KFLDX)  

endif

CALL CLOUDSC &
  & (KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV,&
  & TSPHY,&
  & state%T, state%q , tendency_cml, tendency_tmp,tendency_loc, &
  & tendency_vdf%a, tendency_vdf%cld(:,:,NCLDQL),tendency_vdf%cld(:,:,NCLDQI), tendency_dyn%a, &
  &     tendency_dyn%cld(:,:,NCLDQL), tendency_dyn%cld(:,:,NCLDQI),&
  & PRAD%PHRSW,    PRAD%PHRLW,&
  & PAUX%PVERVEL,  PAUX%PRSF1,    PAUX%PRS1,&
  & PSURF%PSD_VF(:,YSD_VF%YLSM%MP),     LLKEYS%LLCUM,    PDIAG%ITYPE,&
  & PDIAG%ZLU,     PDIAG%ZLUDE,    PDIAG%ZSNDE,   PDIAG%PMFU,     PDIAG%PMFD,&
  & LDSLPHY,  LDMAINCALL,&
  & state%a, &
  & state%cld, &
  & PSAVTENDSAT,&
!-- arrays for aerosol-cloud interactions
!!-- aerosol climatology     & ZQAER, 6, &
  & AUXL%ZLCRIT_AER,AUXL%ZICRIT_AER, &     
  & AUXL%ZRE_ICE, &     
  & AUXL%ZCCN,     AUXL%ZNICE,&
!--
  & PDIAG%PCOVPTOT,ZRAINFRAC_TOPRFZ,&
  & FLUX%PFCSQL,   FLUX%PFCSQN ,  FLUX%PFCQNNG,  FLUX%PFCQLNG,&
  & FLUX%PFSQRF,   FLUX%PFSQSF ,  FLUX%PFCQRNG,  FLUX%PFCQSNG,&
  & FLUX%PFSQLTUR, FLUX%PFSQITUR , &
  & FLUX%PFPLSL,   FLUX%PFPLSN,   FLUX%PFHPSL,   FLUX%PFHPSN,&
  & PSURF%PSD_XA, KDIM%KFLDX)  

!-------------------------------------------------------------------------------
!*         2.  Diagnose surface precipitation type and accumulate freezing rain
!-------------------------------------------------------------------------------
! WMO code table (4.201) for precipitation type
! 0 Reserved
! 1 Rain
! 2 Thunderstorm
! 3 Freezing Rain
! 4 Mixed/Ice
! 5 Snow
! 6 Wet snow
! 7 Melting snow (sleet)
! 8 Ice pellets

DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZRAIN = FLUX%PFPLCL(JL,KDIM%KLEV) + FLUX%PFPLSL(JL,KDIM%KLEV)
  !ZSNOW = FLUX%PFPLCN(JL,KDIM%KLEV) + FLUX%PFPLSN(JL,KDIM%KLEV)
  !ZTOT  = ZRAIN + ZSNOW
  ZTOT  = ZRAIN
  ! Default precipitation type this timestep
  PDIAG%PPRECTYPE(JL) = 0._JPRB
  ! Default freezing rain accumulation this timestep 
  PDIAG%PFZRA(JL) = 0._JPRB   

  IF (ZTOT > 1.E-10_JPRB) THEN   ! If there is precipitation
    ZRAINFRAC_SFC = ZRAIN/ZTOT

    ! If 2m T above freezing (dry bulb at the moment although should use wet bulb)
    IF (PSURF%PSD_VD(JL,YSD_VD%Y2T%MP) > RTT) THEN
    
      IF (ZRAINFRAC_SFC > 0.8_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 1._JPRB  ! Rain
      ELSEIF (ZRAINFRAC_SFC > 0.2_JPRB .AND. ZRAINFRAC_SFC  <=0.8_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 7._JPRB  ! Melting snow (sleet)
      ELSEIF (ZRAINFRAC_SFC <= 0.2_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 6._JPRB  ! Wet snow
      ELSEIF (ZRAINFRAC_SFC <= 0.01_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 5._JPRB  ! Dry snow
      ENDIF
    
    ! T below freezing
    ELSE

      IF (ZRAINFRAC_SFC <= 0.2_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 5._JPRB  ! Dry snow
      ELSEIF (ZRAINFRAC_SFC > 0.2 .AND. ZRAINFRAC_TOPRFZ(JL) > 0.8_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 3._JPRB  ! Freezing rain 
        PDIAG%PFZRA(JL)     = ZTOT     ! Freezing rain rate (assume = all precip)
      ELSEIF (ZRAINFRAC_SFC > 0.2 .AND. ZRAINFRAC_TOPRFZ(JL) <= 0.8_JPRB) THEN
        PDIAG%PPRECTYPE(JL) = 8._JPRB  ! Ice pellets
      ENDIF    

    ENDIF
  ENDIF
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('CLOUD_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUD_LAYER
