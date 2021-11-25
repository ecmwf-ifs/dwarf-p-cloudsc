! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTERFACE
SUBROUTINE CLOUDSC_IN &
 & (iu,    KLON,    KLEV,&
 & PTSPHY,&
 & PT, PQ, tendency_cml,tendency_tmp, &
 & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
 & PHRSW,    PHRLW,&
 & PVERVEL,  PAP,      PAPH,&
 & PLSM,     LDCUM,    KTYPE, &
 & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD,&
 & LDSLPHY,  LDMAINCALL, &
 !---prognostic fields
 & PA,&
 & PCLV,  &
 & PSUPSAT, &
!-- arrays for aerosol-cloud interactions
!!! & PQAER,    KAER, &
 & PLCRIT_AER,PICRIT_AER,&
 & PRE_ICE,&
 & PCCN,     PNICE,&
 & PEXTRA,   KFLDX)  
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMPHYDER ,ONLY : STATE_TYPE
USE YOECLDP  , ONLY : NCLV
!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------

INTEGER(KIND=JPIM), INTENT(in) :: iu
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels

REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLCRIT_AER(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PICRIT_AER(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PRE_ICE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)

REAL(KIND=JPRB)   ,INTENT(OUT)    :: PTSPHY            ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PT(KLON,KLEV)    ! T at start of callpar
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PQ(KLON,KLEV)    ! Q at start of callpar
TYPE (STATE_TYPE) , INTENT (OUT)  :: tendency_cml   ! cumulative tendency used for final output
TYPE (STATE_TYPE) , INTENT (OUT)  :: tendency_tmp   ! cumulative tendency used as input
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNL(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNI(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVERVEL(KLON,KLEV) !Vertical velocity
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLSM(KLON)       ! Land fraction (0-1) 
LOGICAL           ,INTENT(OUT)    :: LDCUM(KLON)      ! Convection active
INTEGER(KIND=JPIM),INTENT(OUT)    :: KTYPE(KLON)      ! Convection type 0,1,2
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLU(KLON,KLEV)   ! Conv. condensate
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLUDE(KLON,KLEV) ! Conv. detrained water 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSNDE(KLON,KLEV) ! Conv. detrained snow
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
LOGICAL           ,INTENT(OUT)    :: LDSLPHY 
LOGICAL           ,INTENT(OUT)    :: LDMAINCALL       ! T if main call to cloudsc
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields

REAL(KIND=JPRB)   ,INTENT(OUT)    :: PCLV(KLON,KLEV,NCLV) 

 ! Supersat clipped at previous time level in SLTEND
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSUPSAT(KLON,KLEV)

END SUBROUTINE CLOUDSC_IN
END INTERFACE
