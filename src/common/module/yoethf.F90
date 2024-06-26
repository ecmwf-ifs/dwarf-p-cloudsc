! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOETHF

USE PARKIND1,    ONLY : JPIM, JPRB
USE FILE_IO_MOD, ONLY : LOAD_SCALAR

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: R2ES
REAL(KIND=JPRB) :: R3LES
REAL(KIND=JPRB) :: R3IES
REAL(KIND=JPRB) :: R4LES
REAL(KIND=JPRB) :: R4IES
REAL(KIND=JPRB) :: R5LES
REAL(KIND=JPRB) :: R5IES
REAL(KIND=JPRB) :: RVTMP2
REAL(KIND=JPRB) :: RHOH2O
REAL(KIND=JPRB) :: R5ALVCP
REAL(KIND=JPRB) :: R5ALSCP
REAL(KIND=JPRB) :: RALVDCP
REAL(KIND=JPRB) :: RALSDCP
REAL(KIND=JPRB) :: RALFDCP
REAL(KIND=JPRB) :: RTWAT
REAL(KIND=JPRB) :: RTBER
REAL(KIND=JPRB) :: RTBERCU
REAL(KIND=JPRB) :: RTICE
REAL(KIND=JPRB) :: RTICECU
REAL(KIND=JPRB) :: RTWAT_RTICE_R
REAL(KIND=JPRB) :: RTWAT_RTICECU_R
REAL(KIND=JPRB) :: RKOOP1
REAL(KIND=JPRB) :: RKOOP2

TYPE :: TOETHF
REAL(KIND=JPRB) :: R2ES
REAL(KIND=JPRB) :: R3LES
REAL(KIND=JPRB) :: R3IES
REAL(KIND=JPRB) :: R4LES
REAL(KIND=JPRB) :: R4IES
REAL(KIND=JPRB) :: R5LES
REAL(KIND=JPRB) :: R5IES
REAL(KIND=JPRB) :: RVTMP2
REAL(KIND=JPRB) :: RHOH2O
REAL(KIND=JPRB) :: R5ALVCP
REAL(KIND=JPRB) :: R5ALSCP
REAL(KIND=JPRB) :: RALVDCP
REAL(KIND=JPRB) :: RALSDCP
REAL(KIND=JPRB) :: RALFDCP
REAL(KIND=JPRB) :: RTWAT
REAL(KIND=JPRB) :: RTBER
REAL(KIND=JPRB) :: RTBERCU
REAL(KIND=JPRB) :: RTICE
REAL(KIND=JPRB) :: RTICECU
REAL(KIND=JPRB) :: RTWAT_RTICE_R
REAL(KIND=JPRB) :: RTWAT_RTICECU_R
REAL(KIND=JPRB) :: RKOOP1
REAL(KIND=JPRB) :: RKOOP2
END TYPE TOETHF

TYPE(TOETHF), ALLOCATABLE :: YRTHF

!     J.-J. MORCRETTE                   91/07/14  ADAPTED TO I.F.S.

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *R__ES*   REAL      *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                         MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                         ICE(*R_IES*).
!     *RVTMP2*  REAL      *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  REAL      *DENSITY OF LIQUID WATER.   (RATM/100.)
!     *R5ALVCP* REAL      *R5LES*RLVTT/RCPD
!     *R5ALSCP* REAL      *R5IES*RLSTT/RCPD
!     *RALVDCP* REAL      *RLVTT/RCPD
!     *RALSDCP* REAL      *RLSTT/RCPD
!     *RALFDCP* REAL      *RLMLT/RCPD
!     *RTWAT*   REAL      *RTWAT=RTT
!     *RTBER*   REAL      *RTBER=RTT-0.05
!     *RTBERCU  REAL      *RTBERCU=RTT-5.0
!     *RTICE*   REAL      *RTICE=RTT-0.1
!     *RTICECU* REAL      *RTICECU=RTT-23.0
!     *RKOOP?   REAL      *CONSTANTS TO DESCRIBE KOOP FORM FOR NUCLEATION
!     *RTWAT_RTICE_R*   REAL      *RTWAT_RTICE_R=1./(RTWAT-RTICE)
!     *RTWAT_RTICECU_R* REAL      *RTWAT_RTICECU_R=1./(RTWAT-RTICECU)

!$acc declare copyin(r2es, r3les, r3ies, r4les, r4ies, r5les, r5ies, &
!$acc   r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
!$acc   rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)

!$omp declare target(r2es, r3les, r3ies, r4les, r4ies, r5les, r5ies)
!$omp declare target(  r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu)
!$omp declare target(  rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)

!       ----------------------------------------------------------------

CONTAINS

  SUBROUTINE YOETHF_LOAD_PARAMETERS()
    CALL LOAD_SCALAR('R2ES', R2ES)
    CALL LOAD_SCALAR('R3LES', R3LES)
    CALL LOAD_SCALAR('R3IES', R3IES)
    CALL LOAD_SCALAR('R4LES', R4LES)
    CALL LOAD_SCALAR('R4IES', R4IES)
    CALL LOAD_SCALAR('R5LES', R5LES)
    CALL LOAD_SCALAR('R5IES', R5IES)
    CALL LOAD_SCALAR('R5ALVCP', R5ALVCP)
    CALL LOAD_SCALAR('R5ALSCP', R5ALSCP)
    CALL LOAD_SCALAR('RALVDCP', RALVDCP)
    CALL LOAD_SCALAR('RALSDCP', RALSDCP)
    CALL LOAD_SCALAR('RALFDCP', RALFDCP)
    CALL LOAD_SCALAR('RTWAT', RTWAT)
    CALL LOAD_SCALAR('RTICE', RTICE)
    CALL LOAD_SCALAR('RTICECU', RTICECU)
    CALL LOAD_SCALAR('RTWAT_RTICE_R', RTWAT_RTICE_R)
    CALL LOAD_SCALAR('RTWAT_RTICECU_R', RTWAT_RTICECU_R)
    CALL LOAD_SCALAR('RKOOP1', RKOOP1)
    CALL LOAD_SCALAR('RKOOP2', RKOOP2)
    CALL YRTHF_COPY_PARAMETERS()
!$acc update device(r2es, r3les, r3ies, r4les, r4ies, r5les, r5ies, &
!$acc   r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
!$acc   rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)
  END SUBROUTINE YOETHF_LOAD_PARAMETERS

  SUBROUTINE YRTHF_COPY_PARAMETERS()
    IF(.NOT.ALLOCATED(YRTHF)) ALLOCATE(YRTHF)
    YRTHF%R2ES            = R2ES
    YRTHF%R3LES           = R3LES
    YRTHF%R3IES           = R3IES
    YRTHF%R4LES           = R4LES
    YRTHF%R4IES           = R4IES
    YRTHF%R5LES           = R5LES
    YRTHF%R5IES           = R5IES
    YRTHF%R5ALVCP         = R5ALVCP
    YRTHF%R5ALSCP         = R5ALSCP
    YRTHF%RALVDCP         = RALVDCP
    YRTHF%RALSDCP         = RALSDCP
    YRTHF%RALFDCP         = RALFDCP
    YRTHF%RTWAT           = RTWAT
    YRTHF%RTICE           = RTICE
    YRTHF%RTICECU         = RTICECU
    YRTHF%RTWAT_RTICE_R   = RTWAT_RTICE_R
    YRTHF%RTWAT_RTICECU_R = RTWAT_RTICECU_R
    YRTHF%RKOOP1          = RKOOP1
    YRTHF%RKOOP2          = RKOOP2
  END SUBROUTINE YRTHF_COPY_PARAMETERS

END MODULE YOETHF
