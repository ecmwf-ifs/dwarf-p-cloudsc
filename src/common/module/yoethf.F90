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
USE FILE_IO_MOD, ONLY : LOAD_SCALAR, load_scalar_real

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

REAL(KIND=JPRB) :: R2ES_H, R3LES_H, R3IES_H, R4LES_H, R4IES_H, R5LES_H, R5IES_H, RVTMP2_H, RHOH2O_H, R5ALVCP_H, R5ALSCP_H, &
        & RALVDCP_H, RALSDCP_H, RALFDCP_H, RTWAT_H, RTBER_H, RTBERCU_H, RTICE_H, RTICECU_H, RTWAT_RTICE_R_H, &
        & RTWAT_RTICECU_R_H, RKOOP1_H, RKOOP2_H

ATTRIBUTES(CONSTANT) :: R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, RVTMP2, RHOH2O, R5ALVCP, R5ALSCP, &
        & RALVDCP, RALSDCP, RALFDCP, RTWAT, RTBER, RTBERCU, RTICE, RTICECU, RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
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

!$omp declare target(r2es, r3les, r3ies, r4les, r4ies, r5les, r5ies, &
!$omp   r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
!$omp   rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)

!       ----------------------------------------------------------------

CONTAINS

  SUBROUTINE YOETHF_LOAD_PARAMETERS()
    CALL LOAD_SCALAR_real('R2ES', R2ES_h)
    CALL LOAD_SCALAR_real('R3LES', R3LES_h)
    CALL LOAD_SCALAR_real('R3IES', R3IES_h)
    CALL LOAD_SCALAR_real('R4LES', R4LES_h)
    CALL LOAD_SCALAR_real('R4IES', R4IES_h)
    CALL LOAD_SCALAR_real('R5LES', R5LES_h)
    CALL LOAD_SCALAR_real('R5IES', R5IES_h)
    CALL LOAD_SCALAR_real('R5ALVCP', R5ALVCP_h)
    CALL LOAD_SCALAR_real('R5ALSCP', R5ALSCP_h)
    CALL LOAD_SCALAR_real('RALVDCP', RALVDCP_h)
    CALL LOAD_SCALAR_real('RALSDCP', RALSDCP_h)
    CALL LOAD_SCALAR_real('RALFDCP', RALFDCP_h)
    CALL LOAD_SCALAR_real('RTWAT', RTWAT_h)
    CALL LOAD_SCALAR_real('RTICE', RTICE_h)
    CALL LOAD_SCALAR_real('RTICECU', RTICECU_h)
    CALL LOAD_SCALAR_real('RTWAT_RTICE_R', RTWAT_RTICE_R_h)
    CALL LOAD_SCALAR_real('RTWAT_RTICECU_R', RTWAT_RTICECU_R_h)
    CALL LOAD_SCALAR_real('RKOOP1', RKOOP1_h)
    CALL LOAD_SCALAR_real('RKOOP2', RKOOP2_h)

    R2ES=R2ES_H; R3LES=R3LES_H; R3IES=R3IES_H; R4LES=R4LES_H; R4IES=R4IES_H; R5LES=R5LES_H; R5IES=R5IES_H; 
    RVTMP2=RVTMP2_H; RHOH2O=RHOH2O_H; R5ALVCP=R5ALVCP_H; R5ALSCP=R5ALSCP_H; RALVDCP=RALVDCP_H; RALSDCP=RALSDCP_H;
    RALFDCP=RALFDCP_H; RTWAT=RTWAT_H; RTBER=RTBER_H; RTBERCU=RTBERCU_H; RTICE=RTICE_H; RTICECU=RTICECU_H;
    RTWAT_RTICE_R=RTWAT_RTICE_R_H; RTWAT_RTICECU_R=RTWAT_RTICECU_R_H; RKOOP1=RKOOP1_H; RKOOP2=RKOOP2_H

  END SUBROUTINE YOETHF_LOAD_PARAMETERS

END MODULE YOETHF
