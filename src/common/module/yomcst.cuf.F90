! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCST_CUF

USE PARKIND1,    ONLY : JPRB
USE YOMCST, ONLY: RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV

IMPLICIT NONE

REAL(KIND=JPRB), CONSTANT :: RG_D
REAL(KIND=JPRB), CONSTANT :: RD_D
REAL(KIND=JPRB), CONSTANT :: RV_D
REAL(KIND=JPRB), CONSTANT :: RCPD_D
REAL(KIND=JPRB), CONSTANT :: RETV_D
REAL(KIND=JPRB), CONSTANT :: RTT_D
REAL(KIND=JPRB), CONSTANT :: RLVTT_D
REAL(KIND=JPRB), CONSTANT :: RLSTT_D
REAL(KIND=JPRB), CONSTANT :: RLMLT_D


!    ------------------------------------------------------------------

CONTAINS

ATTRIBUTES(HOST)  SUBROUTINE YOMCST_UPDATE_DEVICE()
  RG_D=RG; RD_D=RD; RCPD_D=RCPD; RETV_D=RETV; RLVTT_D=RLVTT;
  RLSTT_D=RLSTT; RLMLT_D=RLMLT; RTT_D=RTT; RV_D=RV
END SUBROUTINE YOMCST_UPDATE_DEVICE


END MODULE YOMCST_CUF
