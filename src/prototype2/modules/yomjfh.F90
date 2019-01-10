MODULE YOMJFH

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

! Use of MASS library
!  N_VMASS: < or = 0 if not using MASS library vector routines
!           > 0      if using MASS library vector routines

INTEGER(KIND=JPIM) :: N_VMASS=0

!    -----------------------------------------------------------------

END MODULE YOMJFH
