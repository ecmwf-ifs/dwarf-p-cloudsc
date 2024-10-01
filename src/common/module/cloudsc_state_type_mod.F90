MODULE CLOUDSC_STATE_TYPE_MOD

! USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

type state_type
  !$loki dimension(klon,klev)
  REAL(KIND=JPRB), dimension(:,:), pointer :: u,v,T   ! GMV fields
  !$loki dimension(klon,klev)
  REAL(KIND=JPRB), dimension(:,:), pointer :: o3,q,a  ! GFL fields
  !$loki dimension(klon,klev,5)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: cld   ! composed cloud array
  !REAL(KIND=JPRB), dimension(:,:), pointer :: qsat    ! spec. humidity at saturation

#ifdef  USE_FIELD_API
  CLASS(FIELD_3RB), POINTER :: F_T, F_A, F_Q
  CLASS(FIELD_4RB), POINTER :: F_CLD
#endif
end type state_type

