MODULE FCCLD_MOD
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
  USE PARKIND1,   ONLY : JPIM, JPRB
  USE YOETHF,     ONLY : RKOOP1, RKOOP2
  USE FCTTRE_MOD, ONLY : FOEELIQ, FOEEICE

  IMPLICIT NONE
  CONTAINS

  ELEMENTAL FUNCTION FOKOOP(PTARE)
    REAL(KIND=JPRB) :: FOKOOP
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOKOOP = MIN(RKOOP1-RKOOP2*PTARE,FOEELIQ(PTARE)/FOEEICE(PTARE))
  END FUNCTION FOKOOP

END MODULE FCCLD_MOD
