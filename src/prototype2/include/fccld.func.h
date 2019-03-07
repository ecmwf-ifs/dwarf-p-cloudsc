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
FOKOOP (PTARE) = MIN(RKOOP1-RKOOP2*PTARE,FOEELIQ(PTARE)/FOEEICE(PTARE))
