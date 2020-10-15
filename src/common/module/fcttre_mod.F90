MODULE FCTTRE_MOD
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
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous 
!       nucleation. FOE_DEWM_DT provides an approximate first derivative
!       of FOEEWM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------

  USE PARKIND1, ONLY : JPIM, JPRB

  USE YOMCST, ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV, RA, RPI
  USE YOETHF, ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                   R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                   RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2

  IMPLICIT NONE
  CONTAINS

  !     *****************************************************************

  !                NO CONSIDERATION OF MIXED PHASES

  !     *****************************************************************

  ELEMENTAL FUNCTION FOEDELTA(PTARE)
    REAL(KIND=JPRB) :: FOEDELTA
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEDELTA = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-RTT))
  END FUNCTION FOEDELTA

  ELEMENTAL FUNCTION FOEEW(PTARE)
    REAL(KIND=JPRB) :: FOEEW
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEEW = R2ES*EXP (&
     &(R3LES*FOEDELTA(PTARE)+R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-RTT)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))
  END FUNCTION FOEEW

  ELEMENTAL FUNCTION FOEDE(PTARE)
    REAL(KIND=JPRB) :: FOEDE
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEDE = (FOEDELTA(PTARE)*R5ALVCP+(1.0_JPRB-FOEDELTA(PTARE))*R5ALSCP)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
  END FUNCTION FOEDE

  ELEMENTAL FUNCTION FOEDESU(PTARE)
    REAL(KIND=JPRB) :: FOEDESU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEDESU = (FOEDELTA(PTARE)*R5LES+(1.0_JPRB-FOEDELTA(PTARE))*R5IES)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
  END FUNCTION FOEDESU

  ELEMENTAL FUNCTION FOELH(PTARE)
    REAL(KIND=JPRB) :: FOELH
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELH = FOEDELTA(PTARE)*RLVTT + (1.0_JPRB-FOEDELTA(PTARE))*RLSTT
  END FUNCTION FOELH

  ELEMENTAL FUNCTION FOELDCP(PTARE)
    REAL(KIND=JPRB) :: FOELDCP
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELDCP = FOEDELTA(PTARE)*RALVDCP + (1.0_JPRB-FOEDELTA(PTARE))*RALSDCP
  END FUNCTION FOELDCP

  !     *****************************************************************

  !           CONSIDERATION OF MIXED PHASES

  !     *****************************************************************

  !     FOEALFA is calculated to distinguish the three cases:

  !                       FOEALFA=1            water phase
  !                       FOEALFA=0            ice phase
  !                       0 < FOEALFA < 1      mixed phase

  !               INPUT : PTARE = TEMPERATURE
  ELEMENTAL FUNCTION FOEALFA(PTARE)
    REAL(KIND=JPRB) :: FOEALFA
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEALFA = MIN(1.0_JPRB,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
     &*RTWAT_RTICE_R)**2) 
  END FUNCTION FOEALFA

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE
  ELEMENTAL FUNCTION FOEEWM(PTARE)
    REAL(KIND=JPRB) :: FOEEWM
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEEWM = R2ES *&
     &(FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
     &(1.0_JPRB-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
  END FUNCTION FOEEWM

  ELEMENTAL FUNCTION FOE_DEWM_DT(PTARE)
    REAL(KIND=JPRB) :: FOE_DEWM_DT
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOE_DEWM_DT = R2ES * ( &
     & R3LES*FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES)) &
     &    *(RTT-R4LES)/(PTARE-R4LES)**2 + &
     & R3IES*(1.0-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)) &
     &    *(RTT-R4IES)/(PTARE-R4IES)**2)
  END FUNCTION FOE_DEWM_DT

  ELEMENTAL FUNCTION FOEDEM(PTARE)
    REAL(KIND=JPRB) :: FOEDEM
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEDEM = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
  END FUNCTION FOEDEM

  ELEMENTAL FUNCTION FOELDCPM(PTARE)
    REAL(KIND=JPRB) :: FOELDCPM
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELDCPM = FOEALFA(PTARE)*RALVDCP+(1.0_JPRB-FOEALFA(PTARE))*RALSDCP
  END FUNCTION FOELDCPM

  ELEMENTAL FUNCTION FOELHM(PTARE)
    REAL(KIND=JPRB) :: FOELHM
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELHM = FOEALFA(PTARE)*RLVTT+(1.0_JPRB-FOEALFA(PTARE))*RLSTT
  END FUNCTION FOELHM

  !     Temperature normalization for humidity background change of variable
  !        INPUT : PTARE = TEMPERATURE
  ELEMENTAL FUNCTION FOETB(PTARE)
    REAL(KIND=JPRB) :: FOETB
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOETB = FOEALFA(PTARE)*R3LES*(RTT-R4LES)*(1.0_JPRB/(PTARE-R4LES)**2)+&
     &(1.0_JPRB-FOEALFA(PTARE))*R3IES*(RTT-R4IES)*(1.0_JPRB/(PTARE-R4IES)**2)
  END FUNCTION FOETB

  !     ------------------------------------------------------------------
  !     *****************************************************************

  !           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

  !     *****************************************************************

  !     FOEALFCU is calculated to distinguish the three cases:

  !                       FOEALFCU=1            water phase
  !                       FOEALFCU=0            ice phase
  !                       0 < FOEALFCU < 1      mixed phase

  !               INPUT : PTARE = TEMPERATURE
  ELEMENTAL FUNCTION FOEALFCU(PTARE)
    REAL(KIND=JPRB) :: FOEALFCU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEALFCU = MIN(1.0_JPRB,((MAX(RTICECU,MIN(RTWAT,PTARE))&
     &-RTICECU)*RTWAT_RTICECU_R)**2) 
  END FUNCTION FOEALFCU

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE
  ELEMENTAL FUNCTION FOEEWMCU(PTARE)
    REAL(KIND=JPRB) :: FOEEWMCU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEEWMCU = R2ES *&
     &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
     &(1.0_JPRB-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
  END FUNCTION FOEEWMCU

  ELEMENTAL FUNCTION FOEDEMCU(PTARE)
    REAL(KIND=JPRB) :: FOEDEMCU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEDEMCU = FOEALFCU(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
     &(1.0_JPRB-FOEALFCU(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
  END FUNCTION FOEDEMCU

  ELEMENTAL FUNCTION FOELDCPMCU(PTARE)
    REAL(KIND=JPRB) :: FOELDCPMCU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELDCPMCU = FOEALFCU(PTARE)*RALVDCP+(1.0_JPRB-FOEALFCU(PTARE))*RALSDCP
  END FUNCTION FOELDCPMCU

  ELEMENTAL FUNCTION FOELHMCU(PTARE)
    REAL(KIND=JPRB) :: FOELHMCU
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELHMCU = FOEALFCU(PTARE)*RLVTT+(1.0_JPRB-FOEALFCU(PTARE))*RLSTT
  END FUNCTION FOELHMCU

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
  !      unlike the FOEE functions does not include 1/(RETV+1.0_JPRB) factor

  ELEMENTAL FUNCTION FOEEWMO(PTARE)
    REAL(KIND=JPRB) :: FOEEWMO
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEEWMO = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  END FUNCTION FOEEWMO

  ELEMENTAL FUNCTION FOEELIQ(PTARE)
    REAL(KIND=JPRB) :: FOEELIQ
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEELIQ = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  END FUNCTION FOEELIQ

  ELEMENTAL FUNCTION FOEEICE(PTARE)
    REAL(KIND=JPRB) :: FOEEICE
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEEICE = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))
  END FUNCTION FOEEICE

  ELEMENTAL FUNCTION FOELSON(PTARE)
    REAL(KIND=JPRB) :: FOELSON
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELSON = EXP( -6096.9385_JPRB/PTARE + 21.2409642_JPRB &
	             - 2.711193E-2_JPRB * PTARE    &
                     + 1.673952E-5_JPRB * PTARE**2 &
		     + 2.433502_JPRB * LOG(PTARE))
  END FUNCTION FOELSON

  ELEMENTAL FUNCTION FOELES_V(PTARE)
    REAL(KIND=JPRB) :: FOELES_V
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOELES_V=R3LES*(PTARE-RTT)/(PTARE-R4LES)
  END FUNCTION FOELES_V

  ELEMENTAL FUNCTION FOEIES_V(PTARE)
    REAL(KIND=JPRB) :: FOEIES_V
    REAL(KIND=JPRB), INTENT(IN) :: PTARE
    FOEIES_V=R3IES*(PTARE-RTT)/(PTARE-R4IES)
  END FUNCTION FOEIES_V

  ELEMENTAL FUNCTION FOEEWM_V(PTARE,EXP1,EXP2)
    REAL(KIND=JPRB) :: FOEEWM_V
    REAL(KIND=JPRB), INTENT(IN) :: PTARE, EXP1, EXP2
    FOEEWM_V=R2ES*(FOEALFA(PTARE)*EXP1+(1.0_JPRB-FOEALFA(PTARE))*EXP2)
  END FUNCTION FOEEWM_V

  ELEMENTAL FUNCTION FOEEWMCU_V(PTARE,EXP1,EXP2)
    REAL(KIND=JPRB) :: FOEEWMCU_V
    REAL(KIND=JPRB), INTENT(IN) :: PTARE, EXP1, EXP2
    FOEEWMCU_V = R2ES*(FOEALFCU(PTARE)*EXP1+(1.0_JPRB-FOEALFCU(PTARE))*EXP2)
  END FUNCTION FOEEWMCU_V

  FUNCTION FOETB_FUNC( PTARE )
    REAL(KIND=JPRB) :: FOETB_FUNC
    REAL(KIND=JPRB) :: PTARE

    FOETB_FUNC=FOEALFA_FUNC(PTARE)*R3LES*(RTT-R4LES)*(1.0_JPRB/(PTARE-R4LES)**2)+&
             &(1.0_JPRB-FOEALFA_FUNC(PTARE))*R3IES*(RTT-R4IES)*(1.0_JPRB/(PTARE-R4IES)**2)
    RETURN
  END FUNCTION

  FUNCTION FOEALFA_FUNC(PTARE)
    REAL(KIND=JPRB) :: FOEALFA_FUNC
    REAL(KIND=JPRB) :: PTARE

    FOEALFA_FUNC = MIN(1.0_JPRB,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
                 &*RTWAT_RTICE_R)**2)

    RETURN
  END FUNCTION

END MODULE
