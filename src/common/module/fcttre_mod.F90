! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

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

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEDELTA(PTARE)
    REAL(KIND=JPRB) :: FOEDELTA
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEDELTA = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-RTT))
  END FUNCTION FOEDELTA

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEW(PTARE)
    REAL(KIND=JPRB) :: FOEEW
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEEW = R2ES*EXP (&
     &(R3LES*FOEDELTA(PTARE)+R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-RTT)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))
  END FUNCTION FOEEW

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEDE(PTARE)
    REAL(KIND=JPRB) :: FOEDE
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEDE = (FOEDELTA(PTARE)*R5ALVCP+(1.0_JPRB-FOEDELTA(PTARE))*R5ALSCP)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
  END FUNCTION FOEDE

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEDESU(PTARE)
    REAL(KIND=JPRB) :: FOEDESU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEDESU = (FOEDELTA(PTARE)*R5LES+(1.0_JPRB-FOEDELTA(PTARE))*R5IES)&
     &/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
  END FUNCTION FOEDESU

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELH(PTARE)
    REAL(KIND=JPRB) :: FOELH
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELH = FOEDELTA(PTARE)*RLVTT + (1.0_JPRB-FOEDELTA(PTARE))*RLSTT
  END FUNCTION FOELH

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELDCP(PTARE)
    REAL(KIND=JPRB) :: FOELDCP
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

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
ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEALFA(PTARE)
    REAL(KIND=JPRB) :: FOEALFA
    REAL(KIND=JPRB), VALUE  :: PTARE
    !$acc routine seq

    FOEALFA = MIN(1.0_JPRB,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
     &*RTWAT_RTICE_R)**2) 
  END FUNCTION FOEALFA

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE
ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEWM(PTARE)
    REAL(KIND=JPRB) :: FOEEWM
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEEWM = R2ES *&
     &(FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
     &(1.0_JPRB-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
  END FUNCTION FOEEWM

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOE_DEWM_DT(PTARE)
    REAL(KIND=JPRB) :: FOE_DEWM_DT
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOE_DEWM_DT = R2ES * ( &
     & R3LES*FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES)) &
     &    *(RTT-R4LES)/(PTARE-R4LES)**2 + &
     & R3IES*(1.0-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)) &
     &    *(RTT-R4IES)/(PTARE-R4IES)**2)
  END FUNCTION FOE_DEWM_DT

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEDEM(PTARE)
    REAL(KIND=JPRB) :: FOEDEM
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEDEM = FOEALFA(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
  END FUNCTION FOEDEM

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELDCPM(PTARE)
    REAL(KIND=JPRB) :: FOELDCPM
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELDCPM = FOEALFA(PTARE)*RALVDCP+(1.0_JPRB-FOEALFA(PTARE))*RALSDCP
  END FUNCTION FOELDCPM

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELHM(PTARE)
    REAL(KIND=JPRB) :: FOELHM
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELHM = FOEALFA(PTARE)*RLVTT+(1.0_JPRB-FOEALFA(PTARE))*RLSTT
  END FUNCTION FOELHM

  !     Temperature normalization for humidity background change of variable
  !        INPUT : PTARE = TEMPERATURE
ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOETB(PTARE)
    REAL(KIND=JPRB) :: FOETB
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

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
ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEALFCU(PTARE)
    REAL(KIND=JPRB) :: FOEALFCU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEALFCU = MIN(1.0_JPRB,((MAX(RTICECU,MIN(RTWAT,PTARE))&
     &-RTICECU)*RTWAT_RTICECU_R)**2) 
  END FUNCTION FOEALFCU

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE
ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEWMCU(PTARE)
    REAL(KIND=JPRB) :: FOEEWMCU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEEWMCU = R2ES *&
     &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
     &(1.0_JPRB-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
  END FUNCTION FOEEWMCU

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEDEMCU(PTARE)
    REAL(KIND=JPRB) :: FOEDEMCU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEDEMCU = FOEALFCU(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
     &(1.0_JPRB-FOEALFCU(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
  END FUNCTION FOEDEMCU

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELDCPMCU(PTARE)
    REAL(KIND=JPRB) :: FOELDCPMCU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELDCPMCU = FOEALFCU(PTARE)*RALVDCP+(1.0_JPRB-FOEALFCU(PTARE))*RALSDCP
  END FUNCTION FOELDCPMCU

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELHMCU(PTARE)
    REAL(KIND=JPRB) :: FOELHMCU
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

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

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEWMO(PTARE)
    REAL(KIND=JPRB) :: FOEEWMO
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEEWMO = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  END FUNCTION FOEEWMO

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEELIQ(PTARE)
    REAL(KIND=JPRB) :: FOEELIQ
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEELIQ = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  END FUNCTION FOEELIQ

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEICE(PTARE)
    REAL(KIND=JPRB) :: FOEEICE
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEEICE = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))
  END FUNCTION FOEEICE

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELSON(PTARE)
    REAL(KIND=JPRB) :: FOELSON
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELSON = EXP( -6096.9385_JPRB/PTARE + 21.2409642_JPRB &
                     - 2.711193E-2_JPRB * PTARE    &
                     + 1.673952E-5_JPRB * PTARE**2 &
                     + 2.433502_JPRB * LOG(PTARE))
  END FUNCTION FOELSON

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOELES_V(PTARE)
    REAL(KIND=JPRB) :: FOELES_V
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOELES_V=R3LES*(PTARE-RTT)/(PTARE-R4LES)
  END FUNCTION FOELES_V

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEIES_V(PTARE)
    REAL(KIND=JPRB) :: FOEIES_V
    REAL(KIND=JPRB), VALUE :: PTARE
    !$acc routine seq

    FOEIES_V=R3IES*(PTARE-RTT)/(PTARE-R4IES)
  END FUNCTION FOEIES_V

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEWM_V(PTARE,EXP1,EXP2)
    REAL(KIND=JPRB) :: FOEEWM_V
    REAL(KIND=JPRB), VALUE :: PTARE, EXP1, EXP2
    !$acc routine seq

    FOEEWM_V=R2ES*(FOEALFA(PTARE)*EXP1+(1.0_JPRB-FOEALFA(PTARE))*EXP2)
  END FUNCTION FOEEWM_V

ATTRIBUTES(DEVICE)  PURE ELEMENTAL FUNCTION FOEEWMCU_V(PTARE,EXP1,EXP2)
    REAL(KIND=JPRB) :: FOEEWMCU_V
    REAL(KIND=JPRB), VALUE :: PTARE, EXP1, EXP2
    !$acc routine seq

    FOEEWMCU_V = R2ES*(FOEALFCU(PTARE)*EXP1+(1.0_JPRB-FOEALFCU(PTARE))*EXP2)
  END FUNCTION FOEEWMCU_V

END MODULE
