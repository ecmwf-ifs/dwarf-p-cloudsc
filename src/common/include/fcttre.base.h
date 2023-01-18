!*
!     ------------------------------------------------------------------

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
!       also a mix of both for the temperature range  _PREFIX2_ RTICE < T <  _PREFIX2_ RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous 
!       nucleation. FOE_DEWM_DT provides an approximate first derivative
!       of FOEEWM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------
!     *****************************************************************

!                NO CONSIDERATION OF MIXED PHASES

!     *****************************************************************
REAL(KIND=JPRB) :: FOEDELTA
REAL(KIND=JPRB) :: PTARE
FOEDELTA (PTARE) = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE- _PREFIX1_ RTT))

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice

!     THERMODYNAMICAL FUNCTIONS .

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEW,FOEDE,FOEDESU,FOELH,FOELDCP
FOEEW ( PTARE ) =  _PREFIX2_ R2ES*EXP (&
  &( _PREFIX2_ R3LES*FOEDELTA(PTARE)+ _PREFIX2_ R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE- _PREFIX1_ RTT)&
&/ (PTARE-( _PREFIX2_ R4LES*FOEDELTA(PTARE)+ _PREFIX2_ R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))

FOEDE ( PTARE ) = &
  &(FOEDELTA(PTARE)* _PREFIX2_ R5ALVCP+(1.0_JPRB-FOEDELTA(PTARE))* _PREFIX2_ R5ALSCP)&
&/ (PTARE-( _PREFIX2_ R4LES*FOEDELTA(PTARE)+ _PREFIX2_ R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOEDESU ( PTARE ) = &
  &(FOEDELTA(PTARE)* _PREFIX2_ R5LES+(1.0_JPRB-FOEDELTA(PTARE))* _PREFIX2_ R5IES)&
&/ (PTARE-( _PREFIX2_ R4LES*FOEDELTA(PTARE)+ _PREFIX2_ R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOELH ( PTARE ) =&
         &FOEDELTA(PTARE)* _PREFIX1_ RLVTT + (1.0_JPRB-FOEDELTA(PTARE))* _PREFIX1_ RLSTT

FOELDCP ( PTARE ) = &
         &FOEDELTA(PTARE)* _PREFIX2_ RALVDCP + (1.0_JPRB-FOEDELTA(PTARE))* _PREFIX2_ RALSDCP

!     *****************************************************************

!           CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!     FOEALFA is calculated to distinguish the three cases:

!                       FOEALFA=1            water phase
!                       FOEALFA=0            ice phase
!                       0 < FOEALFA < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFA
FOEALFA (PTARE) = MIN(1.0_JPRB,((MAX( _PREFIX2_ RTICE,MIN( _PREFIX2_ RTWAT,PTARE))- _PREFIX2_ RTICE)&
 &* _PREFIX2_ RTWAT_RTICE_R)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWM,FOEDEM,FOELDCPM,FOELHM,FOE_DEWM_DT
FOEEWM ( PTARE ) =  _PREFIX2_ R2ES *&
     &(FOEALFA(PTARE)*EXP( _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES))+&
  &(1.0_JPRB-FOEALFA(PTARE))*EXP( _PREFIX2_ R3IES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4IES)))

FOE_DEWM_DT( PTARE ) =  _PREFIX2_ R2ES * ( &
     &  _PREFIX2_ R3LES*FOEALFA(PTARE)*EXP( _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES)) &
     &    *( _PREFIX1_ RTT- _PREFIX2_ R4LES)/(PTARE- _PREFIX2_ R4LES)**2 + &
     &  _PREFIX2_ R3IES*(1.0-FOEALFA(PTARE))*EXP( _PREFIX2_ R3IES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4IES)) &
     &    *( _PREFIX1_ RTT- _PREFIX2_ R4IES)/(PTARE- _PREFIX2_ R4IES)**2)

FOEDEM ( PTARE ) = FOEALFA(PTARE)* _PREFIX2_ R5ALVCP*(1.0_JPRB/(PTARE- _PREFIX2_ R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))* _PREFIX2_ R5ALSCP*(1.0_JPRB/(PTARE- _PREFIX2_ R4IES)**2)

FOELDCPM ( PTARE ) = FOEALFA(PTARE)* _PREFIX2_ RALVDCP+&
            &(1.0_JPRB-FOEALFA(PTARE))* _PREFIX2_ RALSDCP

FOELHM ( PTARE ) =&
         &FOEALFA(PTARE)* _PREFIX1_ RLVTT+(1.0_JPRB-FOEALFA(PTARE))* _PREFIX1_ RLSTT


!     Temperature normalization for humidity background change of variable
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOETB
FOETB ( PTARE )=FOEALFA(PTARE)* _PREFIX2_ R3LES*( _PREFIX1_ RTT- _PREFIX2_ R4LES)*(1.0_JPRB/(PTARE- _PREFIX2_ R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))* _PREFIX2_ R3IES*( _PREFIX1_ RTT- _PREFIX2_ R4IES)*(1.0_JPRB/(PTARE- _PREFIX2_ R4IES)**2)

!     ------------------------------------------------------------------
!     *****************************************************************

!           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

!     *****************************************************************

!     FOEALFCU is calculated to distinguish the three cases:

!                       FOEALFCU=1            water phase
!                       FOEALFCU=0            ice phase
!                       0 < FOEALFCU < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFCU 
FOEALFCU (PTARE) = MIN(1.0_JPRB,((MAX( _PREFIX2_ RTICECU,MIN( _PREFIX2_ RTWAT,PTARE))&
&- _PREFIX2_ RTICECU)* _PREFIX2_ RTWAT_RTICECU_R)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWMCU,FOEDEMCU,FOELDCPMCU,FOELHMCU
FOEEWMCU ( PTARE ) =  _PREFIX2_ R2ES *&
     &(FOEALFCU(PTARE)*EXP( _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES))+&
  &(1.0_JPRB-FOEALFCU(PTARE))*EXP( _PREFIX2_ R3IES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4IES)))

FOEDEMCU ( PTARE )=FOEALFCU(PTARE)* _PREFIX2_ R5ALVCP*(1.0_JPRB/(PTARE- _PREFIX2_ R4LES)**2)+&
             &(1.0_JPRB-FOEALFCU(PTARE))* _PREFIX2_ R5ALSCP*(1.0_JPRB/(PTARE- _PREFIX2_ R4IES)**2)

FOELDCPMCU ( PTARE ) = FOEALFCU(PTARE)* _PREFIX2_ RALVDCP+&
            &(1.0_JPRB-FOEALFCU(PTARE))* _PREFIX2_ RALSDCP

FOELHMCU ( PTARE ) =&
         &FOEALFCU(PTARE)* _PREFIX1_ RLVTT+(1.0_JPRB-FOEALFCU(PTARE))* _PREFIX1_ RLSTT
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
!      unlike the FOEE functions does not include 1/( _PREFIX1_ RETV+1.0_JPRB) factor

REAL(KIND=JPRB) :: FOEEWMO, FOEELIQ, FOEEICE, FOELSON 
FOEEWMO( PTARE ) =  _PREFIX2_ R2ES*EXP( _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES))
FOEELIQ( PTARE ) =  _PREFIX2_ R2ES*EXP( _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES))
FOEEICE( PTARE ) =  _PREFIX2_ R2ES*EXP( _PREFIX2_ R3IES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4IES))
FOELSON( PTARE ) = EXP( -6096.9385_JPRB/PTARE + 21.2409642_JPRB &
	             - 2.711193E-2_JPRB * PTARE    &
                     + 1.673952E-5_JPRB * PTARE**2 &
		     + 2.433502_JPRB * LOG(PTARE))

REAL(KIND=JPRB) :: FOEEWM_V,FOEEWMCU_V,FOELES_V,FOEIES_V
REAL(KIND=JPRB) :: EXP1,EXP2
      FOELES_V(PTARE)= _PREFIX2_ R3LES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4LES)
      FOEIES_V(PTARE)= _PREFIX2_ R3IES*(PTARE- _PREFIX1_ RTT)/(PTARE- _PREFIX2_ R4IES)
      FOEEWM_V( PTARE,EXP1,EXP2 )= _PREFIX2_ R2ES*(FOEALFA(PTARE)*EXP1+ &
          & (1.0_JPRB-FOEALFA(PTARE))*EXP2)
      FOEEWMCU_V ( PTARE,EXP1,EXP2 ) =  _PREFIX2_ R2ES*(FOEALFCU(PTARE)*EXP1+&
          &(1.0_JPRB-FOEALFCU(PTARE))*EXP2)

