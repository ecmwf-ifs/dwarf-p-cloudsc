! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_INIT_MOD

CONTAINS

SUBROUTINE CLOUDSC_INIT &
 !---input
 & (KIDIA,    KFDIA,    KLON,    KLEV,&
 & TENDENCY_LOC_T, TENDENCY_LOC_Q, TENDENCY_LOC_A, TENDENCY_LOC_CLD)

!===============================================================================
!**** *CLOUDSC* -  ROUTINE FOR PARAMATERIZATION OF CLOUD PROCESSES
!                  FOR PROGNOSTIC CLOUD SCHEME
!!
!     M.Tiedtke, C.Jakob, A.Tompkins, R.Forbes     (E.C.M.W.F.)
!!
!     PURPOSE
!     -------
!          THIS ROUTINE UPDATES THE CONV/STRAT CLOUD FIELDS.
!          THE FOLLOWING PROCESSES ARE CONSIDERED:
!        - Detrainment of cloud water from convective updrafts
!        - Evaporation/condensation of cloud water in connection
!           with heating/cooling such as by subsidence/ascent
!        - Erosion of clouds by turbulent mixing of cloud air
!           with unsaturated environmental air
!        - Deposition onto ice when liquid water present (Bergeron-Findeison) 
!        - Conversion of cloud water into rain (collision-coalescence)
!        - Conversion of cloud ice to snow (aggregation)
!        - Sedimentation of rain, snow and ice
!        - Evaporation of rain and snow
!        - Melting of snow and ice
!        - Freezing of liquid and rain
!        Note: Turbulent transports of s,q,u,v at cloud tops due to
!           buoyancy fluxes and lw radiative cooling are treated in 
!           the VDF scheme
!!
!     INTERFACE.
!     ----------
!          *CLOUDSC* IS CALLED FROM *CALLPAR*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!     T,Q,L,PHI AND DETRAINMENT OF CLOUD WATER FROM THE
!     CONVECTIVE CLOUDS (MASSFLUX CONVECTION SCHEME), BOUNDARY
!     LAYER TURBULENT FLUXES OF HEAT AND MOISTURE, RADIATIVE FLUXES,
!     OMEGA.
!     IT RETURNS ITS OUTPUT TO:
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES T AND Q
!        AS WELL AS CLOUD VARIABLES L AND C
!      2.GENERATES PRECIPITATION FLUXES FROM STRATIFORM CLOUDS
!!
!     EXTERNALS.
!     ----------
!          NONE
!!
!     MODIFICATIONS.
!     -------------
!      M. TIEDTKE    E.C.M.W.F.     8/1988, 2/1990
!     CH. JAKOB      E.C.M.W.F.     2/1994 IMPLEMENTATION INTO IFS
!     A.TOMPKINS     E.C.M.W.F.     2002   NEW NUMERICS
!        01-05-22 : D.Salmond   Safety modifications
!        02-05-29 : D.Salmond   Optimisation
!        03-01-13 : J.Hague     MASS Vector Functions  J.Hague
!        03-10-01 : M.Hamrud    Cleaning
!        04-12-14 : A.Tompkins  New implicit solver and physics changes
!        04-12-03 : A.Tompkins & M.Ko"hler  moist PBL
!     G.Mozdzynski  09-Jan-2006  EXP security fix
!        19-01-09 : P.Bechtold  Changed increased RCLDIFF value for KTYPE=2
!        07-07-10 : A.Tompkins/R.Forbes  4-Phase flexible microphysics
!        01-03-11 : R.Forbes    Mixed phase changes and tidy up
!        01-10-11 : R.Forbes    Melt ice to rain, allow rain to freeze
!        01-10-11 : R.Forbes    Limit supersat to avoid excessive values
!        31-10-11 : M.Ahlgrimm  Add rain, snow and PEXTRA to DDH output
!        17-02-12 : F.Vana      Simplified/optimized LU factorization
!        18-05-12 : F.Vana      Cleaning + better support of sequential physics
!        N.Semane+P.Bechtold     04-10-2012 Add RVRFACTOR factor for small planet
!        01-02-13 : R.Forbes    New params of autoconv/acc,rain evap,snow riming
!        15-03-13 : F. Vana     New dataflow + more tendencies from the first call
!        K. Yessad (July 2014): Move some variables.
!        F. Vana  05-Mar-2015  Support for single precision
!        15-01-15 : R.Forbes    Added new options for snow evap & ice deposition
!        10-01-15 : R.Forbes    New physics for rain freezing
!        23-10-14 : P. Bechtold remove zeroing of convection arrays
!
!     SWITCHES.
!     --------
!!
!     MODEL PARAMETERS
!     ----------------
!     RCLDIFF:    PARAMETER FOR EROSION OF CLOUDS
!     RCLCRIT_SEA:  THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER SEA
!     RCLCRIT_LAND: THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER LAND
!     RLCRITSNOW: THRESHOLD VALUE FOR SNOW AUTOCONVERSION
!     RKCONV:     PARAMETER FOR AUTOCONVERSION OF CLOUDS (KESSLER)
!     RCLDMAX:    MAXIMUM POSSIBLE CLW CONTENT (MASON,1971)
!!
!     REFERENCES.
!     ----------
!     TIEDTKE MWR 1993
!     JAKOB PhD 2000
!     GREGORY ET AL. QJRMS 2000
!     TOMPKINS ET AL. QJRMS 2007
!!
!===============================================================================

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMPHYDER ,ONLY : STATE_TYPE
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV  
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
USE YOECLDP  , ONLY : TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
#ifndef CLOUDSC_STMT_FUNC
USE FCTTRE_MOD, ONLY: FOEDELTA, FOEALFA, FOEEWM, FOEEICE, FOEELIQ, FOELDCP, FOELDCPM, FOEDEM
USE FCCLD_MOD, ONLY : FOKOOP
#endif

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
! PLCRIT_AER : critical liquid mmr for rain autoconversion process
! PICRIT_AER : critical liquid mmr for snow autoconversion process
! PRE_LIQ : liq Re
! PRE_ICE : ice Re
! PCCN    : liquid cloud condensation nuclei
! PNICE   : ice number concentration (cf. CCN)

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: TENDENCY_LOC_T(KLON, KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: TENDENCY_LOC_Q(KLON, KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: TENDENCY_LOC_A(KLON, KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: TENDENCY_LOC_CLD(KLON, KLEV, NCLV)

INTEGER(KIND=JPIM) :: JM, JK, JL

#ifdef CLOUDSC_STMT_FUNC
#include "fcttre.func.h"
#include "fccld.func.h"
#endif

! -----------------------------------------------
! INITIALIZATION OF OUTPUT TENDENCIES
! -----------------------------------------------
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    TENDENCY_LOC_T(JL,JK)=0.0_JPRB
    TENDENCY_LOC_Q(JL,JK)=0.0_JPRB
    TENDENCY_LOC_A(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JM=1,NCLV-1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      TENDENCY_LOC_CLD(JL,JK,JM)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE CLOUDSC_INIT

END MODULE CLOUDSC_INIT_MOD
