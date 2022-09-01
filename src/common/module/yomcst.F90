! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCST

USE PARKIND1,    ONLY : JPRB

IMPLICIT NONE

!!** SAVE

!     ------------------------------------------------------------------

!*    Common of physical constants
!     You will find the meanings in the annex 1 of the documentation

! A1.0 Fundamental constants
! * RPI          : number Pi
! * RCLUM        : light velocity
! * RHPLA        : Planck constant
! * RKBOL        : Bolzmann constant
! * RNAVO        : Avogadro number
!**REAL(KIND=JPRB) :: RPI
!**REAL(KIND=JPRB) :: RCLUM
!**REAL(KIND=JPRB) :: RHPLA
!**REAL(KIND=JPRB) :: RKBOL
!**REAL(KIND=JPRB) :: RNAVO
!**
!**! A1.1 Astronomical constants
!**! * RDAY         : duration of the solar day
!**! * RDAYI        : invariant time unit of 86400s
!**! * RHOUR        : duration of the solar hour 
!**! * REA          : astronomical unit (mean distance Earth-sun)
!**! * REPSM        : polar axis tilting angle
!**! * RSIYEA       : duration of the sideral year
!**! * RSIDAY       : duration of the sideral day
!**! * ROMEGA       : angular velocity of the Earth rotation
!**REAL(KIND=JPRB) :: RDAY
!**REAL(KIND=JPRB) :: RDAYI
!**REAL(KIND=JPRB) :: RHOUR
!**REAL(KIND=JPRB) :: REA
!**REAL(KIND=JPRB) :: REPSM
!**REAL(KIND=JPRB) :: RSIYEA
!**REAL(KIND=JPRB) :: RSIDAY
!**REAL(KIND=JPRB) :: ROMEGA
!**
!**! A1.2 Geoide
!**! * RA           : Earth radius
!**! * RG           : gravity constant
!**! * R1SA         : 1/RA
!**REAL(KIND=JPRB) :: RA
REAL(KIND=JPRB) :: RG, rg_h
!**REAL(KIND=JPRB) :: R1SA
!**
!**! A1.3 Radiation
!**! * RSIGMA       : Stefan-Bolzman constant
!**! * RI0          : solar constant
!**REAL(KIND=JPRB) :: RSIGMA
!**REAL(KIND=JPRB) :: RI0
!**
!**! A1.4 Thermodynamic gas phase
!**! * R            : perfect gas constant
!**! * RMD          : dry air molar mass
!**! * RMV          : vapour water molar mass
!**! * RMO3         : ozone molar mass
!**! * RD           : R_dry (dry air constant)
!**! * RV           : R_vap (vapour water constant)
!**! * RCPD         : Cp_dry (dry air calorific capacity at constant pressure)
!**! * RCPV         : Cp_vap (vapour calorific capacity at constant pressure)
!**! * RCVD         : Cv_dry (dry air calorific capacity at constant volume)
!**! * RCVV         : Cv_vap (vapour calorific capacity at constant volume)
!**! * RKAPPA       : Kappa = R_dry/Cp_dry
!**! * RETV         : R_vap/R_dry - 1
!**! * RMCO2        : CO2 (carbon dioxyde) molar mass
!**! * RMCH4        : CH4 (methane) molar mass
!**! * RMN2O        : N2O molar mass
!**! * RMCO         : CO (carbon monoxyde) molar mass
!**! * RMHCHO       : HCHO molar mass
!**! * RMNO2        : NO2 (nitrogen dioxyde) molar mass
!**! * RMSO2        : SO2 (sulfur dioxyde) molar mass
!**! * RMSO4        : SO4 (sulphate) molar mass
!**REAL(KIND=JPRB) :: R
!**REAL(KIND=JPRB) :: RMD
!**REAL(KIND=JPRB) :: RMV
!**REAL(KIND=JPRB) :: RMO3
REAL(KIND=JPRB) :: RD, rd_h
REAL(KIND=JPRB) :: RV, rv_h
REAL(KIND=JPRB) :: RCPD, rcpd_h
!**REAL(KIND=JPRB) :: RCPV
!**REAL(KIND=JPRB) :: RCVD
!**REAL(KIND=JPRB) :: RCVV
!**REAL(KIND=JPRB) :: RKAPPA
REAL(KIND=JPRB) :: RETV, retv_h
!**REAL(KIND=JPRB) :: RMCO2
!**REAL(KIND=JPRB) :: RMCH4
!**REAL(KIND=JPRB) :: RMN2O
!**REAL(KIND=JPRB) :: RMCO
!**REAL(KIND=JPRB) :: RMHCHO
!**REAL(KIND=JPRB) :: RMNO2
!**REAL(KIND=JPRB) :: RMSO2
!**REAL(KIND=JPRB) :: RMSO4
!**
!**! A1.5,6 Thermodynamic liquid,solid phases
!**! * RCW          : Cw (calorific capacity of liquid water)
!**! * RCS          : Cs (calorific capacity of solid water)
!**REAL(KIND=JPRB) :: RCW
!**REAL(KIND=JPRB) :: RCS
!**
!**! A1.7 Thermodynamic transition of phase
!**! * RATM         : pre_n = "normal" pressure
!**! * RTT          : Tt = temperature of water fusion at "pre_n"
!**! * RLVTT        : RLvTt = vaporisation latent heat at T=Tt
!**! * RLSTT        : RLsTt = sublimation latent heat at T=Tt
!**! * RLVZER       : RLv0 = vaporisation latent heat at T=0K
!**! * RLSZER       : RLs0 = sublimation latent heat at T=0K
!**! * RLMLT        : RLMlt = melting latent heat at T=Tt
!**! * RDT          : Tt - Tx(ew-ei)
!**REAL(KIND=JPRB) :: RATM
REAL(KIND=JPRB) :: RTT, rtt_h
REAL(KIND=JPRB) :: RLVTT, rlvtt_h
REAL(KIND=JPRB) :: RLSTT, rlstt_h
!**REAL(KIND=JPRB) :: RLVZER
!**REAL(KIND=JPRB) :: RLSZER
REAL(KIND=JPRB) :: RLMLT, rlmlt_h
!**REAL(KIND=JPRB) :: RDT

ATTRIBUTES(DEVICE) :: RG, RCPD, RD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV 
! A1.8 Curve of saturation
! * RESTT        : es(Tt) = saturation vapour tension at T=Tt
! * RGAMW        : Rgamw = (Cw-Cp_vap)/R_vap
! * RBETW        : Rbetw = RLvTt/R_vap + Rgamw*Tt
! * RALPW        : Ralpw = log(es(Tt)) + Rbetw/Tt + Rgamw*log(Tt)
! * RGAMS        : Rgams = (Cs-Cp_vap)/R_vap
! * RBETS        : Rbets = RLsTt/R_vap + Rgams*Tt
! * RALPS        : Ralps = log(es(Tt)) + Rbets/Tt + Rgams*log(Tt)
! * RALPD        : Ralpd = Ralps - Ralpw
! * RBETD        : Rbetd = Rbets - Rbetw
! * RGAMD        : Rgamd = Rgams - Rgamw
!**REAL(KIND=JPRB) :: RESTT
!**REAL(KIND=JPRB) :: RGAMW
!**REAL(KIND=JPRB) :: RBETW
!**REAL(KIND=JPRB) :: RALPW
!**REAL(KIND=JPRB) :: RGAMS
!**REAL(KIND=JPRB) :: RBETS
!**REAL(KIND=JPRB) :: RALPS
!**REAL(KIND=JPRB) :: RALPD
!**REAL(KIND=JPRB) :: RBETD
!**REAL(KIND=JPRB) :: RGAMD
!**
!**! NaN value
!**! CHARACTER(LEN=8), PARAMETER :: CSNAN = &
!**!   & CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(244)//CHAR(127)
!**REAL(KIND=JPRB) :: RSNAN
!**
!$acc declare copyin(rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv)
!$omp declare target(rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv)


!    ------------------------------------------------------------------

!!CONTAINS
!!
!!attributes(host)  SUBROUTINE YOMCST_LOAD_PARAMETERS()
!!USE FILE_IO_MOD, ONLY : LOAD_SCALAR_real
!!    CALL LOAD_SCALAR_real('RG', RG_h)
!!    CALL LOAD_SCALAR_real('RD', RD_h)
!!    CALL LOAD_SCALAR_real('RCPD', RCPD_h)
!!    CALL LOAD_SCALAR_real('RETV', RETV_h)
!!    CALL LOAD_SCALAR_real('RLVTT', RLVTT_h)
!!    CALL LOAD_SCALAR_real('RLSTT', RLSTT_h)
!!    CALL LOAD_SCALAR_real('RLMLT', RLMLT_h)
!!    CALL LOAD_SCALAR_real('RTT', RTT_h)
!!    CALL LOAD_SCALAR_real('RV', RV_h)

!!    rg=rg_h; rd=rd_h; rcpd=rcpd_h; retv=retv_h; rlvtt=rlvtt_h; rlstt=rlstt_h; rlmlt=rlmlt_h; rtt=rtt_h; rv=rv_h
!!  END SUBROUTINE YOMCST_LOAD_PARAMETERS

END MODULE YOMCST
