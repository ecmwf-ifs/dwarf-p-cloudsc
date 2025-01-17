! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, &
     & KFLDX, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
     & YDOMCST, YDOETHF, YDECLDP )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    USE YOECLDP  , ONLY : TECLDP
    USE YOMCST   , ONLY : TOMCST
    USE YOETHF   , ONLY : TOETHF

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    INTEGER(KIND=JPIM), INTENT(IN)    :: KFLDX
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_TMP(:) ! cumulative tendency used as input
    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFA(:,:,:)  ! CC from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFL(:,:,:)  ! Liq from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFI(:,:,:)  ! Ice from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNA(:,:,:) ! CC from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNL(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNI(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRSW(:,:,:) ! Short-wave heating rate
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRLW(:,:,:) ! Long-wave heating rate
    REAL(KIND=JPRB),    INTENT(IN)    :: PVERVEL(:,:,:) !Vertical velocity
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLSM(:,:)    ! Land fraction (0-1)
    LOGICAL        ,    INTENT(IN)    :: LDCUM(:,:)   ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN)    :: KTYPE(:,:)   ! Convection type 0,1,2
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PSNDE(:,:,:) ! Conv. detrained snow
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: PLCRIT_AER(:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PICRIT_AER(:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PRE_ICE(:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB),    INTENT(IN)    :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PRAINFRAC_TOPRFZ(:,:) 
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLF(:,:,:)  ! Flux of liquid
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQIF(:,:,:)  ! Flux of ice
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQLNG(:,:,:) ! -ve corr for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQNNG(:,:,:) ! -ve corr for ice
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQRF(:,:,:)  ! Flux diagnostics
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQSF(:,:,:)  !    for DDH, generic
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQRNG(:,:,:) ! rain
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQSNG(:,:,:) ! snow
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLTUR(:,:,:) ! liquid flux due to VDF
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQITUR(:,:,:) ! ice flux due to VDF
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    TYPE(TOMCST)    :: YDOMCST
    TYPE(TOETHF)    :: YDOETHF
    TYPE(TECLDP)    :: YDECLDP

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         CALL CLOUDSC &
              & (    1,    ICEND,    NPROMA,  NLEV,&
              & PTSPHY,&
              & PT(:,:,IBL), PQ(:,:,IBL), &
              & TENDENCY_TMP(IBL)%T, TENDENCY_TMP(IBL)%Q, TENDENCY_TMP(IBL)%A, TENDENCY_TMP(IBL)%CLD, &
              & TENDENCY_LOC(IBL)%T, TENDENCY_LOC(IBL)%Q, TENDENCY_LOC(IBL)%A, TENDENCY_LOC(IBL)%CLD, &
              & PVFA(:,:,IBL), PVFL(:,:,IBL), PVFI(:,:,IBL), PDYNA(:,:,IBL), PDYNL(:,:,IBL), PDYNI(:,:,IBL), &
              & PHRSW(:,:,IBL),    PHRLW(:,:,IBL),&
              & PVERVEL(:,:,IBL),  PAP(:,:,IBL),      PAPH(:,:,IBL),&
              & PLSM(:,IBL),       LDCUM(:,IBL),      KTYPE(:,IBL), &
              & PLU(:,:,IBL),      PLUDE(:,:,IBL),    PSNDE(:,:,IBL),    PMFU(:,:,IBL),     PMFD(:,:,IBL),&
              !---prognostic fields
              & PA(:,:,IBL),       PCLV(:,:,:,IBL),   PSUPSAT(:,:,IBL),&
              !-- arrays for aerosol-cloud interactions
              & PLCRIT_AER(:,:,IBL),PICRIT_AER(:,:,IBL),&
              & PRE_ICE(:,:,IBL),&
              & PCCN(:,:,IBL),     PNICE(:,:,IBL),&
              !---diagnostic output
              & PCOVPTOT(:,:,IBL), PRAINFRAC_TOPRFZ(:,IBL),&
              !---resulting fluxes
              & PFSQLF(:,:,IBL),   PFSQIF (:,:,IBL),  PFCQNNG(:,:,IBL),  PFCQLNG(:,:,IBL),&
              & PFSQRF(:,:,IBL),   PFSQSF (:,:,IBL),  PFCQRNG(:,:,IBL),  PFCQSNG(:,:,IBL),&
              & PFSQLTUR(:,:,IBL), PFSQITUR (:,:,IBL), &
              & PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL),   PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL),&
              & KFLDX, &
              & YDOMCST, YDOETHF, YDECLDP)

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)
    
  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
