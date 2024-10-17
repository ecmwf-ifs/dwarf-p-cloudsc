! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_FIELD_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  USE EC_PMON_MOD, ONLY: EC_PMON
  USE CLOUDSC_FIELD_STATE_MOD, ONLY: CLOUDSC_AUX_TYPE, CLOUDSC_FLUX_TYPE, CLOUDSC_STATE_TYPE

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_FIELD( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, KFLDX, PTSPHY, PAUX, FLUX, &
     & TENDENCY_TMP, TENDENCY_LOC, YDOMCST, YDOETHF, YDECLDP)
    ! Driver routine that invokes the optimized CLAW-based CLOUDSC GPU kernel
    
    USE YOECLDP  , ONLY : TECLDP
    USE YOMCST   , ONLY : TOMCST
    USE YOETHF   , ONLY : TOETHF

    INTEGER(KIND=JPIM)        ,INTENT(IN) :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    INTEGER(KIND=JPIM)        ,INTENT(IN) :: KFLDX 
    REAL(KIND=JPRB)           ,INTENT(IN) :: PTSPHY       ! PHYSICS TIMESTEP
    TYPE(CLOUDSC_AUX_TYPE)    ,INTENT(IN) :: PAUX
    TYPE(CLOUDSC_FLUX_TYPE)   ,INTENT(IN) :: FLUX
    TYPE(CLOUDSC_STATE_TYPE)  ,INTENT(IN) :: TENDENCY_TMP
    TYPE(CLOUDSC_STATE_TYPE)  ,INTENT(INOUT) :: TENDENCY_LOC
    TYPE(TOMCST)              , INTENT(INOUT) :: YDOMCST
    TYPE(TOETHF)              , INTENT(INOUT) :: YDOETHF
    TYPE(TECLDP)              , INTENT(INOUT) :: YDECLDP

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND, NGPBLKS
    
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    INTEGER(KIND=JPIB) :: ENERGY, POWER, POWER_TOTAL, POWER_MAX, POWER_COUNT
    LOGICAL            :: LEC_PMON = .FALSE.
    CHARACTER(LEN=1)   :: CLEC_PMON

    CALL GET_ENVIRONMENT_VARIABLE('EC_PMON', CLEC_PMON)
    IF (CLEC_PMON == '1') LEC_PMON = .TRUE.

    POWER_MAX = 0_JPIB
    POWER_TOTAL = 0_JPIB
    POWER_COUNT = 0_JPIB

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$loki data

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID,energy,power) &
    !$omp& num_threads(NUMOMP) firstprivate(PAUX, FLUX, TENDENCY_TMP, TENDENCY_LOC)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime) reduction(+:power_total,power_count) reduction(max:power_max)
    DO JKGLO=1,NGPTOT,NPROMA
        IBL=(JKGLO-1)/NPROMA+1
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

        CALL PAUX%UPDATE_VIEW(IBL)
        CALL FLUX%UPDATE_VIEW(IBL)
        CALL TENDENCY_LOC%UPDATE_VIEW(IBL)
        CALL TENDENCY_TMP%UPDATE_VIEW(IBL)
        
        !-- These were uninitialized : meaningful only when we compare error differences
        PAUX%PCOVPTOT = 0.0_JPRB
        TENDENCY_LOC%CLD(:,:,NCLV) = 0.0_JPRB



        CALL CLOUDSC(    1,    ICEND,    NPROMA,  NLEV, & ! These could also be accessed through FIELD_STATE
              & PTSPHY,&
              & PAUX%PT, PAUX%PQ, &
              & TENDENCY_TMP%T, TENDENCY_TMP%Q, TENDENCY_TMP%A, TENDENCY_TMP%CLD, &
              & TENDENCY_LOC%T, TENDENCY_LOC%Q, TENDENCY_LOC%A, TENDENCY_LOC%CLD, &
              & PAUX%PVFA, PAUX%PVFL, PAUX%PVFI, PAUX%PDYNA, PAUX%PDYNL, PAUX%PDYNI, &
              & PAUX%PHRSW,    PAUX%PHRLW,&
              & PAUX%PVERVEL,  PAUX%PAP,      PAUX%PAPH,&
              & PAUX%PLSM,     PAUX%LDCUM,    PAUX%KTYPE, &
              & PAUX%PLU,      PAUX%PLUDE,    PAUX%PSNDE,    PAUX%PMFU,     PAUX%PMFD,&
              !---prognostic fields
              & PAUX%PA,&
              & PAUX%PCLV,  &
              & PAUX%PSUPSAT,&
!             -- arrays for aerosol-cloud interactions
!             !! & PQAER,    KAER, &
              & PAUX%PLCRIT_AER,PAUX%PICRIT_AER,&
              & PAUX%PRE_ICE,&
              & PAUX%PCCN,     PAUX%PNICE,&
              !---diagnostic output
              & PAUX%PCOVPTOT, PAUX%PRAINFRAC_TOPRFZ,&
              !---resulting fluxes
              & FLUX%PFSQLF,   FLUX%PFSQIF ,  FLUX%PFCQNNG,  FLUX%PFCQLNG,&
              & FLUX%PFSQRF,   FLUX%PFSQSF ,  FLUX%PFCQRNG,  FLUX%PFCQSNG,&
              & FLUX%PFSQLTUR, FLUX%PFSQITUR , &
              & FLUX%PFPLSL,   FLUX%PFPLSN,   FLUX%PFHPSL,   FLUX%PFHPSN, KFLDX, &
              & YDOMCST, YDOETHF, YDECLDP)

        IF (LEC_PMON) THEN
          ! Sample power consuption
          IF (MOD(IBL, 100) == 0) THEN
            CALL EC_PMON(ENERGY, POWER)
            POWER_MAX = MAX(POWER_MAX, POWER)
            POWER_TOTAL = POWER_TOTAL + POWER
            POWER_COUNT = POWER_COUNT + 1
          END IF
        END IF

        ! Log number of columns processed by this thread
        CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      !$loki end data

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

      IF (LEC_PMON) THEN
        print *, "Power usage (sampled):: max: ", POWER_MAX, "avg:", &
         & (REAL(POWER_TOTAL, KIND=JPRD) / REAL(POWER_COUNT, KIND=JPRD)), &
         & "count:", POWER_COUNT
      END IF
    
  END SUBROUTINE CLOUDSC_DRIVER_FIELD

END MODULE CLOUDSC_DRIVER_FIELD_MOD

