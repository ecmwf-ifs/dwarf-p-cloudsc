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
  USE EC_PMON_MOD, ONLY: EC_PMON
  USE CLOUDSC_GLOBAL_ATLAS_STATE_MOD, ONLY: CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW
  
  USE ATLAS_MODULE
  USE, INTRINSIC :: ISO_C_BINDING
  USE ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS_MODULE

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOTG, KFLDX, PTSPHY)
  
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NGPTOTG, KFLDX
    REAL(KIND=JPRB), INTENT(IN)   :: PTSPHY       ! Physics timestep

    TYPE(CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW) :: FBLOCK
    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    TYPE(ATLAS_FIELD) :: FIELD
    INTEGER(KIND=JPIM)    :: NPROMA, NLEV, NGPTOT

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    INTEGER(KIND=JPIB) :: ENERGY, POWER, POWER_TOTAL, POWER_MAX, POWER_COUNT
    LOGICAL            :: LEC_PMON = .FALSE.
    CHARACTER(LEN=1)   :: CLEC_PMON
    TYPE(ATLAS_TRACE)  :: TRACE

    TRACE = ATLAS_TRACE("cloudsc_driver_mod.F90", __LINE__, "CLOUDSC_DRIVER","COMPUTE")

    CALL GET_ENVIRONMENT_VARIABLE('EC_PMON', CLEC_PMON)
    IF (CLEC_PMON == '1') LEC_PMON = .TRUE.

    POWER_MAX = 0_JPIB
    POWER_TOTAL = 0_JPIB
    POWER_COUNT = 0_JPIB

    FIELD = FSET%FIELD("PCLV")
    FSPACE = FIELD%FUNCTIONSPACE()
    NPROMA = FIELD%SHAPE(1)
    NLEV = FSPACE%LEVELS()
    NGPTOT = FSPACE%SIZE()

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGTOT=', i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOT, NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID,energy,power,FBLOCK) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime) reduction(+:power_total,power_count) reduction(max:power_max)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         ! get block views
         call FBLOCK%GET_BLOCK(FSET, IBL)

         !-- These were uninitialized : meaningful only when we compare error differences
         FBLOCK%PCOVPTOT(:,:) = 0.0_JPRB
         FBLOCK%TENDENCY_LOC%cld(:,:,NCLV) = 0.0_JPRB

         CALL CLOUDSC &
              & (    1,    ICEND,    NPROMA,  NLEV,&
              & PTSPHY,&
              & FBLOCK%PT, FBLOCK%PQ, FBLOCK%TENDENCY_CML, FBLOCK%TENDENCY_TMP, FBLOCK%TENDENCY_LOC, &
              & FBLOCK%PVFA, FBLOCK%PVFL, FBLOCK%PVFI, FBLOCK%PDYNA, FBLOCK%PDYNL, FBLOCK%PDYNI, &
              & FBLOCK%PHRSW, FBLOCK%PHRLW,&
              & FBLOCK%PVERVEL, FBLOCK%PAP, FBLOCK%PAPH,&
              & FBLOCK%PLSM, FBLOCK%LDCUM, FBLOCK%KTYPE, &
              & FBLOCK%PLU, FBLOCK%PLUDE, &
              & FBLOCK%PSNDE, FBLOCK%PMFU, FBLOCK%PMFD,&
              !---prognostic fields
              & FBLOCK%PA, FBLOCK%PCLV, FBLOCK%PSUPSAT,&
              !-- arrays for aerosol-cloud interactions
              & FBLOCK%PLCRIT_AER, FBLOCK%PICRIT_AER,&
              & FBLOCK%PRE_ICE,&
              & FBLOCK%PCCN, FBLOCK%PNICE,&
              !---diagnostic output
              & FBLOCK%PCOVPTOT, FBLOCK%PRAINFRAC_TOPRFZ,&
              !---resulting fluxes
              & FBLOCK%PFSQLF,   FBLOCK%PFSQIF,  FBLOCK%PFCQNNG,  FBLOCK%PFCQLNG,&
              & FBLOCK%PFSQRF,   FBLOCK%PFSQSF,  FBLOCK%PFCQRNG,  FBLOCK%PFCQSNG,&
              & FBLOCK%PFSQLTUR, FBLOCK%PFSQITUR, &
              & FBLOCK%PFPLSL,   FBLOCK%PFPLSN,  FBLOCK%PFHPSL,   FBLOCK%PFHPSN,&
              & KFLDX)

          !--- end of a future plan to replace the call to CLOUDSC ------ 

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

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

      IF (LEC_PMON) THEN
        print *, "Power usage (sampled):: max: ", POWER_MAX, "avg:", &
         & (REAL(POWER_TOTAL, KIND=JPRD) / REAL(POWER_COUNT, KIND=JPRD)), &
         & "count:", POWER_COUNT
      END IF

      CALL FIELD%FINAL()
      CALL FSPACE%FINAL()
      CALL TRACE%FINAL()

  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
