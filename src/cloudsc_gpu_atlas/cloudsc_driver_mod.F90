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

#INCLUDE "cloudsc.F90"

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

    !$acc routine(CLOUDSC)

    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PT(:,:,:)       ! T at start of callpar
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PQ(:,:,:)       ! Q at start of callpar
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PVFA(:,:,:)     ! CC from VDF scheme
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PVFL(:,:,:)     ! Liq from VDF scheme
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PVFI(:,:,:)     ! Ice from VDF scheme
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PDYNA(:,:,:)    ! CC from Dynamics
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PDYNL(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PDYNI(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PHRSW(:,:,:)    ! Short-wave heating rate
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PHRLW(:,:,:)    ! Long-wave heating rate
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PVERVEL(:,:,:)  ! Vertical velocity
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PAP(:,:,:)      ! Pressure on full levels
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PAPH(:,:,:)     ! Pressure on half levels
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLSM(:,:)       ! Land fraction (0-1)
    LOGICAL, POINTER, CONTIGUOUS          :: LDCUM(:,:)      ! Convection active
    INTEGER(KIND=JPIM),POINTER, CONTIGUOUS:: KTYPE(:,:)      ! Convection type 0,1,2
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLU(:,:,:)      ! Conv. condensate
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLUDE(:,:,:)    ! Conv. detrained water
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PSNDE(:,:,:)    ! Conv. detrained snow
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PMFU(:,:,:)     ! Conv. mass flux up
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PMFD(:,:,:)     ! Conv. mass flux down
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PA(:,:,:)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCLV(:,:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLCRIT_AER(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PICRIT_AER(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PRE_ICE(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCOVPTOT(:,:,:)    ! Precip fraction
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PRAINFRAC_TOPRFZ(:,:)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQLF(:,:,:)    ! Flux of liquid
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQIF(:,:,:)    ! Flux of ice
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFCQLNG(:,:,:)   ! -ve corr for liq
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFCQNNG(:,:,:)   ! -ve corr for ice
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQRF(:,:,:)    ! Flux diagnostics
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQSF(:,:,:)    !    for DDH, generic
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFCQRNG(:,:,:)   ! rain
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFCQSNG(:,:,:)   ! snow
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQLTUR(:,:,:)  ! liquid flux due to VDF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFSQITUR(:,:,:)  ! ice flux due to VDF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFPLSL(:,:,:)    ! liq+rain sedim flux
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFPLSN(:,:,:)    ! ice+snow sedim flux
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFHPSL(:,:,:)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PFHPSN(:,:,:)    ! ice number concentration (cf. CCN)

    TYPE(STATE_TYPE), POINTER, DIMENSION(:) :: TEND_CML
    TYPE(STATE_TYPE), POINTER, DIMENSION(:) :: TEND_LOC
    TYPE(STATE_TYPE), POINTER, DIMENSION(:) :: TEND_TMP

    REAL(KIND=JPRB), POINTER :: TMP3D(:,:,:,:)

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

    CALL FSET%DATA(1, PLCRIT_AER)
    CALL FSET%DATA(2, PICRIT_AER)
    CALL FSET%DATA(3, PRE_ICE)
    CALL FSET%DATA(4, PCCN)
    CALL FSET%DATA(5, PNICE)
    CALL FSET%DATA(6, PT)
    CALL FSET%DATA(7, PQ)
    CALL FSET%DATA(8, PVFA)
    CALL FSET%DATA(9, PVFL)
    CALL FSET%DATA(10, PVFI)
    CALL FSET%DATA(11, PDYNA)
    CALL FSET%DATA(12, PDYNL)
    CALL FSET%DATA(13, PDYNI)
    CALL FSET%DATA(14, PHRSW)
    CALL FSET%DATA(15, PHRLW)
    CALL FSET%DATA(16, PVERVEL)
    CALL FSET%DATA(17, PAP)
    CALL FSET%DATA(18, PLU)
    CALL FSET%DATA(19, PLUDE)
    CALL FSET%DATA(20, PSNDE)
    CALL FSET%DATA(21, PMFU)
    CALL FSET%DATA(22, PMFD)
    CALL FSET%DATA(23, PA)
    CALL FSET%DATA(24, PSUPSAT)
    CALL FSET%DATA(25, PLSM)
    CALL FSET%DATA(26, LDCUM)
    CALL FSET%DATA(27, KTYPE)
    CALL FSET%DATA(28, PAPH)
    CALL FSET%DATA(29, PCLV)

    CALL FSET%DATA(30, TMP3D)

    DO JKGLO = 1, FSPACE%NBLKS()
        TEND_CML(JKGLO)%T   => TMP3D(:,:,1,JKGLO)
        TEND_CML(JKGLO)%A   => TMP3D(:,:,2,JKGLO)
        TEND_CML(JKGLO)%Q   => TMP3D(:,:,3,JKGLO)
        TEND_CML(JKGLO)%CLD => TMP3D(:,:,4:,JKGLO)
    END DO
    CALL FSET%DATA(31, TMP3D)
    DO JKGLO = 1, FSPACE%NBLKS()
        TEND_TMP(JKGLO)%T   => TMP3D(:,:,1,JKGLO)
        TEND_TMP(JKGLO)%A   => TMP3D(:,:,2,JKGLO)
        TEND_TMP(JKGLO)%Q   => TMP3D(:,:,3,JKGLO)
        TEND_TMP(JKGLO)%CLD => TMP3D(:,:,4:,JKGLO)
    END DO
    CALL FSET%DATA(32, TMP3D)
    DO JKGLO = 1, FSPACE%NBLKS()
        TEND_LOC(JKGLO)%T   => TMP3D(:,:,1,JKGLO)
        TEND_LOC(JKGLO)%A   => TMP3D(:,:,2,JKGLO)
        TEND_LOC(JKGLO)%Q   => TMP3D(:,:,3,JKGLO)
        TEND_LOC(JKGLO)%CLD => TMP3D(:,:,4:,JKGLO)
    END DO

    CALL FSET%DATA(33, PFSQLF)
    CALL FSET%DATA(34, PFSQIF)
    CALL FSET%DATA(35, PFCQLNG)
    CALL FSET%DATA(36, PFCQNNG)
    CALL FSET%DATA(37, PFSQRF)
    CALL FSET%DATA(38, PFSQSF)
    CALL FSET%DATA(39, PFCQRNG)
    CALL FSET%DATA(40, PFCQSNG)
    CALL FSET%DATA(41, PFSQLTUR)
    CALL FSET%DATA(42, PFSQITUR)
    CALL FSET%DATA(43, PFPLSL)
    CALL FSET%DATA(44, PFPLSN)
    CALL FSET%DATA(45, PFHPSL)
    CALL FSET%DATA(46, PFHPSN)
    CALL FSET%DATA(47, PCOVPTOT)
    CALL FSET%DATA(48, PRAINFRAC_TOPRFZ)

    DO JKGLO = 1, 48
!        IF (JKGLO > 29 .OR. JKGLO < 33) CONTINUE
        FIELD = FSET%FIELD(JKGLO)
        CALL FIELD%UPDATE_DEVICE()
    END DO

!!    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID,energy,power,FBLOCK) &
!!    !$omp& num_threads(NUMOMP)
!$acc data deviceptr(PT, PQ,TEND_CML,TEND_LOC,&
!$acc & TEND_TMP, PVFA, PVFL, PVFI, &
!$acc & PDYNA, PDYNL, PDYNI, PHRSW,    PHRLW,&
!$acc & PVERVEL,  PAP,      PAPH,&
!$acc & PLSM,       LDCUM,      KTYPE, &
!$acc & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD,&
!$acc & PA,       PCLV,   PSUPSAT,&
!$acc & PLCRIT_AER,PICRIT_AER,&
!$acc & PRE_ICE, PCCN,     PNICE,&
!$acc & PCOVPTOT, PRAINFRAC_TOPRFZ,&
!$acc & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,&
!$acc & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG,&
!$acc & PFSQLTUR, PFSQITUR , &
!$acc & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN)

    ! Local timer for each thread
    !TID = GET_THREAD_NUM()
    !CALL TIMER%THREAD_START(TID)

!!    !$omp do schedule(runtime) reduction(+:power_total,power_count) reduction(max:power_max)
!$acc parallel loop gang vector_length(NPROMA)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         ! get block views
         !call FBLOCK%GET_BLOCK(FSET, IBL)

         !-- These were uninitialized : meaningful only when we compare error differences
         !FBLOCK%PCOVPTOT(:,:) = 0.0_JPRB
         !FBLOCK%TENDENCY_LOC%cld(:,:,NCLV) = 0.0_JPRB

          CALL CLOUDSC &
                & (1, ICEND, NPROMA, NLEV, PTSPHY,&
                & PT(:,:,IBL), PQ(:,:,IBL), &
                & TEND_CML(IBL), TEND_TMP(IBL), TEND_LOC(IBL), &
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
                & KFLDX )

          !--- end of a future plan to replace the call to CLOUDSC ------ 

         IF (LEC_PMON) THEN
           ! Sample power consuption
           IF (MOD(IBL, 100) == 0) THEN
             !CALL EC_PMON(ENERGY, POWER)
             POWER_MAX = MAX(POWER_MAX, POWER)
             POWER_TOTAL = POWER_TOTAL + POWER
             POWER_COUNT = POWER_COUNT + 1
           END IF
         END IF

         ! Log number of columns processed by this thread
         !CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO
!$acc end parallel loop

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
!!     !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

!!      !$omp end parallel
!$acc end data

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
