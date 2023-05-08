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

  SUBROUTINE CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOT, NGPTOTG, KFLDX, PTSPHY)
  
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NGPTOT, NGPTOTG, KFLDX
    REAL(KIND=JPRB), INTENT(IN)   :: PTSPHY       ! Physics timestep

    TYPE(CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW) :: FBLOCK
    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    TYPE(ATLAS_FIELD) :: FIELD
    INTEGER(KIND=JPIM)    :: NPROMA, NLEV
    REAL(KIND=JPRB), POINTER   :: PT(:,:)    ! T at start of callpar
    REAL(KIND=JPRB), POINTER   :: PQ(:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE), POINTER   :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE), POINTER   :: TENDENCY_TMP(:) ! cumulative tendency used as input
    TYPE(STATE_TYPE), POINTER   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB), POINTER:: PVFA(:,:)  ! CC from VDF scheme
    REAL(KIND=JPRB), POINTER:: PVFL(:,:)  ! Liq from VDF scheme
    REAL(KIND=JPRB), POINTER:: PVFI(:,:)  ! Ice from VDF scheme
    REAL(KIND=JPRB), POINTER:: PDYNA(:,:) ! CC from Dynamics
    REAL(KIND=JPRB), POINTER:: PDYNL(:,:) ! Liq from Dynamics
    REAL(KIND=JPRB), POINTER:: PDYNI(:,:) ! Liq from Dynamics
    REAL(KIND=JPRB), POINTER:: PHRSW(:,:) ! Short-wave heating rate
    REAL(KIND=JPRB), POINTER:: PHRLW(:,:) ! Long-wave heating rate
    REAL(KIND=JPRB), POINTER:: PVERVEL(:,:) !Vertical velocity
    REAL(KIND=JPRB), POINTER:: PAP(:,:)   ! Pressure on full levels
    REAL(KIND=JPRB), POINTER:: PAPH(:,:)  ! Pressure on half levels
    REAL(KIND=JPRB), POINTER:: PLSM(:)    ! Land fraction (0-1)
    LOGICAL, POINTER           :: LDCUM(:)   ! Convection active
    INTEGER(KIND=JPIM), POINTER  :: KTYPE(:)   ! Convection type 0,1,2
    REAL(KIND=JPRB), POINTER:: PLU(:,:)   ! Conv. condensate
    REAL(KIND=JPRB), POINTER:: PLUDE(:,:) ! Conv. detrained water
    REAL(KIND=JPRB), POINTER:: PSNDE(:,:) ! Conv. detrained snow
    REAL(KIND=JPRB), POINTER:: PMFU(:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB), POINTER:: PMFD(:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB), POINTER:: PA(:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB), POINTER:: PCLV(:,:,:) 
    REAL(KIND=JPRB), POINTER:: PSUPSAT(:,:)
    REAL(KIND=JPRB), POINTER:: PLCRIT_AER(:,:) 
    REAL(KIND=JPRB), POINTER:: PICRIT_AER(:,:) 
    REAL(KIND=JPRB), POINTER:: PRE_ICE(:,:) 
    REAL(KIND=JPRB), POINTER:: PCCN(:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), POINTER:: PNICE(:,:)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB), POINTER:: PCOVPTOT(:,:) ! Precip fraction
    REAL(KIND=JPRB), POINTER:: PRAINFRAC_TOPRFZ(:) 
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), POINTER   :: PFSQLF(:,:)  ! Flux of liquid
    REAL(KIND=JPRB), POINTER   :: PFSQIF(:,:)  ! Flux of ice
    REAL(KIND=JPRB), POINTER   :: PFCQLNG(:,:) ! -ve corr for liq
    REAL(KIND=JPRB), POINTER   :: PFCQNNG(:,:) ! -ve corr for ice
    REAL(KIND=JPRB), POINTER   :: PFSQRF(:,:)  ! Flux diagnostics
    REAL(KIND=JPRB), POINTER   :: PFSQSF(:,:)  !    for DDH, generic
    REAL(KIND=JPRB), POINTER   :: PFCQRNG(:,:) ! rain
    REAL(KIND=JPRB), POINTER   :: PFCQSNG(:,:) ! snow
    REAL(KIND=JPRB), POINTER   :: PFSQLTUR(:,:) ! liquid flux due to VDF
    REAL(KIND=JPRB), POINTER   :: PFSQITUR(:,:) ! ice flux due to VDF
    REAL(KIND=JPRB), POINTER   :: PFPLSL(:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB), POINTER   :: PFPLSN(:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB), POINTER   :: PFHPSL(:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB), POINTER   :: PFHPSN(:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

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

    FIELD = FSET%FIELD("PEXTRA")
    FSPACE = FIELD%FUNCTIONSPACE()
    NPROMA = FSPACE%BLOCK_SIZE(1)
    NLEV = FSPACE%LEVELS()
    print *, "NPROMA, NUMOMP ", NPROMA, NUMOMP

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID,energy,power) &
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

         !--- a future plan to replace the call to CLOUDSC ------ 
         !
         ! type( block_state_t )
         !   real(c_double), pointer :: PT(:,:)
         !   type(state_type) :: tendency_LOC
         !   type(state_type) :: tendency_TMP
         !   type(state_type) :: tendency_CML
         ! end type
         ! call extract_block( FSET, IBL, config, block_state )
         !     call FSET%FIELD("PT")%BLOCK_DATA(IBL,PT,CONFIG)
         !     call FSET%FIELD("PQ")%BLOCK_DATA(IBL,PQ,CONFIG)
         ! call cloudsc_atlas ( FSET, IBL, config )

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
    
  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
