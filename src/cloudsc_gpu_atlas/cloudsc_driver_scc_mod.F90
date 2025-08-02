! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_SCC_MOD
  USE PARKIND1, ONLY: JPIM, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_GLOBAL_ATLAS_STATE_MOD, ONLY: OUT_VAR_NAMES
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM, FTIMER
  USE, INTRINSIC :: IEEE_EXCEPTIONS, ONLY : IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_INVALID

#ifdef CLOUDSC_GPU_SCC
  USE CLOUDSC_GPU_SCC_MOD, ONLY: CLOUDSC => CLOUDSC_SCC
#endif
#ifdef CLOUDSC_GPU_SCC_HOIST
  USE CLOUDSC_GPU_SCC_HOIST_MOD, ONLY: CLOUDSC => CLOUDSC_SCC_HOIST
#endif
#ifdef CLOUDSC_GPU_SCC_K_CACHING
  USE CLOUDSC_GPU_SCC_K_CACHING_MOD, ONLY: CLOUDSC => CLOUDSC_SCC_K_CACHING
#endif

  USE ATLAS_MODULE
  USE PLUTO_MODULE, ONLY: PLUTO, PLUTO_ALLOCATOR
  USE, INTRINSIC :: ISO_C_BINDING

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_KERNEL( TIMER, NGPTOT, NGPBLKS, NLEV, NPROMA, PTSPHY, &
      & PLCRIT_AER, PICRIT_AER, PRE_ICE,    PCCN,       PNICE, &
      & PT,         PQ,         PVFA,       PVFL,       PVFI,  & 
      & PDYNA,      PDYNL,      PDYNI,      PHRSW,      PHRLW, &
      & PVERVEL,    PAP,        PLU,        PLUDE,      PSNDE, &
      & PMFU,       PMFD,       PA,         PSUPSAT,           &
      & TENDENCY_CML_T, TENDENCY_CML_A, TENDENCY_CML_Q, &
      & TENDENCY_TMP_T, TENDENCY_TMP_A, TENDENCY_TMP_Q, &
      & PLSM,       LDCUM,      KTYPE,      PAPH,       PCLV,     &
      & TENDENCY_CML_CLD,           TENDENCY_TMP_CLD, &
      & PFSQLF,     PFSQIF,     PFCQNNG,    PFCQLNG,    PFSQRF,   &
      & PFSQSF ,    PFCQRNG,    PFCQSNG,    PFSQLTUR,   PFSQITUR, &
      & PFPLSL,     PFPLSN,     PFHPSL,     PFHPSN,     PCOVPTOT, &
      & TENDENCY_LOC_T, TENDENCY_LOC_A, TENDENCY_LOC_Q, &
      & PRAINFRAC_TOPRFZ, TENDENCY_LOC_CLD )

    TYPE(PERFORMANCE_TIMER), INTENT(INOUT) :: TIMER
    INTEGER(KIND=JPIM), INTENT(IN) :: NGPTOT
    INTEGER(KIND=JPIM), INTENT(IN) :: NGPBLKS
    INTEGER(KIND=JPIM), INTENT(IN) :: NLEV
    INTEGER(KIND=JPIM), INTENT(IN) :: NPROMA
    REAL(KIND=JPRB), INTENT(IN)   :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PLCRIT_AER(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PICRIT_AER(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PRE_ICE(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PT(:,:,:)       ! T at start of callpar
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PQ(:,:,:)       ! Q at start of callpar
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PVFA(:,:,:)     ! CC from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PVFL(:,:,:)     ! Liq from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PVFI(:,:,:)     ! Ice from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PDYNA(:,:,:)    ! CC from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PDYNL(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PDYNI(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PHRSW(:,:,:)    ! Short-wave heating rate
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PHRLW(:,:,:)    ! Long-wave heating rate
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PVERVEL(:,:,:)  ! Vertical velocity
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PAP(:,:,:)      ! Pressure on full levels
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PLU(:,:,:)      ! Conv. condensate
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(INOUT)  :: PLUDE(:,:,:)    ! Conv. detrained water
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PSNDE(:,:,:)    ! Conv. detrained snow
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PMFU(:,:,:)     ! Conv. mass flux up
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PMFD(:,:,:)     ! Conv. mass flux down
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PA(:,:,:)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_CML_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_CML_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_CML_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_TMP_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_TMP_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_TMP_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PLSM(:,:)       ! Land fraction (0-1)
    LOGICAL,         CONTIGUOUS, INTENT(IN)  :: LDCUM(:,:)      ! Convection active
    INTEGER(KIND=JPIM),CONTIGUOUS, INTENT(IN):: KTYPE(:,:)      ! Convection type 0,1,2
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PAPH(:,:,:)     ! Pressure on half levels
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: PCLV(:,:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_CML_CLD(:,:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(IN)  :: TENDENCY_TMP_CLD(:,:,:,:)

    ! output variables
    ! Note: flux diagnostics for DDH budget are named P*
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQLF(:,:,:)    ! Flux of liquid
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQIF(:,:,:)    ! Flux of ice
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFCQLNG(:,:,:)   ! -ve corr for liq
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFCQNNG(:,:,:)   ! -ve corr for ice
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQRF(:,:,:)    ! Flux diagnostics
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQSF(:,:,:)    !    for DDH, generic
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFCQRNG(:,:,:)   ! rain
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFCQSNG(:,:,:)   ! snow
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQLTUR(:,:,:)  ! liquid flux due to VDF
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFSQITUR(:,:,:)  ! ice flux due to VDF
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFPLSL(:,:,:)    ! liq+rain sedim flux
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFPLSN(:,:,:)    ! ice+snow sedim flux
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFHPSL(:,:,:)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PFHPSN(:,:,:)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PCOVPTOT(:,:,:)    ! Precip fraction
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: TENDENCY_LOC_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: TENDENCY_LOC_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: TENDENCY_LOC_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: PRAINFRAC_TOPRFZ(:,:)
    REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT)  :: TENDENCY_LOC_CLD(:,:,:,:)

    INTEGER :: JKGLO, IBL, ICEND, HOIST_POOL_STR_LEN, JL
    CHARACTER(32) :: HOIST_POOL_STR
    LOGICAL :: LL_HALT_INVALID

#ifndef CLOUDSC_GPU_SCC
    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP
#endif

#ifdef CLOUDSC_GPU_SCC_HOIST
    ! Local declarations of promoted temporaries via atlas::pluto
    TYPE(pluto_allocator) :: device_allocator
    REAL(KIND=JPRB), POINTER :: ZFOEALFA(:,:,:), ZTP1(:,:,:), &
        & ZLI(:,:,:), ZA(:,:,:), ZAORIG(:,:,:), ZLIQFRAC(:,:,:), &
        & ZICEFRAC(:,:,:), ZQSMIX(:,:,:), ZQSLIQ(:,:,:), &
        & ZQSICE(:,:,:), ZFOEEWMT(:,:,:), ZFOEEW(:,:,:), ZFOEELIQT(:,:,:)
    REAL(KIND=JPRB), POINTER :: ZQX(:,:,:,:), ZQX0(:,:,:,:), &
        & ZPFPLSX(:,:,:,:), ZLNEG(:,:,:,:), ZQXN2D(:,:,:,:)
    INTEGER(KIND=JPIM) :: HOIST_POOL = 0

    CALL GET_ENVIRONMENT_VARIABLE("HOIST_POOL", VALUE=HOIST_POOL_STR, LENGTH=HOIST_POOL_STR_LEN)
    IF (HOIST_POOL_STR_LEN > 0 ) THEN
      READ(HOIST_POOL_STR(1:HOIST_POOL_STR_LEN), *) HOIST_POOL
    ENDIF
    IF (ALL(HOIST_POOL/=[0,1])) THEN
        PRINT *, "HOIST_POOL must be 0 or 1. Exiting."
        STOP
    END IF
#endif

#ifndef CLOUDSC_GPU_SCC
    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP
#else
!$acc data copyin(YRECLDP)
#endif

#ifdef CLOUDSC_GPU_SCC_HOIST
 if (HOIST_POOL == 1) then
    device_allocator = pluto%make_allocator(pluto%device_pool_resource())
  else
    device_allocator = pluto%make_allocator(pluto%device_resource())
  endif
  call device_allocator%allocate(ZFOEALFA, [NPROMA, NLEV+1, NGPBLKS])
  call device_allocator%allocate(ZTP1, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZLI, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZA, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZAORIG, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZLIQFRAC, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZICEFRAC, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZQX, [NPROMA, NLEV, NCLV, NGPBLKS])
  call device_allocator%allocate(ZQX0, [NPROMA, NLEV, NCLV, NGPBLKS])
  call device_allocator%allocate(ZPFPLSX, [NPROMA, NLEV+1, NCLV, NGPBLKS])
  call device_allocator%allocate(ZLNEG, [NPROMA, NLEV, NCLV, NGPBLKS])
  call device_allocator%allocate(ZQXN2D, [NPROMA, NLEV, NCLV, NGPBLKS])
  call device_allocator%allocate(ZQSMIX, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZQSLIQ, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZQSICE, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZFOEEWMT, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZFOEEW, [NPROMA, NLEV, NGPBLKS])
  call device_allocator%allocate(ZFOEELIQT, [NPROMA, NLEV, NGPBLKS])
#endif

! temporary disable floating-point-trapping of FE_INVALID caused by
! too aggressive optimisation within the CLOUDSC function
  CALL IEEE_GET_HALTING_MODE(IEEE_INVALID, LL_HALT_INVALID)
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)

!$acc data deviceptr(&
!$acc & PLCRIT_AER, PICRIT_AER, PRE_ICE,    PCCN,       PNICE, &
!$acc & PT,         PQ,         PVFA,       PVFL,       PVFI,  &
!$acc & PDYNA,      PDYNL,      PDYNI,      PHRSW,      PHRLW, &
!$acc & PVERVEL,    PAP,        PLU,        PLUDE,      PSNDE, &
!$acc & PMFU,       PMFD,       PA,         PSUPSAT,           &
!$acc & TENDENCY_CML_T, TENDENCY_CML_A, TENDENCY_CML_Q, &
!$acc & TENDENCY_TMP_T, TENDENCY_TMP_A, TENDENCY_TMP_Q, &
!$acc & PLSM,       LDCUM,      KTYPE,      PAPH,       PCLV,     &
!$acc & TENDENCY_CML_CLD,           TENDENCY_TMP_CLD, &
!$acc & PFSQLF,     PFSQIF,     PFCQNNG,    PFCQLNG,    PFSQRF,   &
!$acc & PFSQSF ,    PFCQRNG,    PFCQSNG,    PFSQLTUR,   PFSQITUR, &
!$acc & PFPLSL,     PFPLSN,     PFHPSL,     PFHPSN,     PCOVPTOT, &
!$acc & TENDENCY_LOC_T, TENDENCY_LOC_A, TENDENCY_LOC_Q, &
!$acc & PRAINFRAC_TOPRFZ, TENDENCY_LOC_CLD )

    ! Local timer for each thread. The main code part does not run threaded, so TID=0
    CALL TIMER%THREAD_START(TID=0)

!$acc parallel loop gang vector_length(NPROMA)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

#ifndef CLOUDSC_GPU_SCC
!$acc loop vector
      DO JL=1,ICEND
#endif
        CALL CLOUDSC &
        & (1, ICEND, NPROMA, NLEV, PTSPHY,&
        & PT(:,:,IBL), PQ(:,:,IBL), &
        & TENDENCY_TMP_T(:,:,IBL), TENDENCY_TMP_Q(:,:,IBL), TENDENCY_TMP_A(:,:,IBL), TENDENCY_TMP_CLD(:,:,:,IBL), &
        & TENDENCY_LOC_T(:,:,IBL), TENDENCY_LOC_Q(:,:,IBL), TENDENCY_LOC_A(:,:,IBL), TENDENCY_LOC_CLD(:,:,:,IBL), &
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
#if defined CLOUDSC_GPU_SCC
        & YRECLDP=YRECLDP)
#else
        & LOCAL_YRECLDP, &
#if defined CLOUDSC_GPU_SCC_HOIST
        & ZFOEALFA(:,:,IBL), ZTP1(:,:,IBL), ZLI(:,:,IBL), ZA(:,:,IBL), ZAORIG(:,:,IBL), &
        & ZLIQFRAC(:,:,IBL), ZICEFRAC(:,:,IBL), ZQX(:,:,:,IBL), ZQX0(:,:,:,IBL), ZPFPLSX(:,:,:,IBL), &
        & ZLNEG(:,:,:,IBL), ZQXN2D(:,:,:,IBL), ZQSMIX(:,:,IBL), ZQSLIQ(:,:,IBL), ZQSICE(:,:,IBL), &
        & ZFOEEWMT(:,:,IBL), ZFOEEW(:,:,IBL), ZFOEELIQT(:,:,IBL), &
#endif
        JL=JL)
#endif
#ifndef CLOUDSC_GPU_SCC
      ENDDO
#endif
    ENDDO
!$acc end parallel loop
!$acc end data

    ! The main code part does not run threaded, so TID=0
    CALL TIMER%THREAD_END(TID=0)

    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, LL_HALT_INVALID)

#ifdef CLOUDSC_GPU_SCC
!$acc end data
#endif

#ifdef CLOUDSC_GPU_SCC_HOIST
    call device_allocator%deallocate(ZFOEALFA)
    call device_allocator%deallocate(ZTP1)
    call device_allocator%deallocate(ZLI)
    call device_allocator%deallocate(ZA)
    call device_allocator%deallocate(ZAORIG)
    call device_allocator%deallocate(ZLIQFRAC)
    call device_allocator%deallocate(ZICEFRAC)
    call device_allocator%deallocate(ZQX)
    call device_allocator%deallocate(ZQX0)
    call device_allocator%deallocate(ZPFPLSX)
    call device_allocator%deallocate(ZLNEG)
    call device_allocator%deallocate(ZQXN2D)
    call device_allocator%deallocate(ZQSMIX)
    call device_allocator%deallocate(ZQSLIQ)
    call device_allocator%deallocate(ZQSICE)
    call device_allocator%deallocate(ZFOEEWMT)
    call device_allocator%deallocate(ZFOEEW)
    call device_allocator%deallocate(ZFOEELIQT)
#endif
  END SUBROUTINE CLOUDSC_KERNEL


  SUBROUTINE CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOTG, PTSPHY)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NGPTOTG
    REAL(KIND=JPRB), INTENT(IN)   :: PTSPHY       ! Physics timestep

    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    TYPE(ATLAS_FIELD) :: FIELD
    INTEGER(KIND=JPIM)    :: NPROMA, NLEV, NGPTOT

    INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    TYPE(ATLAS_TRACE)  :: TRACE

    ! input variables
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLCRIT_AER(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PICRIT_AER(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PRE_ICE(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)
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
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLU(:,:,:)      ! Conv. condensate
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLUDE(:,:,:)    ! Conv. detrained water
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PSNDE(:,:,:)    ! Conv. detrained snow
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PMFU(:,:,:)     ! Conv. mass flux up
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PMFD(:,:,:)     ! Conv. mass flux down
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PA(:,:,:)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_CML_T(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_CML_A(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_CML_Q(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_TMP_T(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_TMP_A(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_TMP_Q(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PLSM(:,:)       ! Land fraction (0-1)
    LOGICAL, POINTER, CONTIGUOUS          :: LDCUM(:,:)      ! Convection active
    INTEGER(KIND=JPIM),POINTER, CONTIGUOUS:: KTYPE(:,:)      ! Convection type 0,1,2
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PAPH(:,:,:)     ! Pressure on half levels
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCLV(:,:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_CML_CLD(:,:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_TMP_CLD(:,:,:,:)

    ! output variables
    ! Note: flux diagnostics for DDH budget are named P*
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
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PCOVPTOT(:,:,:)    ! Precip fraction
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_LOC_T(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_LOC_A(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_LOC_Q(:,:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: PRAINFRAC_TOPRFZ(:,:)
    REAL(KIND=JPRB), POINTER, CONTIGUOUS  :: TENDENCY_LOC_CLD(:,:,:,:)

    REAL(KIND=JPRB), POINTER :: TMP3D(:,:,:,:)
    REAL(KIND=JPRD) :: TDIFF

    TRACE = ATLAS_TRACE("cloudsc_driver_mod.F90", __LINE__, "CLOUDSC_DRIVER","COMPUTE")

    FIELD = FSET%FIELD("PCLV")
    FSPACE = FIELD%FUNCTIONSPACE()
    NPROMA = FIELD%SHAPE(1)
    NLEV = FSPACE%LEVELS()
    NGPTOT = FSPACE%SIZE()
    NGPBLKS = FIELD%SHAPE(4)

1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGTOT=', i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOT, NGPTOTG,NPROMA,NGPBLKS
    end if

    ! timing of the device allocate
    TDIFF = FTIMER()
    CALL FSET%ALLOCATE_DEVICE()
    TDIFF = FTIMER() - TDIFF
    PRINT *, "Device allocate [ms]: ", INT(TDIFF * 1000.0_JPRD)

    CALL TIMER%START(NUMOMP)

    CALL FSET%UPDATE_DEVICE([(JKGLO, JKGLO=1,37)])

    ! input
    CALL FSET%DEVICE_DATA(1, PLCRIT_AER)
    CALL FSET%DEVICE_DATA(2, PICRIT_AER)
    CALL FSET%DEVICE_DATA(3, PRE_ICE)
    CALL FSET%DEVICE_DATA(4, PCCN)
    CALL FSET%DEVICE_DATA(5, PNICE)
    CALL FSET%DEVICE_DATA(6, PT)
    CALL FSET%DEVICE_DATA(7, PQ)
    CALL FSET%DEVICE_DATA(8, PVFA)
    CALL FSET%DEVICE_DATA(9, PVFL)
    CALL FSET%DEVICE_DATA(10, PVFI)
    CALL FSET%DEVICE_DATA(11, PDYNA)
    CALL FSET%DEVICE_DATA(12, PDYNL)
    CALL FSET%DEVICE_DATA(13, PDYNI)
    CALL FSET%DEVICE_DATA(14, PHRSW)
    CALL FSET%DEVICE_DATA(15, PHRLW)
    CALL FSET%DEVICE_DATA(16, PVERVEL)
    CALL FSET%DEVICE_DATA(17, PAP)
    CALL FSET%DEVICE_DATA(18, PLU)
    CALL FSET%DEVICE_DATA(19, PLUDE) ! in/out
    CALL FSET%DEVICE_DATA(20, PSNDE)
    CALL FSET%DEVICE_DATA(21, PMFU)
    CALL FSET%DEVICE_DATA(22, PMFD)
    CALL FSET%DEVICE_DATA(23, PA)
    CALL FSET%DEVICE_DATA(24, PSUPSAT)
    CALL FSET%DEVICE_DATA(25, TENDENCY_CML_T)
    CALL FSET%DEVICE_DATA(26, TENDENCY_CML_A)
    CALL FSET%DEVICE_DATA(27, TENDENCY_CML_Q)
    CALL FSET%DEVICE_DATA(28, TENDENCY_TMP_T)
    CALL FSET%DEVICE_DATA(29, TENDENCY_TMP_A)
    CALL FSET%DEVICE_DATA(30, TENDENCY_TMP_Q)
    CALL FSET%DEVICE_DATA(31, PLSM)
    CALL FSET%DEVICE_DATA(32, LDCUM)
    CALL FSET%DEVICE_DATA(33, KTYPE)
    CALL FSET%DEVICE_DATA(34, PAPH)
    CALL FSET%DEVICE_DATA(35, PCLV)
    CALL FSET%DEVICE_DATA(36, TENDENCY_CML_CLD)
    CALL FSET%DEVICE_DATA(37, TENDENCY_TMP_CLD)

    ! output
    CALL FSET%DEVICE_DATA(38, PFSQLF)
    CALL FSET%DEVICE_DATA(39, PFSQIF)
    CALL FSET%DEVICE_DATA(40, PFCQLNG)
    CALL FSET%DEVICE_DATA(41, PFCQNNG)
    CALL FSET%DEVICE_DATA(42, PFSQRF)
    CALL FSET%DEVICE_DATA(43, PFSQSF)
    CALL FSET%DEVICE_DATA(44, PFCQRNG)
    CALL FSET%DEVICE_DATA(45, PFCQSNG)
    CALL FSET%DEVICE_DATA(46, PFSQLTUR)
    CALL FSET%DEVICE_DATA(47, PFSQITUR)
    CALL FSET%DEVICE_DATA(48, PFPLSL)
    CALL FSET%DEVICE_DATA(49, PFPLSN)
    CALL FSET%DEVICE_DATA(50, PFHPSL)
    CALL FSET%DEVICE_DATA(51, PFHPSN)
    CALL FSET%DEVICE_DATA(52, PCOVPTOT)
    CALL FSET%DEVICE_DATA(53, PRAINFRAC_TOPRFZ)
    CALL FSET%DEVICE_DATA(54, TENDENCY_LOC_T)
    CALL FSET%DEVICE_DATA(55, TENDENCY_LOC_A)
    CALL FSET%DEVICE_DATA(56, TENDENCY_LOC_Q)
    CALL FSET%DEVICE_DATA(57, TENDENCY_LOC_CLD)

    CALL CLOUDSC_KERNEL( TIMER, NGPTOT, NGPBLKS, NLEV, NPROMA, PTSPHY, &
      & PLCRIT_AER, PICRIT_AER, PRE_ICE,    PCCN,       PNICE, &
      & PT,         PQ,         PVFA,       PVFL,       PVFI,  & 
      & PDYNA,      PDYNL,      PDYNI,      PHRSW,      PHRLW, &
      & PVERVEL,    PAP,        PLU,        PLUDE,      PSNDE, &
      & PMFU,       PMFD,       PA,         PSUPSAT,           &
      & TENDENCY_CML_T, TENDENCY_CML_A, TENDENCY_CML_Q, &
      & TENDENCY_TMP_T, TENDENCY_TMP_A, TENDENCY_TMP_Q, &
      & PLSM,       LDCUM,      KTYPE,      PAPH,       PCLV,     &
      & TENDENCY_CML_CLD,           TENDENCY_TMP_CLD, &
      & PFSQLF,     PFSQIF,     PFCQNNG,    PFCQLNG,    PFSQRF,   &
      & PFSQSF ,    PFCQRNG,    PFCQSNG,    PFSQLTUR,   PFSQITUR, &
      & PFPLSL,     PFPLSN,     PFHPSL,     PFHPSN,     PCOVPTOT, &
      & TENDENCY_LOC_T, TENDENCY_LOC_A, TENDENCY_LOC_Q, &
      & PRAINFRAC_TOPRFZ, TENDENCY_LOC_CLD &
    )

    CALL TIMER%END()

    CALL FSET%UPDATE_HOST(["PLUDE"])
    CALL FSET%UPDATE_HOST(OUT_VAR_NAMES)

    CALL FSET%DEALLOCATE_DEVICE()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=0, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

    CALL FIELD%FINAL()
    CALL FSPACE%FINAL()
    CALL TRACE%FINAL()

  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_SCC_MOD
