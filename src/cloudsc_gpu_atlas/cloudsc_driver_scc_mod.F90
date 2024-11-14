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
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM, FTIMER

#ifdef CLOUDSC_GPU_SCC
  USE CLOUDSC_GPU_SCC_MOD, ONLY: CLOUDSC_SCC
#endif
#ifdef CLOUDSC_GPU_SCC_HOIST
  USE CLOUDSC_GPU_SCC_HOIST_MOD, ONLY: CLOUDSC_SCC_HOIST
#endif
#ifdef CLOUDSC_GPU_SCC_K_CACHING
  USE CLOUDSC_GPU_SCC_K_CACHING_MOD, ONLY: CLOUDSC_SCC_K_CACHING
#endif

  USE ATLAS_MODULE
  USE, INTRINSIC :: ISO_C_BINDING
  USE ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS_MODULE

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
    REAL(KIND=JPRB), CONTIGUOUS  :: PLCRIT_AER(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: PICRIT_AER(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: PRE_ICE(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), CONTIGUOUS  :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB), CONTIGUOUS  :: PT(:,:,:)       ! T at start of callpar
    REAL(KIND=JPRB), CONTIGUOUS  :: PQ(:,:,:)       ! Q at start of callpar
    REAL(KIND=JPRB), CONTIGUOUS  :: PVFA(:,:,:)     ! CC from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS  :: PVFL(:,:,:)     ! Liq from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS  :: PVFI(:,:,:)     ! Ice from VDF scheme
    REAL(KIND=JPRB), CONTIGUOUS  :: PDYNA(:,:,:)    ! CC from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS  :: PDYNL(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS  :: PDYNI(:,:,:)    ! Liq from Dynamics
    REAL(KIND=JPRB), CONTIGUOUS  :: PHRSW(:,:,:)    ! Short-wave heating rate
    REAL(KIND=JPRB), CONTIGUOUS  :: PHRLW(:,:,:)    ! Long-wave heating rate
    REAL(KIND=JPRB), CONTIGUOUS  :: PVERVEL(:,:,:)  ! Vertical velocity
    REAL(KIND=JPRB), CONTIGUOUS  :: PAP(:,:,:)      ! Pressure on full levels
    REAL(KIND=JPRB), CONTIGUOUS  :: PLU(:,:,:)      ! Conv. condensate
    REAL(KIND=JPRB), CONTIGUOUS  :: PLUDE(:,:,:)    ! Conv. detrained water
    REAL(KIND=JPRB), CONTIGUOUS  :: PSNDE(:,:,:)    ! Conv. detrained snow
    REAL(KIND=JPRB), CONTIGUOUS  :: PMFU(:,:,:)     ! Conv. mass flux up
    REAL(KIND=JPRB), CONTIGUOUS  :: PMFD(:,:,:)     ! Conv. mass flux down
    REAL(KIND=JPRB), CONTIGUOUS  :: PA(:,:,:)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), CONTIGUOUS  :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_CML_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_CML_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_CML_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_TMP_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_TMP_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_TMP_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: PLSM(:,:)       ! Land fraction (0-1)
    LOGICAL, CONTIGUOUS          :: LDCUM(:,:)      ! Convection active
    INTEGER(KIND=JPIM),CONTIGUOUS:: KTYPE(:,:)      ! Convection type 0,1,2
    REAL(KIND=JPRB), CONTIGUOUS  :: PAPH(:,:,:)     ! Pressure on half levels
    REAL(KIND=JPRB), CONTIGUOUS  :: PCLV(:,:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_CML_CLD(:,:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_TMP_CLD(:,:,:,:)

    ! output variables
    ! Note: flux diagnostics for DDH budget are named P*
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQLF(:,:,:)    ! Flux of liquid
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQIF(:,:,:)    ! Flux of ice
    REAL(KIND=JPRB), CONTIGUOUS  :: PFCQLNG(:,:,:)   ! -ve corr for liq
    REAL(KIND=JPRB), CONTIGUOUS  :: PFCQNNG(:,:,:)   ! -ve corr for ice
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQRF(:,:,:)    ! Flux diagnostics
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQSF(:,:,:)    !    for DDH, generic
    REAL(KIND=JPRB), CONTIGUOUS  :: PFCQRNG(:,:,:)   ! rain
    REAL(KIND=JPRB), CONTIGUOUS  :: PFCQSNG(:,:,:)   ! snow
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQLTUR(:,:,:)  ! liquid flux due to VDF
    REAL(KIND=JPRB), CONTIGUOUS  :: PFSQITUR(:,:,:)  ! ice flux due to VDF
    REAL(KIND=JPRB), CONTIGUOUS  :: PFPLSL(:,:,:)    ! liq+rain sedim flux
    REAL(KIND=JPRB), CONTIGUOUS  :: PFPLSN(:,:,:)    ! ice+snow sedim flux
    REAL(KIND=JPRB), CONTIGUOUS  :: PFHPSL(:,:,:)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), CONTIGUOUS  :: PFHPSN(:,:,:)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB), CONTIGUOUS  :: PCOVPTOT(:,:,:)    ! Precip fraction
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_LOC_T(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_LOC_A(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_LOC_Q(:,:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: PRAINFRAC_TOPRFZ(:,:)
    REAL(KIND=JPRB), CONTIGUOUS  :: TENDENCY_LOC_CLD(:,:,:,:)

    INTEGER :: JKGLO, IBL, ICEND
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

#ifdef CLOUDSC_GPU_SCC_HOIST
    ! Local declarations of promoted temporaries
    REAL(KIND=JPRB) :: ZFOEALFA(NPROMA, NLEV+1, NGPBLKS)
    REAL(KIND=JPRB) :: ZTP1(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZLI(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZA(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZAORIG(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZLIQFRAC(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZICEFRAC(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQX(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQX0(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: ZPFPLSX(NPROMA, NLEV+1, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: ZLNEG(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQXN2D(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQSMIX(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQSLIQ(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZQSICE(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZFOEEWMT(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZFOEEW(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: ZFOEELIQT(NPROMA, NLEV, NGPBLKS)
#endif

#if (defined CLOUDSC_GPU_SCC_HOIST) || (defined CLOUDSC_GPU_SCC_K_CACHING)
    INTEGER(KIND=JPIM) :: JL
    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP
#else
!$acc data copyin(YRECLDP)
#endif

#ifdef CLOUDSC_GPU_SCC_HOIST
!$acc enter data create(ZFOEALFA, ZTP1, ZLI, ZA, ZAORIG, ZLIQFRAC, ZICEFRAC, ZQX, ZQX0,  &
!$acc &   ZPFPLSX, ZLNEG, ZQXN2D, ZQSMIX, ZQSLIQ, ZQSICE, ZFOEEWMT,  &
!$acc &   ZFOEEW, ZFOEELIQT)
#endif

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

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

!$acc parallel loop gang vector_length(NPROMA)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

#ifdef CLOUDSC_GPU_SCC
       CALL CLOUDSC_SCC &
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
        & YRECLDP=YRECLDP)
#elif CLOUDSC_GPU_SCC_HOIST
!$acc loop vector
      DO JL=1,ICEND
        CALL CLOUDSC_SCC_HOIST &
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
        & LOCAL_YRECLDP, &
        & ZFOEALFA(:,:,IBL), ZTP1(:,:,IBL), ZLI(:,:,IBL), ZA(:,:,IBL), ZAORIG(:,:,IBL), &
        & ZLIQFRAC(:,:,IBL), ZICEFRAC(:,:,IBL), ZQX(:,:,:,IBL), ZQX0(:,:,:,IBL), ZPFPLSX(:,:,:,IBL), &
        & ZLNEG(:,:,:,IBL), ZQXN2D(:,:,:,IBL), ZQSMIX(:,:,IBL), ZQSLIQ(:,:,IBL), ZQSICE(:,:,IBL), &
        & ZFOEEWMT(:,:,IBL), ZFOEEW(:,:,IBL), ZFOEELIQT(:,:,IBL), JL=JL)
      ENDDO
#elif CLOUDSC_GPU_SCC_K_CACHING
!$acc loop vector
      DO JL=1,ICEND
       CALL CLOUDSC_SCC_K_CACHING &
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
        & YRECLDP=LOCAL_YRECLDP, JL=JL)
      ENDDO
#endif

  ENDDO
!$acc end parallel loop

    CALL TIMER%THREAD_END(TID)
!$acc end data
#ifdef CLOUDSC_GPU_SCC
!$acc end data
#endif
  END SUBROUTINE CLOUDSC_KERNEL


  SUBROUTINE CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOTG, KFLDX, PTSPHY)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NGPTOTG, KFLDX
    REAL(KIND=JPRB), INTENT(IN)   :: PTSPHY       ! Physics timestep

    !TYPE(CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW) :: FBLOCK
    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    TYPE(ATLAS_FIELD) :: FIELD
    INTEGER(KIND=JPIM)    :: NPROMA, NLEV, NGPTOT

    INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER, TIMER_DEVICE_ALLOC
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
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

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
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
    CALL FSET%DEVICE_DATA(19, PLUDE) !in/out
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
    CALL FSET%DEVICE_DATA(53, TENDENCY_LOC_T)
    CALL FSET%DEVICE_DATA(54, TENDENCY_LOC_A)
    CALL FSET%DEVICE_DATA(55, TENDENCY_LOC_Q)
    CALL FSET%DEVICE_DATA(56, PRAINFRAC_TOPRFZ)
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

    CALL FSET%UPDATE_HOST([character(64) :: "PLUDE", "PCOVPTOT", "PRAINFRAC_TOPRFZ", "TENDENCY_LOC_T", &
      & "TENDENCY_LOC_A", "TENDENCY_LOC_Q", "TENDENCY_LOC_CLD", "PFSQLF", "PFSQIF", "PFCQLNG", "PFCQNNG", &
      & "PFSQRF", "PFSQSF", "PFCQRNG", "PFCQSNG", "PFSQLTUR", "PFSQITUR", "PFPLSL", "PFPLSN", "PFHPSL", "PFHPSN"])
    CALL FSET%DEALLOCATE_DEVICE()

    CALL TIMER%END()
    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

    CALL FIELD%FINAL()
    CALL FSPACE%FINAL()
    CALL TRACE%FINAL()

  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_SCC_MOD
