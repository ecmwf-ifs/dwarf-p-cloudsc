MODULE CLOUDSC_DRIVER_GPU_NAIVE_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE CLOUDSC_ACC_KERNELS_MOD, ONLY: CLOUDSC_ACC_KERNELS

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_NAIVE( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, KFLDX, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
     & BUFFER_CML, BUFFER_TMP, BUFFER_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & LDSLPHY,  LDMAINCALL, PA, &
     & PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
     & PEXTRA    )
    ! Driver routine that invokes the naive CLOUDSC GPU kernel

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    INTEGER(KIND=JPIM)                                    :: KFLDX 
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_TMP(:) ! cumulative tendency used as input
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: BUFFER_CML(:,:,:,:) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: BUFFER_TMP(:,:,:,:) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: BUFFER_LOC(:,:,:,:) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFA(:,:,:)  ! CC from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFL(:,:,:)  ! Liq from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFI(:,:,:)  ! Ice from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNA(:,:,:) ! CC from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNL(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNI(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PHRSW(:,:,:) ! Short-wave heating rate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PHRLW(:,:,:) ! Long-wave heating rate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVERVEL(:,:,:) !Vertical velocity
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLSM(:,:)    ! Land fraction (0-1) 
    LOGICAL           ,POINTER, CONTIGUOUS, INTENT(IN)    :: LDCUM(:,:)   ! Convection active
    INTEGER(KIND=JPIM),POINTER, CONTIGUOUS, INTENT(IN)    :: KTYPE(:,:)   ! Convection type 0,1,2
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PSNDE(:,:,:) ! Conv. detrained snow
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    LOGICAL                                               :: LDSLPHY 
    LOGICAL                                               :: LDMAINCALL   ! T if main call to cloudsc
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: PEXTRA(:,:,:,:) ! extra fields
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLCRIT_AER(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PICRIT_AER(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PRE_ICE(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PRAINFRAC_TOPRFZ(:,:) 
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQLF(:,:,:)  ! Flux of liquid
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQIF(:,:,:)  ! Flux of ice
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQLNG(:,:,:) ! -ve corr for liq
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQNNG(:,:,:) ! -ve corr for ice
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQRF(:,:,:)  ! Flux diagnostics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQSF(:,:,:)  !    for DDH, generic
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQRNG(:,:,:) ! rain
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQSNG(:,:,:) ! snow
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQLTUR(:,:,:) ! liquid flux due to VDF
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQITUR(:,:,:) ! ice flux due to VDF
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

!$acc data &
!$acc copyin( &
!$acc   pt,pq,buffer_cml,buffer_tmp,pvfa, &
!$acc   pvfl,pvfi,pdyna,pdynl,pdyni,phrsw,phrlw,pvervel, &
!$acc   pap,paph,plsm,ldcum,ktype,plu,psnde, &
!$acc   pmfu,pmfd,pa,pclv,psupsat,plcrit_aer,picrit_aer, &
!$acc   pre_ice,pccn,pnice, yrecldp) &
!$acc copy( &
!$acc   buffer_loc,plude,pcovptot,prainfrac_toprfz) &
!$acc copyout( &
!$acc   pfsqlf,pfsqif,pfcqnng, &
!$acc   pfcqlng ,pfsqrf,pfsqsf,pfcqrng,pfcqsng,pfsqltur, &
!$acc   pfsqitur,pfplsl,pfplsn,pfhpsl,pfhpsn)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

       !-- These were uninitialized : meaningful only when we compare error differences
       PCOVPTOT(:,:,IBL) = 0.0_JPRB
       BUFFER_LOC(:,:,6+NCLV,IBL) = 0.0_JPRB

       CALL CLOUDSC_ACC_KERNELS &
        & (    1,    ICEND,    NPROMA,  NLEV,&
        & PTSPHY,&
        & PT(:,:,IBL), PQ(:,:,IBL), &
        & BUFFER_CML(:,:,3,IBL), BUFFER_CML(:,:,6,IBL), BUFFER_CML(:,:,5,IBL), BUFFER_CML(:,:,7:11,IBL), &
        & BUFFER_TMP(:,:,3,IBL), BUFFER_TMP(:,:,6,IBL), BUFFER_TMP(:,:,5,IBL), BUFFER_TMP(:,:,7:11,IBL), &
        & BUFFER_LOC(:,:,3,IBL), BUFFER_LOC(:,:,6,IBL), BUFFER_LOC(:,:,5,IBL), BUFFER_LOC(:,:,7:11,IBL), &
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
        & KFLDX, YRECLDP)

       ! Log number of columns processed by this thread
       CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
    ENDDO

    CALL TIMER%THREAD_END(TID)

!$acc end data

    CALL TIMER%END()

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_NAIVE

END MODULE CLOUDSC_DRIVER_GPU_NAIVE_MOD
