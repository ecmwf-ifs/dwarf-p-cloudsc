MODULE CLOUDSC_DRIVER_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE


CONTAINS
  SUBROUTINE CLOUDSC_DRIVER( &
     & NPROMA, NLEV, NGPTOT, KFLDX, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
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
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    INTEGER(KIND=JPIM)                                    :: NPROMA, NLEV, NGPTOT
    INTEGER(KIND=JPIM)                                    :: KFLDX 
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_TMP(:) ! cumulative tendency used as input
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
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
              & PT(:,:,IBL), PQ(:,:,IBL), TENDENCY_CML(IBL), TENDENCY_TMP(IBL), TENDENCY_LOC(IBL), &
              & PVFA(:,:,IBL), PVFL(:,:,IBL), PVFI(:,:,IBL), PDYNA(:,:,IBL), PDYNL(:,:,IBL), PDYNI(:,:,IBL), &
              & PHRSW(:,:,IBL),    PHRLW(:,:,IBL),&
              & PVERVEL(:,:,IBL),  PAP(:,:,IBL),      PAPH(:,:,IBL),&
              & PLSM(:,IBL),       LDCUM(:,IBL),      KTYPE(:,IBL), &
              & PLU(:,:,IBL),      PLUDE(:,:,IBL),    PSNDE(:,:,IBL),    PMFU(:,:,IBL),     PMFD(:,:,IBL),&
              & LDSLPHY,  LDMAINCALL, &
              !---prognostic fields
              & PA(:,:,IBL),&
              & PCLV(:,:,:,IBL),  &
              & PSUPSAT(:,:,IBL),&
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
              & PEXTRA(:,:,:,IBL),   KFLDX)

      END DO
    
  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
