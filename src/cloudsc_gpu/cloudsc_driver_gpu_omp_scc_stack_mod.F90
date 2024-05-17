! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
MODULE CLOUDSC_DRIVER_GPU_OMP_SCC_STACK_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY: NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY: PERFORMANCE_TIMER, GET_THREAD_NUM
  
  USE CLOUDSC_GPU_OMP_SCC_STACK_MOD, ONLY: CLOUDSC_SCC_STACK
  
  IMPLICIT NONE
  
  CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_STACK (NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, KFLDX, PTSPHY, PT, PQ, &
  & TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
  & BUFFER_CML, BUFFER_TMP, BUFFER_LOC, PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, PHRSW, PHRLW, PVERVEL, PAP, PAPH,  &
  & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER, PICRIT_AER, PRE_ICE, PCCN, PNICE,  &
  & PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG, PFCQSNG, PFSQLTUR, PFSQITUR, PFPLSL,  &
  & PFPLSN, PFHPSL, PFHPSN)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel
    
    USE parkind1, ONLY: JPRB
    USE, intrinsic :: ISO_C_BINDING, ONLY: C_SIZEOF
    INTEGER(KIND=JPIM), INTENT(IN) :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
    REAL(KIND=JPRB), INTENT(IN) :: PTSPHY    ! Physics timestep
    REAL(KIND=JPRB), INTENT(IN) :: PT(NPROMA, NLEV, NGPBLKS)    ! T at start of callpar
    REAL(KIND=JPRB), INTENT(IN) :: PQ(NPROMA, NLEV, NGPBLKS)    ! Q at start of callpar
    TYPE(state_type), INTENT(IN) :: TENDENCY_CML(NGPBLKS)    ! cumulative tendency used for final output
    TYPE(state_type), INTENT(IN) :: TENDENCY_TMP(NGPBLKS)    ! cumulative tendency used as input
    TYPE(state_type), INTENT(OUT) :: TENDENCY_LOC(NGPBLKS)    ! local tendency from cloud scheme
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA, NLEV, 3 + NCLV, NGPBLKS)    ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_TMP(NPROMA, NLEV, 3 + NCLV, NGPBLKS)    ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA, NLEV, 3 + NCLV, NGPBLKS)    ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB), INTENT(IN) :: PVFA(NPROMA, NLEV, NGPBLKS)    ! CC from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PVFL(NPROMA, NLEV, NGPBLKS)    ! Liq from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PVFI(NPROMA, NLEV, NGPBLKS)    ! Ice from VDF scheme
    REAL(KIND=JPRB), INTENT(IN) :: PDYNA(NPROMA, NLEV, NGPBLKS)    ! CC from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PDYNL(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PDYNI(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN) :: PHRSW(NPROMA, NLEV, NGPBLKS)    ! Short-wave heating rate
    REAL(KIND=JPRB), INTENT(IN) :: PHRLW(NPROMA, NLEV, NGPBLKS)    ! Long-wave heating rate
    REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(NPROMA, NLEV, NGPBLKS)    !Vertical velocity
    REAL(KIND=JPRB), INTENT(IN) :: PAP(NPROMA, NLEV, NGPBLKS)    ! Pressure on full levels
    REAL(KIND=JPRB), INTENT(IN) :: PAPH(NPROMA, NLEV + 1, NGPBLKS)    ! Pressure on half levels
    REAL(KIND=JPRB), INTENT(IN) :: PLSM(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    LOGICAL, INTENT(IN) :: LDCUM(NPROMA, NGPBLKS)    ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    REAL(KIND=JPRB), INTENT(IN) :: PLU(NPROMA, NLEV, NGPBLKS)    ! Conv. condensate
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained water
    REAL(KIND=JPRB), INTENT(IN) :: PSNDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained snow
    REAL(KIND=JPRB), INTENT(IN) :: PMFU(NPROMA, NLEV, NGPBLKS)    ! Conv. mass flux up
    REAL(KIND=JPRB), INTENT(IN) :: PMFD(NPROMA, NLEV, NGPBLKS)    ! Conv. mass flux down
    REAL(KIND=JPRB), INTENT(IN) :: PA(NPROMA, NLEV, NGPBLKS)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB), INTENT(IN) :: PCLV(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN) :: PSUPSAT(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN) :: PLCRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN) :: PICRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN) :: PRE_ICE(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN) :: PCCN(NPROMA, NLEV, NGPBLKS)    ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), INTENT(IN) :: PNICE(NPROMA, NLEV, NGPBLKS)    ! ice number concentration (cf. CCN)
    
    REAL(KIND=JPRB), INTENT(INOUT) :: PCOVPTOT(NPROMA, NLEV, NGPBLKS)    ! Precip fraction
    REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(NPROMA, NLEV + 1, NGPBLKS)    ! Flux of liquid
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(NPROMA, NLEV + 1, NGPBLKS)    ! Flux of ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(NPROMA, NLEV + 1, NGPBLKS)    ! -ve corr for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(NPROMA, NLEV + 1, NGPBLKS)    ! -ve corr for ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(NPROMA, NLEV + 1, NGPBLKS)    ! Flux diagnostics
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(NPROMA, NLEV + 1, NGPBLKS)    !    for DDH, generic
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(NPROMA, NLEV + 1, NGPBLKS)    ! rain
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(NPROMA, NLEV + 1, NGPBLKS)    ! snow
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(NPROMA, NLEV + 1, NGPBLKS)    ! liquid flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(NPROMA, NLEV + 1, NGPBLKS)    ! ice flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(NPROMA, NLEV + 1, NGPBLKS)    ! liq+rain sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(NPROMA, NLEV + 1, NGPBLKS)    ! ice+snow sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(NPROMA, NLEV + 1, NGPBLKS)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(NPROMA, NLEV + 1, NGPBLKS)    ! Enthalp flux for ice
    
    INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, JL
    
    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP
    
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID    ! thread id from 0 .. NUMOMP - 1
    INTEGER(KIND=8) :: ISTSZ
    REAL(KIND=JPRB), ALLOCATABLE :: ZSTACK(:, :)
    INTEGER(KIND=8) :: YLSTACK_L
    ISTSZ = (NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8) + 13*NLEV*NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8) +  &
    & 5*NLEV*NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)*NCLV + NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)*NCLV) /  &
    & MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)
    IF (.not.(MOD(NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8) + 13*NLEV*NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8) +  &
    & 5*NLEV*NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)*NCLV + NPROMA*MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)*NCLV,  &
    & MAX(C_SIZEOF(REAL(1, kind=JPRB)), 8)) == 0)) ISTSZ = ISTSZ + 1
    ALLOCATE (ZSTACK(ISTSZ, NGPBLKS))
!$omp target enter data map(alloc: ZSTACK)

    IBL = 1      ! Useless statement to show the compiler that the sepcification part is over!
    
    IF (irank == 0) THEN
1003  FORMAT(5X, 'NUMPROC=', I0, ', NUMOMP=', I0, ', NGPTOTG=', I0, ', NPROMA=', I0, ', NGPBLKS=', I0)
      WRITE(0, 1003) NUMPROC, NUMOMP, NGPTOTG, NPROMA, NGPBLKS
    END IF
    
    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    
    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP
    
!$omp target data &
!$omp map(to: &
!$omp   pt,pq,buffer_cml,buffer_tmp,pvfa, &
!$omp   pvfl,pvfi,pdyna,pdynl,pdyni,phrsw,phrlw,pvervel, &
!$omp   pap,paph,plsm,ldcum,ktype,plu,psnde, &
!$omp   pmfu,pmfd,pa,pclv,psupsat,plcrit_aer,picrit_aer, &
!$omp   pre_ice,pccn,pnice, yrecldp) &
!$omp map(tofrom: &
!$omp   buffer_loc,plude,pcovptot,prainfrac_toprfz) &
!$omp map(from: &
!$omp   pfsqlf,pfsqif,pfcqnng, &
!$omp   pfcqlng ,pfsqrf,pfsqsf,pfcqrng,pfcqsng,pfsqltur, &
!$omp   pfsqitur,pfplsl,pfplsn,pfhpsl,pfhpsn)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)
    
#ifdef HAVE_OMP_TARGET_LOOP_CONSTRUCT
!$omp target teams loop bind(teams) private(YLSTACK_L) thread_limit(nproma)
#else
!$omp target teams distribute private(YLSTACK_L) thread_limit(nproma)
#endif 
    DO JKGLO=1,NGPTOT,NPROMA
      IBL = (JKGLO - 1) / NPROMA + 1
      YLSTACK_L = LOC(ZSTACK(1, IBL))
      ICEND = MIN(NPROMA, NGPTOT - JKGLO + 1)

#ifdef HAVE_OMP_TARGET_LOOP_CONSTRUCT
!$omp loop bind(parallel)
#elif defined(HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD)
!$omp loop bind(thread) private(YLSTACK_L)
#else
!$omp parallel do
#endif
    DO JL=1,ICEND

      CALL CLOUDSC_SCC_STACK(1, ICEND, NPROMA, NLEV, PTSPHY, PT(:, :, IBL), PQ(:, :, IBL), BUFFER_TMP(:, :, 1, IBL),  &
      & BUFFER_TMP(:, :, 3, IBL), BUFFER_TMP(:, :, 2, IBL), BUFFER_TMP(:, :, 4:8, IBL), BUFFER_LOC(:, :, 1, IBL),  &
      & BUFFER_LOC(:, :, 3, IBL), BUFFER_LOC(:, :, 2, IBL), BUFFER_LOC(:, :, 4:8, IBL), PVFA(:, :, IBL), PVFL(:, :, IBL),  &
      & PVFI(:, :, IBL), PDYNA(:, :, IBL), PDYNL(:, :, IBL), PDYNI(:, :, IBL), PHRSW(:, :, IBL), PHRLW(:, :, IBL),  &
      & PVERVEL(:, :, IBL), PAP(:, :, IBL), PAPH(:, :, IBL), PLSM(:, IBL), LDCUM(:, IBL), KTYPE(:, IBL), PLU(:, :, IBL),  &
      & PLUDE(:, :, IBL), PSNDE(:, :, IBL), PMFU(:, :, IBL), PMFD(:, :, IBL), PA(:, :, IBL), PCLV(:, :, :, IBL),  &
      & PSUPSAT(:, :, IBL), PLCRIT_AER(:, :, IBL), PICRIT_AER(:, :, IBL), PRE_ICE(:, :, IBL), PCCN(:, :, IBL), PNICE(:, :, IBL),  &
      & PCOVPTOT(:, :, IBL), PRAINFRAC_TOPRFZ(:, IBL), PFSQLF(:, :, IBL), PFSQIF(:, :, IBL), PFCQNNG(:, :, IBL),  &
      & PFCQLNG(:, :, IBL), PFSQRF(:, :, IBL), PFSQSF(:, :, IBL), PFCQRNG(:, :, IBL), PFCQSNG(:, :, IBL), PFSQLTUR(:, :, IBL),  &
      & PFSQITUR(:, :, IBL), PFPLSL(:, :, IBL), PFPLSN(:, :, IBL), PFHPSL(:, :, IBL), PFHPSN(:, :, IBL), LOCAL_YRECLDP,  &
      & YDSTACK_L=YLSTACK_L, JL=JL)
      !---prognostic fields
      !-- arrays for aerosol-cloud interactions
      !---diagnostic output
      !---resulting fluxes
    
      
    ENDDO  
#ifdef HAVE_OMP_TARGET_LOOP_CONSTRUCT
!$omp end loop
#else
!$omp end parallel do
#endif
     
      
      
    END DO
    
    !-- The "nowait" is here to get correct local timings (tloc) per thread
    !   i.e. we should not wait for slowest thread to finish before measuring tloc
    
    CALL TIMER%THREAD_END(TID)
    
    
!$omp end target data

    CALL TIMER%END()
    
    
    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)
    
    
    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)
 
    DEALLOCATE (ZSTACK)
  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_STACK
  
END MODULE CLOUDSC_DRIVER_GPU_OMP_SCC_STACK_MOD
