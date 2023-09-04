! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
MODULE CLOUDSC_DRIVER_MOD
  ! USE cudafor
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY: NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY: PERFORMANCE_TIMER, GET_THREAD_NUM

  !!!
  USE YOMCST, ONLY: RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV
  USE YOETHF, ONLY: R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE,  &
    & RTICECU, RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
  USE YOECLDP, ONLY: TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  !!!

  !USE YOMCST_CUF, ONLY: YOMCST_UPDATE_DEVICE
  !USE YOETHF_CUF, ONLY: YOETHF_UPDATE_DEVICE
  
  ! USE CLOUDSC_GPU_OMP_SCC_K_CACHING_DR_LOOP_MOD, ONLY: CLOUDSC_SCC_CUF_K_CACHING
  USE cloudsc_c_k_caching_hip_mod, only: cloudsc_c_hip  
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE CLOUDSC_DRIVER (NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, KFLDX, PTSPHY, PT, PQ, TENDENCY_CML, TENDENCY_TMP,  &
  & TENDENCY_LOC, BUFFER_CML, BUFFER_TMP, BUFFER_LOC, PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, PHRSW, PHRLW, PVERVEL, PAP, PAPH,  &
  & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER, PICRIT_AER, PRE_ICE, PCCN, PNICE,  &
  & PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG, PFCQSNG, PFSQLTUR, PFSQITUR, PFPLSL,  &
  & PFPLSN, PFHPSL, PFHPSN)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel
    
    
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
   
    REAL(KIND=JPRB) :: TENDENCY_TMP_T(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_TMP_Q(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_TMP_A(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_TMP_CLD(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_LOC_T(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_LOC_Q(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_LOC_A(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB) :: TENDENCY_LOC_CLD(NPROMA, NLEV, NCLV, NGPBLKS)

    ! tend_loc_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_loc_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_loc_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_loc_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );
    ! tend_cml_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_cml_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_cml_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_cml_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );
    ! tend_tmp_t   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_tmp_q   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_tmp_a   = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
    ! tend_tmp_cld = (double*) malloc( sizeof(double) * nblocks*nlev*nproma*nclv );

    INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND
    
    TYPE(TECLDP) :: LOCAL_YRECLDP

    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID    ! thread id from 0 .. NUMOMP - 1
    !TYPE(TECLDP), DEVICE :: YRECLDP_d
    INTEGER :: istat
    
    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    LOCAL_YRECLDP = YRECLDP

    TENDENCY_TMP_T = BUFFER_TMP(:,:,1,:) 
    TENDENCY_TMP_Q = BUFFER_TMP(:,:,3,:)
    TENDENCY_TMP_A = BUFFER_TMP(:,:,2,:)
    TENDENCY_TMP_CLD = BUFFER_TMP(:,:,4:8,:)
    TENDENCY_LOC_T = BUFFER_LOC(:,:,1,:)
    TENDENCY_LOC_Q = BUFFER_LOC(:,:,3,:)
    TENDENCY_LOC_A = BUFFER_LOC(:,:,2,:)
    TENDENCY_LOC_CLD = BUFFER_LOC(:,:,4:8,:)

    !$omp target data &
    !$omp map(to: &
    !!! $omp   pt,pq,buffer_tmp,pvfa, &
    !$omp   pt,pq,tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld,pvfa, &
    !$omp   pvfl,pvfi,pdyna,pdynl,pdyni,phrsw,phrlw,pvervel, &
    !$omp   pap,paph,plsm,ldcum,ktype,plu,psnde, &
    !$omp   pmfu,pmfd,pa,pclv,psupsat,plcrit_aer,picrit_aer, &
    !$omp   pre_ice,pccn,pnice, yrecldp) &
    !$omp map(tofrom: &
    !!! $omp   buffer_loc,plude,pcovptot,prainfrac_toprfz) &
    !$omp   tendency_loc_t, tendency_loc_q, tendency_loc_a, tendency_loc_cld, plude,pcovptot,prainfrac_toprfz) &
    !$omp map(from: &
    !$omp   pfsqlf,pfsqif,pfcqnng, &
    !$omp   pfcqlng ,pfsqrf,pfsqsf,pfcqrng,pfcqsng,pfsqltur, &
    !$omp   pfsqitur,pfplsl,pfplsn,pfhpsl,pfhpsn)

    
    IF (irank == 0) THEN
 1003 FORMAT(5X, 'NUMPROC=', I0, ', NUMOMP=', I0, ', NGPTOTG=', I0, ', NPROMA=', I0, ', NGPBLKS=', I0)
      WRITE(0, 1003) NUMPROC, NUMOMP, NGPTOTG, NPROMA, NGPBLKS
    END IF
    
    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)
    
    !BLOCKDIM = DIM3(NPROMA, 1, 1)
    !GRIDDIM = DIM3(1, 1, CEILING(REAL(NGPTOT) / REAL(NPROMA)))
    IBL = (JKGLO - 1) / NPROMA + 1
    ICEND = MIN(NPROMA, NGPTOT - JKGLO + 1)
    
    ! $omp target


    ! BUFFER_TMP(:,:,1,IBL), BUFFER_TMP(:,:,3,IBL), BUFFER_TMP(:,:,2,IBL), BUFFER_TMP(:,:,4:8,IBL), &
    !    & BUFFER_LOC(:,:,1,IBL), BUFFER_LOC(:,:,3,IBL), BUFFER_LOC(:,:,2,IBL), BUFFER_LOC(:,:,4:8,IBL)

    !  CALL CLOUDSC_SCC_CUF_K_CACHING<<<GRIDDIM,BLOCKDIM>>>(1, ICEND, NPROMA, NLEV, PTSPHY, PT_d, PQ_d, BUFFER_TMP_d, BUFFER_LOC_d, PVFA_d,  &
    !  CALL CLOUDSC_SCC_CUF_K_CACHING(1, NGPTOT, NPROMA, ICEND, NPROMA, NLEV, PTSPHY, PT, PQ, BUFFER_TMP, BUFFER_LOC, PVFA,  &
    ! CALL CLOUDSC_C_HIP(1, NGPTOT, NPROMA, ICEND, NPROMA, NLEV, PTSPHY, PT, PQ, &
    ! & TENDENCY_TMP_T, TENDENCY_TMP_Q, TENDENCY_TMP_A, TENDENCY_TMP_CLD, TENDENCY_LOC_T, TENDENCY_LOC_Q, &
    ! & TENDENCY_LOC_A, TENDENCY_LOC_CLD , PVFA,  &
    ! & PVFL, PVFI, PDYNA, PDYNL, PDYNI, PHRSW, PHRLW, PVERVEL, PAP, PAPH, PLSM, LDCUM, KTYPE, PLU,  &
    ! & PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER, PICRIT_AER, PRE_ICE, PCCN, PNICE,  &
    ! & PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG, PFCQSNG,  &
    ! & PFSQLTUR, PFSQITUR, PFPLSL, PFPLSN, PFHPSL, PFHPSN, NGPBLKS, YRECLDP=LOCAL_YRECLDP)
    CALL cloudsc_c_hip(ngptot, nproma, 1, nproma, nproma, nlev, ptsphy, pt, &
          & pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
          & tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
          & pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, &
          & pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
          & pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
          & pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn, & ! yrecldp, &
          & ngpblks, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, r3les, r3ies, r4les, &
          & r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
          & rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)
    !istat = cudaDeviceSynchronize()
   
    ! rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, r3les, r3ies, r4les, &
    !      & r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
    !      & rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2

    ! $omp end target
    
    CALL TIMER%THREAD_END(TID)
   
    !$omp end target data

    BUFFER_TMP(:,:,1,:) = TENDENCY_TMP_T
    BUFFER_TMP(:,:,3,:) = TENDENCY_TMP_Q
    BUFFER_TMP(:,:,2,:) = TENDENCY_TMP_A
    BUFFER_TMP(:,:,4:8,:) = TENDENCY_TMP_CLD
    BUFFER_LOC(:,:,1,:) = TENDENCY_LOC_T
    BUFFER_LOC(:,:,3,:) = TENDENCY_LOC_Q
    BUFFER_LOC(:,:,2,:) = TENDENCY_LOC_A
    BUFFER_LOC(:,:,4:8,:) = TENDENCY_LOC_CLD

    ! $omp end target data

    CALL TIMER%END()
    
    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)
    
    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)
    
  END SUBROUTINE CLOUDSC_DRIVER
  
END MODULE CLOUDSC_DRIVER_MOD
