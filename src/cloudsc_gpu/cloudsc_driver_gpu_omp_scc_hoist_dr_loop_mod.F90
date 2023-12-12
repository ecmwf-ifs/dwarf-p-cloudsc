! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
MODULE CLOUDSC_DRIVER_HOIST_DR_LOOP_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY: NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY: PERFORMANCE_TIMER, GET_THREAD_NUM
  
  ! USE YOMCST_CUF, ONLY: YOMCST_UPDATE_DEVICE
  ! USE YOETHF_CUF, ONLY: YOETHF_UPDATE_DEVICE
  
  USE CLOUDSC_HOIST_DR_LOOP_MOD, ONLY: CLOUDSC_HOIST_DR_LOOP
  
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE CLOUDSC_DRIVER_HOIST_DR_LOOP (NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, KFLDX, PTSPHY, PT, PQ, TENDENCY_CML,  &
  & TENDENCY_TMP, TENDENCY_LOC, BUFFER_CML, BUFFER_TMP, BUFFER_LOC, PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, PHRSW, PHRLW,  &
  & PVERVEL, PAP, PAPH, PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER, PICRIT_AER, PRE_ICE,  &
  & PCCN, PNICE, PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG, PFCQSNG, PFSQLTUR,  &
  & PFSQITUR, PFPLSL, PFPLSN, PFHPSL, PFHPSN)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel
    
    ! THIS IS THE CLOUDSC CUF VERSION
    
    ! USE cudafor
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
    
    INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND
    
    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP
    
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID    ! thread id from 0 .. NUMOMP - 1
    ! TYPE(TECLDP), DEVICE :: YRECLDP
    INTEGER :: istat
    
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

    !@cuf print *, 'executing SCC-CUF type: hoist - host side hoisted local arrays'
    
    IBL = 1      ! Useless statement to show the compiler that the sepcification part is over!
    
    !@cuf CALL YOMCST_UPDATE_DEVICE()
    !@cuf CALL YOETHF_UPDATE_DEVICE()
    
    IF (irank == 0) THEN
1003  FORMAT(5X, 'NUMPROC=', I0, ', NUMOMP=', I0, ', NGPTOTG=', I0, ', NPROMA=', I0, ', NGPBLKS=', I0)
      WRITE(0, 1003) NUMPROC, NUMOMP, NGPTOTG, NPROMA, NGPBLKS
    END IF
    
    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    
    ! $acc enter data create(ZFOEALFA, ZTP1, ZLI, ZA, ZAORIG, ZLIQFRAC, ZICEFRAC, ZQX, ZQX0,  &
    ! $acc &   ZPFPLSX, ZLNEG, ZQXN2D, ZQSMIX, ZQSLIQ, ZQSICE, ZFOEEWMT,  &
    ! $acc &   ZFOEEW, ZFOEELIQT)

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP
    
    ! $acc data &
    ! $acc copyin( &
    ! $acc   pt,pq,buffer_cml,buffer_tmp,pvfa, &
    ! $acc   pvfl,pvfi,pdyna,pdynl,pdyni,phrsw,phrlw,pvervel, &
    ! $acc   pap,paph,plsm,ldcum,ktype,plu,psnde, &
    ! $acc   pmfu,pmfd,pa,pclv,psupsat,plcrit_aer,picrit_aer, &
    ! $acc   pre_ice,pccn,pnice, yrecldp) &
    ! $acc copy( &
    ! $acc   buffer_loc,plude,pcovptot,prainfrac_toprfz) &
    ! $acc copyout( &
    ! $acc   pfsqlf,pfsqif,pfcqnng, &
    ! $acc   pfcqlng ,pfsqrf,pfsqsf,pfcqrng,pfcqsng,pfsqltur, &
    ! $acc   pfsqitur,pfplsl,pfplsn,pfhpsl,pfhpsn)  
    
    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)
    
    ! BLOCKDIM = DIM3(NPROMA, 1, 1)
    ! GRIDDIM = DIM3(1, 1, CEILING(REAL(NGPTOT) / REAL(NPROMA)))
    ! IBL = (JKGLO - 1) / NPROMA + 1
    ICEND = MIN(NPROMA, NGPTOT - JKGLO + 1)
    
    ! CALL CUF_HOIST<<<GRIDDIM,BLOCKDIM>>>(1, ICEND, NPROMA, NLEV, PTSPHY, PT, PQ, BUFFER_TMP, BUFFER_LOC,  &
    CALL CLOUDSC_HOIST_DR_LOOP(1, ICEND, NPROMA, NLEV, PTSPHY, PT, PQ, BUFFER_TMP, BUFFER_LOC,  &
    & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, PHRSW, PHRLW, PVERVEL, PAP, PAPH, PLSM, LDCUM, KTYPE,  &
    & PLU, PLUDE, PSNDE, PMFU, PMFD, PA, PCLV, PSUPSAT, PLCRIT_AER, PICRIT_AER, PRE_ICE, PCCN, PNICE,  &
    & PCOVPTOT, PRAINFRAC_TOPRFZ, PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, PFSQRF, PFSQSF, PFCQRNG, PFCQSNG,  &
    & PFSQLTUR, PFSQITUR, PFPLSL, PFPLSN, PFHPSL, PFHPSN, LOCAL_YRECLDP, NGPBLKS, ZFOEALFA,  &
    & ZTP1, ZLI, ZA, ZAORIG, ZLIQFRAC, ZICEFRAC,  &
    & ZQX, ZQX0, ZPFPLSX, ZLNEG, ZQXN2D, ZQSMIX,  &
    & ZQSLIQ, ZQSICE, ZFOEEWMT, ZFOEEW, ZFOEELIQT, NGPTOT, NPROMA)
    ! istat = cudaDeviceSynchronize()
    !---prognostic fields
    !-- arrays for aerosol-cloud interactions
    !---diagnostic output
    !---resulting fluxes
    
    
   ! $acc end data 
    
    
    
    !-- The "nowait" is here to get correct local timings (tloc) per thread
    !   i.e. we should not wait for slowest thread to finish before measuring tloc
    
    CALL TIMER%THREAD_END(TID)
    
    CALL TIMER%END()
    
    
    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)
    
    
    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)
    
    
  END SUBROUTINE CLOUDSC_DRIVER_HOIST_DR_LOOP
  
END MODULE CLOUDSC_DRIVER_HOIST_DR_LOOP_MOD
