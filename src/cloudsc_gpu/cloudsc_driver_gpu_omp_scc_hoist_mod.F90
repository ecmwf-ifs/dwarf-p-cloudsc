! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_GPU_OMP_SCC_HOIST_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE CLOUDSC_GPU_OMP_SCC_HOIST_MOD, ONLY: CLOUDSC_SCC_HOIST

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_HOIST( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG, KFLDX, PTSPHY, &
     & PT, PQ, &
     & BUFFER_CML, BUFFER_TMP, BUFFER_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & PA, &
     & PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
     & YDOMCST, YDOETHF, YDECLDP )
    ! Driver routine that invokes the optimized CLAW-based CLOUDSC GPU kernel
    USE YOECLDP  , ONLY : TECLDP
    USE YOMCST   , ONLY : TOMCST
    USE YOETHF   , ONLY : TOETHF

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM)                                    :: KFLDX
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB), INTENT(IN)    :: PT(NPROMA, NLEV, NGPBLKS) ! T at start of callpar
    REAL(KIND=JPRB), INTENT(IN)    :: PQ(NPROMA, NLEV, NGPBLKS) ! Q at start of callpar
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_TMP(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB), INTENT(IN)    :: PVFA(NPROMA, NLEV, NGPBLKS)     ! CC from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFL(NPROMA, NLEV, NGPBLKS)     ! Liq from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFI(NPROMA, NLEV, NGPBLKS)     ! Ice from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNA(NPROMA, NLEV, NGPBLKS)    ! CC from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNL(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNI(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PHRSW(NPROMA, NLEV, NGPBLKS)    ! Short-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PHRLW(NPROMA, NLEV, NGPBLKS)    ! Long-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PVERVEL(NPROMA, NLEV, NGPBLKS)  !Vertical velocity
    REAL(KIND=JPRB), INTENT(IN)    :: PAP(NPROMA, NLEV, NGPBLKS)      ! Pressure on full levels
    REAL(KIND=JPRB), INTENT(IN)    :: PAPH(NPROMA, NLEV+1, NGPBLKS) ! Pressure on half levels
    REAL(KIND=JPRB), INTENT(IN)    :: PLSM(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    LOGICAL, INTENT(IN)            :: LDCUM(NPROMA, NGPBLKS)    ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    REAL(KIND=JPRB), INTENT(IN)    :: PLU(NPROMA, NLEV, NGPBLKS)      ! Conv. condensate
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained water
    REAL(KIND=JPRB), INTENT(IN)    :: PSNDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained snow
    REAL(KIND=JPRB), INTENT(IN)    :: PMFU(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux up
    REAL(KIND=JPRB), INTENT(IN)    :: PMFD(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux down
    REAL(KIND=JPRB), INTENT(IN)    :: PA(NPROMA, NLEV, NGPBLKS)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), INTENT(IN)    :: PCLV(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PSUPSAT(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PLCRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PICRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PRE_ICE(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PCCN(NPROMA, NLEV, NGPBLKS)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), INTENT(IN)    :: PNICE(NPROMA, NLEV, NGPBLKS)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB), INTENT(INOUT) :: PCOVPTOT(NPROMA, NLEV, NGPBLKS)    ! Precip fraction
    REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(NPROMA, NLEV+1, NGPBLKS)    ! Flux of liquid
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(NPROMA, NLEV+1, NGPBLKS)    ! Flux of ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(NPROMA, NLEV+1, NGPBLKS)    ! Flux diagnostics
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(NPROMA, NLEV+1, NGPBLKS)    !    for DDH, generic
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(NPROMA, NLEV+1, NGPBLKS)   ! rain
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(NPROMA, NLEV+1, NGPBLKS)   ! snow
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(NPROMA, NLEV+1, NGPBLKS)  ! liquid flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(NPROMA, NLEV+1, NGPBLKS)  ! ice flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(NPROMA, NLEV+1, NGPBLKS)    ! liq+rain sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(NPROMA, NLEV+1, NGPBLKS)    ! ice+snow sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(NPROMA, NLEV+1, NGPBLKS)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(NPROMA, NLEV+1, NGPBLKS)    ! ice number concentration (cf. CCN)

    TYPE(TOMCST), INTENT(IN)   :: YDOMCST
    TYPE(TOETHF), INTENT(IN)   :: YDOETHF
    TYPE(TECLDP), INTENT(IN)   :: YDECLDP

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
    INTEGER(KIND=JPIM) :: JL

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

#ifdef HAVE_OMP_TARGET
!$omp target enter data map(alloc: ZFOEALFA(1:NPROMA,1:NLEV+1,1:NGPBLKS), ZTP1(1:NPROMA,1:NLEV,1:NGPBLKS), ZLI(1:NPROMA,1:NLEV,1:NGPBLKS), ZA(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp &   ZAORIG(1:NPROMA,1:NLEV,1:NGPBLKS), ZLIQFRAC(1:NPROMA,1:NLEV,1:NGPBLKS), ZICEFRAC(1:NPROMA,1:NLEV,1:NGPBLKS), ZQX(1:NPROMA,1:NLEV,1:NCLV,1:NGPBLKS), ZQX0(1:NPROMA,1:NLEV,1:NCLV,1:NGPBLKS),  &
!$omp &   ZPFPLSX(1:NPROMA,1:NLEV+1,1:NCLV,1:NGPBLKS), ZLNEG(1:NPROMA,1:NLEV,1:NCLV,1:NGPBLKS), ZQXN2D(1:NPROMA,1:NLEV,1:NCLV,1:NGPBLKS), ZQSMIX(1:NPROMA,1:NLEV,1:NGPBLKS), ZQSLIQ(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp &   ZQSICE(1:NPROMA,1:NLEV,1:NGPBLKS), ZFOEEWMT(1:NPROMA,1:NLEV,1:NGPBLKS),  &
!$omp &   ZFOEEW(1:NPROMA,1:NLEV,1:NGPBLKS), ZFOEELIQT(1:NPROMA,1:NLEV,1:NGPBLKS))
#endif


#ifdef HAVE_OMP_TARGET
!$omp target data &
!$omp map(to: ptsphy,nproma,nlev,&
!$omp   pt(1:NPROMA,1:NLEV,1:NGPBLKS),pq(1:NPROMA,1:NLEV,1:NGPBLKS),&
!$omp   buffer_cml(1:NPROMA,1:NLEV,1:3+NCLV,1:NGPBLKS),buffer_tmp(1:NPROMA,1:NLEV,1:3+NCLV,1:NGPBLKS),pvfa(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   pvfl(1:NPROMA,1:NLEV,1:NGPBLKS),pvfi(1:NPROMA,1:NLEV,1:NGPBLKS),pdyna(1:NPROMA,1:NLEV,1:NGPBLKS),pdynl(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   pdyni(1:NPROMA,1:NLEV,1:NGPBLKS),phrsw(1:NPROMA,1:NLEV,1:NGPBLKS),phrlw(1:NPROMA,1:NLEV,1:NGPBLKS),pvervel(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   pap(1:NPROMA,1:NLEV,1:NGPBLKS),paph(1:NPROMA,1:NLEV+1,1:NGPBLKS),plsm(1:NPROMA,1:NGPBLKS),ldcum(1:NPROMA,1:NGPBLKS),ktype(1:NPROMA,1:NGPBLKS), &
!$omp   plu(1:NPROMA,1:NLEV,1:NGPBLKS),psnde(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   pmfu(1:NPROMA,1:NLEV,1:NGPBLKS),pmfd(1:NPROMA,1:NLEV,1:NGPBLKS),pa(1:NPROMA,1:NLEV,1:NGPBLKS),pclv(1:NPROMA,1:NLEV,1:NCLV,1:NGPBLKS),psupsat(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   plcrit_aer(1:NPROMA,1:NLEV,1:NGPBLKS),picrit_aer(1:NPROMA,1:NLEV,1:NGPBLKS), &
!$omp   pre_ice(1:NPROMA,1:NLEV,1:NGPBLKS),pccn(1:NPROMA,1:NLEV,1:NGPBLKS),pnice(1:NPROMA,1:NLEV,1:NGPBLKS), ydecldp, ydomcst, ydoethf) &
!$omp map(tofrom: &
!$omp   buffer_loc(1:NPROMA,1:NLEV,1:3+NCLV,1:NGPBLKS),plude(1:NPROMA,1:NLEV,1:NGPBLKS),pcovptot(1:NPROMA,1:NLEV,1:NGPBLKS)) &
!$omp map(from: &
!$omp   pfsqlf(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfsqif(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfcqnng(1:NPROMA,1:NLEV+1,1:NGPBLKS), &
!$omp   pfcqlng(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfsqrf(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfsqsf(1:NPROMA,1:NLEV+1,1:NGPBLKS), &
!$omp   pfcqrng(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfcqsng(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfsqltur(1:NPROMA,1:NLEV+1,1:NGPBLKS), &
!$omp   pfsqitur(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfplsl(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfplsn(1:NPROMA,1:NLEV+1,1:NGPBLKS), &
!$omp   pfhpsl(1:NPROMA,1:NLEV+1,1:NGPBLKS),pfhpsn(1:NPROMA,1:NLEV+1,1:NGPBLKS),prainfrac_toprfz(1:NPROMA,1:NGPBLKS))
#endif

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

!!$omp target teams loop bind(teams) thread_limit(64) private(ibl,icend,jkglo)
!$omp target teams distribute private(ibl,icend,jkglo) thread_limit(64)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

!!$omp loop bind(thread)
!!$omp parallel do simd private(JL) safelen(64)
!$omp parallel do private(JL)
      DO JL=1,ICEND
        CALL CLOUDSC_SCC_HOIST &
         & (1, ICEND, NPROMA, NLEV, PTSPHY,&
         & PT(:,:,IBL), PQ(:,:,IBL), &
         & BUFFER_TMP(:,:,1,IBL), BUFFER_TMP(:,:,3,IBL), BUFFER_TMP(:,:,2,IBL), BUFFER_TMP(:,:,4:8,IBL), &
         & BUFFER_LOC(:,:,1,IBL), BUFFER_LOC(:,:,3,IBL), BUFFER_LOC(:,:,2,IBL), BUFFER_LOC(:,:,4:8,IBL), &
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
         & YDOMCST, YDOETHF, YDECLDP, &
         & ZFOEALFA(:,:,IBL), ZTP1(:,:,IBL), ZLI(:,:,IBL), ZA(:,:,IBL), ZAORIG(:,:,IBL), &
         & ZLIQFRAC(:,:,IBL), ZICEFRAC(:,:,IBL), ZQX(:,:,:,IBL), ZQX0(:,:,:,IBL), ZPFPLSX(:,:,:,IBL), &
         & ZLNEG(:,:,:,IBL), ZQXN2D(:,:,:,IBL), ZQSMIX(:,:,IBL), ZQSLIQ(:,:,IBL), ZQSICE(:,:,IBL), &
         & ZFOEEWMT(:,:,IBL), ZFOEEW(:,:,IBL), ZFOEELIQT(:,:,IBL), JL=JL)
      ENDDO
!!$omp end loop
!$omp end parallel do
    ENDDO
!!$omp end target teams loop
!$omp end target teams distribute


    CALL TIMER%THREAD_END(TID)

#ifdef HAVE_OMP_TARGET
!$omp end target data
#endif

    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_HOIST

END MODULE CLOUDSC_DRIVER_GPU_OMP_SCC_HOIST_MOD
