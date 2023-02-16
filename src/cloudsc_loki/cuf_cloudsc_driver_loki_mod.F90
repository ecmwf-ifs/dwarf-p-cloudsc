! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CUF_CLOUDSC_DRIVER_LOKI_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE YOMCST_CUF,ONLY : YOMCST_UPDATE_DEVICE
  USE YOETHF_CUF,ONLY : YOETHF_UPDATE_DEVICE

  USE CLOUDSC_CUF_MOD, ONLY : CLOUDSC_CUF

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CUF_CLOUDSC_DRIVER( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, KFLDX, PTSPHY, &
     & PT, PQ, &
     & TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
     & BUFFER_CML, BUFFER_TMP, BUFFER_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    ! THIS IS THE CLOUDSC CUF VERSION

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM), INTENT(IN)    :: KFLDX
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(NPROMA,NLEV,NGPBLKS)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(NPROMA,NLEV,NGPBLKS)    ! Q at start of callpar
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(NGPBLKS) ! cumulative tendency used for final output
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_TMP(NGPBLKS) ! cumulative tendency used as input
    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(NGPBLKS) ! local tendency from cloud scheme
    REAL(KIND=JPRB),    INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(INOUT) :: BUFFER_TMP(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB),    INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFA(NPROMA,NLEV,NGPBLKS)  ! CC from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFL(NPROMA,NLEV,NGPBLKS)  ! Liq from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFI(NPROMA,NLEV,NGPBLKS)  ! Ice from VDF scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNA(NPROMA,NLEV,NGPBLKS) ! CC from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNL(NPROMA,NLEV,NGPBLKS) ! Liq from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNI(NPROMA,NLEV,NGPBLKS) ! Liq from Dynamics
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRSW(NPROMA,NLEV,NGPBLKS) ! Short-wave heating rate
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRLW(NPROMA,NLEV,NGPBLKS) ! Long-wave heating rate
    REAL(KIND=JPRB),    INTENT(IN)    :: PVERVEL(NPROMA,NLEV,NGPBLKS) !Vertical velocity
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(NPROMA,NLEV,NGPBLKS)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(NPROMA,NLEV+1,NGPBLKS)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLSM(NPROMA,NGPBLKS)    ! Land fraction (0-1)
    LOGICAL        ,    INTENT(IN)    :: LDCUM(NPROMA,NGPBLKS)   ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN)    :: KTYPE(NPROMA,NGPBLKS)   ! Convection type 0,1,2
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(NPROMA,NLEV,NGPBLKS)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(NPROMA,NLEV,NGPBLKS) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PSNDE(NPROMA,NLEV,NGPBLKS) ! Conv. detrained snow
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(NPROMA,NLEV,NGPBLKS)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(NPROMA,NLEV,NGPBLKS)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(NPROMA,NLEV,NGPBLKS)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(NPROMA,NLEV,NCLV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PLCRIT_AER(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PICRIT_AER(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PRE_ICE(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCCN(NPROMA,NLEV,NGPBLKS)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB),    INTENT(IN)    :: PNICE(NPROMA,NLEV,NGPBLKS)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(NPROMA,NLEV,NGPBLKS) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PRAINFRAC_TOPRFZ(NPROMA,NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLF(NPROMA,NLEV+1,NGPBLKS)  ! Flux of liquid
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQIF(NPROMA,NLEV+1,NGPBLKS)  ! Flux of ice
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQLNG(NPROMA,NLEV+1,NGPBLKS) ! -ve corr for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQNNG(NPROMA,NLEV+1,NGPBLKS) ! -ve corr for ice
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQRF(NPROMA,NLEV+1,NGPBLKS)  ! Flux diagnostics
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQSF(NPROMA,NLEV+1,NGPBLKS)  !    for DDH, generic
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQRNG(NPROMA,NLEV+1,NGPBLKS) ! rain
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQSNG(NPROMA,NLEV+1,NGPBLKS) ! snow
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLTUR(NPROMA,NLEV+1,NGPBLKS) ! liquid flux due to VDF
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQITUR(NPROMA,NLEV+1,NGPBLKS) ! ice flux due to VDF
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(NPROMA,NLEV+1,NGPBLKS) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(NPROMA,NLEV+1,NGPBLKS) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(NPROMA,NLEV+1,NGPBLKS) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(NPROMA,NLEV+1,NGPBLKS) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND

    ! Local copy of cloud parameters for offload
    !TYPE(TECLDP) :: LOCAL_YRECLDP

    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    IBL = 1  ! Useless statement to show the compiler that the sepcification part is over!

    !@cuf CALL YOMCST_UPDATE_DEVICE()
    !@cuf CALL YOETHF_UPDATE_DEVICE()

    if (irank == 0) then
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    !LOCAL_YRECLDP = YRECLDP

    !$loki data

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

      CALL CLOUDSC_CUF &
       & (    1,    ICEND,    NPROMA,  NLEV,&
       & PTSPHY,&
       & PT(:,:,IBL), PQ(:,:,IBL), &
       & BUFFER_TMP(:,:,:,IBL), BUFFER_LOC(:,:,:,IBL), &
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
       & YRECLDP)

#ifndef CLOUDSC_GPU_TIMING
      ! Log number of columns processed by this thread (OpenMP mode)
      CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
#endif
    ENDDO

    !-- The "nowait" is here to get correct local timings (tloc) per thread
    !   i.e. we should not wait for slowest thread to finish before measuring tloc
    !$omp end do nowait

    CALL TIMER%THREAD_END(TID)

    !$omp end parallel

    !$loki end data

    CALL TIMER%END()

#ifdef CLOUDSC_GPU_TIMING
    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER % THREAD_LOG(TID=TID, IGPC=NGPTOT)
#endif

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CUF_CLOUDSC_DRIVER

END MODULE CUF_CLOUDSC_DRIVER_LOKI_MOD
