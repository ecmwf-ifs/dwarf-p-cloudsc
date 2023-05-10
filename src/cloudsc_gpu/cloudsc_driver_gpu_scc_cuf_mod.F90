! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_GPU_SCC_CUF_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  !USE YOMCST_CUF,ONLY : YOMCST_UPDATE_DEVICE
  !USE YOETHF_CUF,ONLY : YOETHF_UPDATE_DEVICE

  USE CLOUDSC_GPU_SCC_CUF_MOD, ONLY: CLOUDSC_SCC_CUF
  USE NLEV_MOD, ONLY : NLEV

  ! USE CUDAFOR

  IMPLICIT NONE


CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_CUF( &
     & NUMOMP, NPROMA, NLEV_IN, NGPTOT, NGPBLKS, NGPTOTG, KFLDX, PTSPHY, &
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
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that invokes the optimized CLAW-based CLOUDSC GPU kernel

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NLEV_IN, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM)                                    :: KFLDX 
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB), INTENT(IN)    :: PT(NPROMA, NLEV_IN, NGPBLKS) ! T at start of callpar
    REAL(KIND=JPRB), INTENT(IN)    :: PQ(NPROMA, NLEV_IN, NGPBLKS) ! Q at start of callpar
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV_IN,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_TMP(NPROMA,NLEV_IN,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV_IN,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB), INTENT(IN)    :: PVFA(NPROMA, NLEV_IN, NGPBLKS)     ! CC from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFL(NPROMA, NLEV_IN, NGPBLKS)     ! Liq from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFI(NPROMA, NLEV_IN, NGPBLKS)     ! Ice from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNA(NPROMA, NLEV_IN, NGPBLKS)    ! CC from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNL(NPROMA, NLEV_IN, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNI(NPROMA, NLEV_IN, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PHRSW(NPROMA, NLEV_IN, NGPBLKS)    ! Short-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PHRLW(NPROMA, NLEV_IN, NGPBLKS)    ! Long-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PVERVEL(NPROMA, NLEV_IN, NGPBLKS)  !Vertical velocity
    REAL(KIND=JPRB), INTENT(IN)    :: PAP(NPROMA, NLEV_IN, NGPBLKS)      ! Pressure on full levels
    REAL(KIND=JPRB), INTENT(IN)    :: PAPH(NPROMA, NLEV_IN+1, NGPBLKS) ! Pressure on half levels
    REAL(KIND=JPRB), INTENT(IN)    :: PLSM(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    LOGICAL, INTENT(IN)            :: LDCUM(NPROMA, NGPBLKS)    ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    REAL(KIND=JPRB), INTENT(IN)    :: PLU(NPROMA, NLEV_IN, NGPBLKS)      ! Conv. condensate
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(NPROMA, NLEV_IN, NGPBLKS)    ! Conv. detrained water
    REAL(KIND=JPRB), INTENT(IN)    :: PSNDE(NPROMA, NLEV_IN, NGPBLKS)    ! Conv. detrained snow
    REAL(KIND=JPRB), INTENT(IN)    :: PMFU(NPROMA, NLEV_IN, NGPBLKS)     ! Conv. mass flux up
    REAL(KIND=JPRB), INTENT(IN)    :: PMFD(NPROMA, NLEV_IN, NGPBLKS)     ! Conv. mass flux down
    REAL(KIND=JPRB), INTENT(IN)    :: PA(NPROMA, NLEV_IN, NGPBLKS)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), INTENT(IN)    :: PCLV(NPROMA, NLEV_IN, NCLV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PSUPSAT(NPROMA, NLEV_IN, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PLCRIT_AER(NPROMA, NLEV_IN, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PICRIT_AER(NPROMA, NLEV_IN, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PRE_ICE(NPROMA, NLEV_IN, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PCCN(NPROMA, NLEV_IN, NGPBLKS)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), INTENT(IN)    :: PNICE(NPROMA, NLEV_IN, NGPBLKS)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB), INTENT(INOUT) :: PCOVPTOT(NPROMA, NLEV_IN, NGPBLKS)    ! Precip fraction
    REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(NPROMA, NLEV_IN+1, NGPBLKS)    ! Flux of liquid
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(NPROMA, NLEV_IN+1, NGPBLKS)    ! Flux of ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(NPROMA, NLEV_IN+1, NGPBLKS)   ! -ve corr for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(NPROMA, NLEV_IN+1, NGPBLKS)   ! -ve corr for ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(NPROMA, NLEV_IN+1, NGPBLKS)    ! Flux diagnostics
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(NPROMA, NLEV_IN+1, NGPBLKS)    !    for DDH, generic
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(NPROMA, NLEV_IN+1, NGPBLKS)   ! rain
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(NPROMA, NLEV_IN+1, NGPBLKS)   ! snow
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(NPROMA, NLEV_IN+1, NGPBLKS)  ! liquid flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(NPROMA, NLEV_IN+1, NGPBLKS)  ! ice flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(NPROMA, NLEV_IN+1, NGPBLKS)    ! liq+rain sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(NPROMA, NLEV_IN+1, NGPBLKS)    ! ice+snow sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(NPROMA, NLEV_IN+1, NGPBLKS)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(NPROMA, NLEV_IN+1, NGPBLKS)    ! ice number concentration (cf. CCN)

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICSTART, ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    INTEGER :: ISTAT

    ! Local copy of cloud parameters for offload
    !TYPE(TECLDP), DEVICE :: LOCAL_YRECLDP
 
    !TYPE(DIM3) :: GRIDDIM, BLOCKDIM

!!  device variables
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PT_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS) ! T at start of callpar
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PQ_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS) ! Q at start of callpar
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_CML_d(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_TMP_d(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_LOC_d(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFA_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! CC from VDF scheme
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFL_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Liq from VDF scheme
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFI_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Ice from VDF scheme
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNA_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! CC from Dynamics
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNL_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNI_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PHRSW_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Short-wave heating rate
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PHRLW_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Long-wave heating rate
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVERVEL_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)  !Vertical velocity
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PAP_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)      ! Pressure on full levels
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PAPH_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS) ! Pressure on half levels
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLSM_d(:,:) !!(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    !LOGICAL,            DEVICE, ALLOCATABLE :: LDCUM_d(:,:) !!(NPROMA, NGPBLKS)    ! Convection active
    !INTEGER(KIND=JPIM), DEVICE, ALLOCATABLE :: KTYPE_d(:,:) !!(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLU_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)      ! Conv. condensate
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLUDE_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained water
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PSNDE_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained snow
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PMFU_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux up
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PMFD_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux down
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PA_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)       ! Original Cloud fraction (t)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PCLV_d(:,:,:,:) !!(NPROMA, NLEV, NCLV, NGPBLKS)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PSUPSAT_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLCRIT_AER_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PICRIT_AER_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PRE_ICE_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PCCN_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! liquid cloud condensation nuclei
    !REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PNICE_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! ice number concentration (cf. CCN)

    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PCOVPTOT_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Precip fraction
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PRAINFRAC_TOPRFZ_d(:,:) !!(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQLF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux of liquid
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQIF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux of ice
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQLNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for liq
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQNNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for ice
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQRF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux diagnostics
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQSF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    !    for DDH, generic
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQRNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! rain
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQSNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! snow
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQLTUR_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)  ! liquid flux due to VDF
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQITUR_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)  ! ice flux due to VDF

    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFPLSL_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! liq+rain sedim flux
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFPLSN_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! ice+snow sedim flux
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFHPSL_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Enthalpy flux for liq
    !REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFHPSN_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! ice number concentration (cf. CCN)

#include "abor1.intfb.h"

    ! Transfer global module-scope parameters to constant device memory
    !CALL YOMCST_UPDATE_DEVICE()
    !CALL YOETHF_UPDATE_DEVICE()

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)

    IF (NLEV_IN /= NLEV) THEN
      CALL ABOR1('ERROR: Vertical dimension NLEV does not equal parametrised constant NLEV=137')
    END IF

    !ALLOCATE( PT_d(NPROMA, NLEV, NGPBLKS),PQ_d(NPROMA, NLEV, NGPBLKS),BUFFER_CML_d(NPROMA,NLEV,3+NCLV,NGPBLKS))
    !ALLOCATE( BUFFER_TMP_d(NPROMA,NLEV,3+NCLV,NGPBLKS), BUFFER_LOC_d(NPROMA,NLEV,3+NCLV,NGPBLKS))
    !ALLOCATE( PVFA_d(NPROMA, NLEV, NGPBLKS), PVFL_d(NPROMA, NLEV, NGPBLKS), PVFI_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PDYNA_d(NPROMA, NLEV, NGPBLKS), PDYNL_d(NPROMA, NLEV, NGPBLKS), PDYNI_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PHRSW_d(NPROMA, NLEV, NGPBLKS), PHRLW_d(NPROMA, NLEV, NGPBLKS), PVERVEL_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PAP_d(NPROMA, NLEV, NGPBLKS), PAPH_d(NPROMA, NLEV+1, NGPBLKS), PLSM_d(NPROMA, NGPBLKS) )
    !ALLOCATE( LDCUM_d(NPROMA, NGPBLKS), KTYPE_d(NPROMA, NGPBLKS), PLU_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PLUDE_d(NPROMA, NLEV, NGPBLKS), PSNDE_d(NPROMA, NLEV, NGPBLKS), PMFU_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PMFD_d(NPROMA, NLEV, NGPBLKS), PA_d(NPROMA, NLEV, NGPBLKS), PCLV_d(NPROMA, NLEV, NCLV, NGPBLKS) )
    !ALLOCATE( PSUPSAT_d(NPROMA, NLEV, NGPBLKS), PLCRIT_AER_d(NPROMA, NLEV, NGPBLKS), PICRIT_AER_d(NPROMA, NLEV, NGPBLKS) )
    !ALLOCATE( PRE_ICE_d(NPROMA, NLEV, NGPBLKS), PCCN_d(NPROMA, NLEV, NGPBLKS),  PNICE_d(NPROMA, NLEV, NGPBLKS) )

    !ALLOCATE(PCOVPTOT_d(NPROMA, NLEV, NGPBLKS), PRAINFRAC_TOPRFZ_d(NPROMA, NGPBLKS) )
    !ALLOCATE(PFSQLF_d(NPROMA, NLEV+1, NGPBLKS), PFSQIF_d(NPROMA, NLEV+1, NGPBLKS), PFCQLNG_d(NPROMA, NLEV+1, NGPBLKS) )
    !ALLOCATE(PFCQNNG_d(NPROMA, NLEV+1, NGPBLKS), PFSQRF_d(NPROMA, NLEV+1, NGPBLKS), PFSQSF_d(NPROMA, NLEV+1, NGPBLKS) )
    !ALLOCATE(PFCQRNG_d(NPROMA, NLEV+1, NGPBLKS), PFCQSNG_d(NPROMA, NLEV+1, NGPBLKS) ,PFSQLTUR_d(NPROMA, NLEV+1, NGPBLKS) )
    !ALLOCATE( PFSQITUR_d(NPROMA, NLEV+1, NGPBLKS))

    !ALLOCATE(PFPLSL_d(NPROMA, NLEV+1, NGPBLKS), PFPLSN_d(NPROMA, NLEV+1, NGPBLKS))
    !ALLOCATE(PFHPSL_d(NPROMA, NLEV+1, NGPBLKS), PFHPSN_d(NPROMA, NLEV+1, NGPBLKS) )


1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    !LOCAL_YRECLDP = YRECLDP

    !pt_d=pt; pq_d=pq; buffer_tmp_d=buffer_tmp; buffer_loc_d = buffer_loc; pvfa_d=pvfa; pvfl_d=pvfl; pvfi_d=pvfi; pdyna_d=pdyna;
    !pdynl_d=pdynl; pdyni_d=pdyni; phrsw_d=phrsw; phrlw_d=phrlw; pvervel_d=pvervel; pap_d=pap; paph_d=paph;
    !plsm_d=plsm; ldcum_d=ldcum; ktype_d=ktype; plu_d=plu; plude_d=plude; psnde_d=psnde; pmfu_d=pmfu; pmfd_d=pmfd;
    !pa_d=pa; pclv_d=pclv; psupsat_d=psupsat; plcrit_aer_d=plcrit_aer; picrit_aer_d=picrit_aer; pre_ice_d=pre_ice;
    !pccn_d=pccn; pnice_d=pnice; pcovptot_d=pcovptot; prainfrac_toprfz_d=prainfrac_toprfz; pfsqlf_d=pfsqlf;
    !pfsqif_d=pfsqif; pfcqnng_d=pfcqnng; pfcqlng_d=pfcqlng; pfsqrf_d=pfsqrf; pfsqsf_d=pfsqsf; pfcqrng_d=pfcqrng; pfcqsng_d=pfcqsng
    !pfsqltur_d=pfsqltur; pfsqitur_d=pfsqitur; pfplsl_d=pfplsl; pfplsn_d=pfplsn; pfhpsl_d=pfhpsl; pfhpsn_d=pfhpsn
    
    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)


    ICSTART=1
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

    !GRIDDIM = DIM3(1,1,CEILING(REAL(NGPTOT)/REAL(NPROMA)))
    !BLOCKDIM = DIM3(NPROMA,1,1)

    ! CALL CLOUDSC_SCC_CUF<<<GRIDDIM,BLOCKDIM >>> &
    CALL CLOUDSC_SCC_CUF &    
        & (NGPTOT, NPROMA, ICSTART, ICEND, NPROMA, NGPBLKS, PTSPHY,&
        & PT,PQ, &
        & BUFFER_TMP, &
        & BUFFER_LOC, &
        & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
        & PHRSW, PHRLW, &
        & PVERVEL, PAP, PAPH, &
        & PLSM, LDCUM, KTYPE, &
        & PLU, PLUDE, PSNDE, PMFU, PMFD, &
        !---prognostic fields
        & PA, PCLV, PSUPSAT, &
        !-- arrays for aerosol-cloud interactions
        & PLCRIT_AER, PICRIT_AER, &
        & PRE_ICE, &
        & PCCN, PNICE, &
        !---diagnostic output
        & PCOVPTOT, PRAINFRAC_TOPRFZ, &
        !---resulting fluxes
        & PFSQLF, PFSQIF, PFCQNNG, PFCQLNG, &
        & PFSQRF, PFSQSF, PFCQRNG, PFCQSNG, &
        & PFSQLTUR, PFSQITUR, &
        & PFPLSL, PFPLSN, PFHPSL, PFHPSN, &
        & YRECLDP )

    !ISTAT = cudaDeviceSynchronize()
    
    CALL TIMER%THREAD_END(TID)

    !buffer_tmp=buffer_tmp_d; buffer_loc = buffer_loc_d; 
    !plude=plude_d;
    !pcovptot=pcovptot_d;
    !pfsqlf=pfsqlf_d; pfsqif=pfsqif_d; pfcqnng=pfcqnng_d; pfcqlng=pfcqlng_d; pfsqrf=pfsqrf_d; pfsqsf=pfsqsf_d; pfcqrng=pfcqrng_d; 
    !pfcqsng=pfcqsng_d; pfsqltur=pfsqltur_d; pfsqitur=pfsqitur_d; pfplsl=pfplsl_d; pfplsn=pfplsn_d; pfhpsl=pfhpsl_d; pfhpsn=pfhpsn_d;

    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_CUF

END MODULE CLOUDSC_DRIVER_GPU_SCC_CUF_MOD
