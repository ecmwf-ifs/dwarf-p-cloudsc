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

  USE CLOUDSC_GPU_SCC_CUF_MOD, ONLY: CLOUDSC_SCC_CUF
  USE NLEV_MOD, ONLY : NLEV

  USE CUDAFOR

  IMPLICIT NONE


CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_CUF( &
     & NUMOMP, NPROMA, NGPTOT, NGPBLKS, NGPTOTG, KFLDX, PTSPHY, &
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

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NGPTOT, NGPBLKS, NGPTOTG
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

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICSTART, ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    INTEGER :: ISTAT

    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP
 
    TYPE(DIM3) :: GRIDDIM, BLOCKDIM

!!  device variables
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PT_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS) ! T at start of callpar
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PQ_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS) ! Q at start of callpar
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_CML_D(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_TMP_D(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: BUFFER_LOC_D(:,:,:,:) !!(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFA_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! CC from VDF scheme
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFL_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Liq from VDF scheme
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVFI_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Ice from VDF scheme
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNA_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! CC from Dynamics
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNL_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PDYNI_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PHRSW_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Short-wave heating rate
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PHRLW_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Long-wave heating rate
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PVERVEL_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)  !Vertical velocity
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PAP_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)      ! Pressure on full levels
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PAPH_D(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS) ! Pressure on half levels
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLSM_D(:,:) !!(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    LOGICAL,            DEVICE, ALLOCATABLE :: LDCUM_D(:,:) !!(NPROMA, NGPBLKS)    ! Convection active
    INTEGER(KIND=JPIM), DEVICE, ALLOCATABLE :: KTYPE_D(:,:) !!(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLU_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)      ! Conv. condensate
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLUDE_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained water
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PSNDE_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained snow
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PMFU_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux up
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PMFD_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux down
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PA_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PCLV_D(:,:,:,:) !!(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PSUPSAT_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PLCRIT_AER_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PICRIT_AER_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PRE_ICE_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PCCN_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB),    DEVICE, ALLOCATABLE :: PNICE_D(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PCOVPTOT_d(:,:,:) !!(NPROMA, NLEV, NGPBLKS)    ! Precip fraction
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PRAINFRAC_TOPRFZ_d(:,:) !!(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQLF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux of liquid
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQIF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux of ice
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQLNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for liq
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQNNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for ice
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQRF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Flux diagnostics
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQSF_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    !    for DDH, generic
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQRNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! rain
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFCQSNG_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)   ! snow
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQLTUR_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)  ! liquid flux due to VDF
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFSQITUR_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)  ! ice flux due to VDF

    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFPLSL_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! liq+rain sedim flux
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFPLSN_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! ice+snow sedim flux
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFHPSL_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: PFHPSN_d(:,:,:) !!(NPROMA, NLEV+1, NGPBLKS)    ! ice number concentration (cf. CCN)

    TYPE(TECLDP),DEVICE :: LOCAL_YRECLDP_d

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)

    ALLOCATE( PT_d(NPROMA, NLEV, NGPBLKS),PQ_d(NPROMA, NLEV, NGPBLKS),BUFFER_CML_d(NPROMA,NLEV,3+NCLV,NGPBLKS))
    ALLOCATE( BUFFER_TMP_d(NPROMA,NLEV,3+NCLV,NGPBLKS), BUFFER_LOC_d(NPROMA,NLEV,3+NCLV,NGPBLKS))
    ALLOCATE( PVFA_d(NPROMA, NLEV, NGPBLKS), PVFL_d(NPROMA, NLEV, NGPBLKS), PVFI_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PDYNA_d(NPROMA, NLEV, NGPBLKS), PDYNL_d(NPROMA, NLEV, NGPBLKS), PDYNI_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PHRSW_d(NPROMA, NLEV, NGPBLKS), PHRLW_d(NPROMA, NLEV, NGPBLKS), PVERVEL_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PAP_d(NPROMA, NLEV, NGPBLKS), PAPH_d(NPROMA, NLEV+1, NGPBLKS), PLSM_d(NPROMA, NGPBLKS) )
    ALLOCATE( LDCUM_d(NPROMA, NGPBLKS), KTYPE_d(NPROMA, NGPBLKS), PLU_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PLUDE_d(NPROMA, NLEV, NGPBLKS), PSNDE_d(NPROMA, NLEV, NGPBLKS), PMFU_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PMFD_d(NPROMA, NLEV, NGPBLKS), PA_d(NPROMA, NLEV, NGPBLKS), PCLV_d(NPROMA, NLEV, NCLV, NGPBLKS) )
    ALLOCATE( PSUPSAT_d(NPROMA, NLEV, NGPBLKS), PLCRIT_AER_d(NPROMA, NLEV, NGPBLKS), PICRIT_AER_d(NPROMA, NLEV, NGPBLKS) )
    ALLOCATE( PRE_ICE_d(NPROMA, NLEV, NGPBLKS), PCCN_d(NPROMA, NLEV, NGPBLKS),  PNICE_d(NPROMA, NLEV, NGPBLKS) )

    ALLOCATE(PCOVPTOT_d(NPROMA, NLEV, NGPBLKS), PRAINFRAC_TOPRFZ_d(NPROMA, NGPBLKS) )
    ALLOCATE(PFSQLF_d(NPROMA, NLEV+1, NGPBLKS), PFSQIF_d(NPROMA, NLEV+1, NGPBLKS), PFCQLNG_d(NPROMA, NLEV+1, NGPBLKS) )
    ALLOCATE(PFCQNNG_d(NPROMA, NLEV+1, NGPBLKS), PFSQRF_d(NPROMA, NLEV+1, NGPBLKS), PFSQSF_d(NPROMA, NLEV+1, NGPBLKS) )
    ALLOCATE(PFCQRNG_d(NPROMA, NLEV+1, NGPBLKS), PFCQSNG_d(NPROMA, NLEV+1, NGPBLKS) ,PFSQLTUR_d(NPROMA, NLEV+1, NGPBLKS) )
    ALLOCATE( PFSQITUR_d(NPROMA, NLEV+1, NGPBLKS))

    ALLOCATE(PFPLSL_d(NPROMA, NLEV+1, NGPBLKS), PFPLSN_d(NPROMA, NLEV+1, NGPBLKS))
    ALLOCATE(PFHPSL_d(NPROMA, NLEV+1, NGPBLKS), PFHPSN_d(NPROMA, NLEV+1, NGPBLKS) )


1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP

    pt_d=pt; pq_d=pq; buffer_tmp_d=buffer_tmp; buffer_loc_d = buffer_loc; pvfa_d=pvfa; pvfl_d=pvfl; pvfi_d=pvfi; pdyna_d=pdyna;
    pdynl_d=pdynl; pdyni_d=pdyni; phrsw_d=phrsw; phrlw_d=phrlw; pvervel_d=pvervel; pap_d=pap; paph_d=paph;
    plsm_d=plsm; ldcum_d=ldcum; ktype_d=ktype; plu_d=plu; plude_d=plude; psnde_d=psnde; pmfu_d=pmfu; pmfd_d=pmfd;
    pa_d=pa; pclv_d=pclv; psupsat_d=psupsat; plcrit_aer_d=plcrit_aer; picrit_aer_d=picrit_aer; pre_ice_d=pre_ice;
    pccn_d=pccn; pnice_d=pnice; pcovptot_d=pcovptot; prainfrac_toprfz_d=prainfrac_toprfz; pfsqlf_d=pfsqlf;
    pfsqif_d=pfsqif; pfcqnng_d=pfcqnng; pfcqlng_d=pfcqlng; pfsqrf_d=pfsqrf; pfsqsf_d=pfsqsf; pfcqrng_d=pfcqrng; pfcqsng_d=pfcqsng
    pfsqltur_d=pfsqltur; pfsqitur_d=pfsqitur; pfplsl_d=pfplsl; pfplsn_d=pfplsn; pfhpsl_d=pfhpsl; pfhpsn_d=pfhpsn
    local_yrecldp_d=local_yrecldp
    print *, 'In driver after copy to device !!!'
    
    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)


       ICSTART=1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

    GRIDDIM = DIM3(1,1,CEILING(REAL(NGPTOT)/REAL(NPROMA)))
    BLOCKDIM = DIM3(NPROMA,1,1)

    CALL CLOUDSC_SCC_CUF<<<GRIDDIM,BLOCKDIM , (nproma+1)*(25*8*3+1)>>> &
        & (ICSTART, ICEND, NPROMA, ngpblks, PTSPHY,&
        & pt_d,pq_d, &
        & buffer_tmp_d, &
        & buffer_loc_d, &
        & pvfa_d, pvfl_d, pvfi_d, pdyna_d, pdynl_d, pdyni_d, &
        & phrsw_d, phrlw_d, &
        & pvervel_d, pap_d, paph_d, &
        & plsm_d, ldcum_d, ktype_d, &
        & plu_d, plude_d, psnde_d, pmfu_d, pmfd_d, &
        !---prognostic fields
        & pa_d, pclv_d, psupsat_d, &
        !-- arrays for aerosol-cloud interactions
        & plcrit_aer_d, picrit_aer_d, &
        & pre_ice_d, &
        & pccn_d, pnice_d, &
        !---diagnostic output
        & pcovptot_d, prainfrac_toprfz_d, &
        !---resulting fluxes
        & pfsqlf_d, pfsqif_d, pfcqnng_d, pfcqlng_d, &
        & pfsqrf_d, pfsqsf_d, pfcqrng_d, pfcqsng_d, &
        & pfsqltur_d, pfsqitur_d, &
        & pfplsl_d, pfplsn_d, pfhpsl_d, pfhpsn_d, &
        &  YRECLDP=LOCAL_YRECLDP_d)

    istat = cudaDeviceSynchronize()
    
    CALL TIMER%THREAD_END(TID)

    print *, 'In driver after cloudsc !!!'
    buffer_tmp=buffer_tmp_d; buffer_loc = buffer_loc_d; 
    plude=plude_d; 
    pcovptot=pcovptot_d;
    pfsqlf=pfsqlf_d; pfsqif=pfsqif_d; pfcqnng=pfcqnng_d; pfcqlng=pfcqlng_d; pfsqrf=pfsqrf_d; pfsqsf=pfsqsf_d; pfcqrng=pfcqrng_d; 
    pfcqsng=pfcqsng_d; pfsqltur=pfsqltur_d; pfsqitur=pfsqitur_d; pfplsl=pfplsl_d; pfplsn=pfplsn_d; pfhpsl=pfhpsl_d; pfhpsn=pfhpsn_d;
    print *, 'In driver after copy back to host !!!'


    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_CUF

END MODULE CLOUDSC_DRIVER_GPU_SCC_CUF_MOD
