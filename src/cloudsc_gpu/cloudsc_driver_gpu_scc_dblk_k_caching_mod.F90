! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_GPU_SCC_DBLK_K_CACHING_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE CLOUDSC_GPU_SCC_K_CACHING_MOD, ONLY: CLOUDSC_SCC_K_CACHING
  USE OPENACC
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_DBLK_K_CACHING( &
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
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that invokes the optimized CLAW-based CLOUDSC GPU kernel

    INTEGER(KIND=JPIM)             :: JL
    INTEGER(KIND=JPIM)             :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM)             :: KFLDX
    REAL(KIND=JPRB)                :: PTSPHY       ! Physics timestep
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

    ! Temporary buffers used for double blocked loop
    ! copyin
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pt_block
	!$acc declare device_resident(pt_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pq_block
	!$acc declare device_resident(pq_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: buffer_tmp_block
	!$acc declare device_resident(buffer_tmp_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfa_block
	!$acc declare device_resident(pvfa_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfl_block
	!$acc declare device_resident(pvfl_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfi_block
	!$acc declare device_resident(pvfi_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdyna_block
	!$acc declare device_resident(pdyna_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdynl_block
	!$acc declare device_resident(pdynl_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdyni_block
	!$acc declare device_resident(pdyni_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: phrsw_block
	!$acc declare device_resident(phrsw_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: phrlw_block
	!$acc declare device_resident(phrlw_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvervel_block
	!$acc declare device_resident(pvervel_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pap_block
	!$acc declare device_resident(pap_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: paph_block
	!$acc declare device_resident(paph_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:) :: plsm_block
	!$acc declare device_resident(plsm_block)
    LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ldcum_block
	!$acc declare device_resident(ldcum_block)
    INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:) :: ktype_block
	!$acc declare device_resident(ktype_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plu_block
	!$acc declare device_resident(plu_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: psnde_block
	!$acc declare device_resident(psnde_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pmfu_block
	!$acc declare device_resident(pmfu_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pmfd_block
	!$acc declare device_resident(pmfd_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pa_block
	!$acc declare device_resident(pa_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: pclv_block
	!$acc declare device_resident(pclv_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: psupsat_block
	!$acc declare device_resident(psupsat_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plcrit_aer_block
	!$acc declare device_resident(plcrit_aer_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: picrit_aer_block
	!$acc declare device_resident(picrit_aer_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pre_ice_block
	!$acc declare device_resident(pre_ice_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pccn_block
	!$acc declare device_resident(pccn_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pnice_block
	!$acc declare device_resident(pnice_block)
    ! copy
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: buffer_loc_block
	!$acc declare device_resident(buffer_loc_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plude_block
	!$acc declare device_resident(plude_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pcovptot_block
	!$acc declare device_resident(pcovptot_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:) :: prainfrac_toprfz_block
	!$acc declare device_resident(prainfrac_toprfz_block)
    ! copyout
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqlf_block
	!$acc declare device_resident(pfsqlf_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqif_block
	!$acc declare device_resident(pfsqif_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqnng_block
	!$acc declare device_resident(pfcqnng_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqlng_block
	!$acc declare device_resident(pfcqlng_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqrf_block
	!$acc declare device_resident(pfsqrf_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqsf_block
	!$acc declare device_resident(pfsqsf_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqrng_block
	!$acc declare device_resident(pfcqrng_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqsng_block
	!$acc declare device_resident(pfcqsng_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqltur_block
	!$acc declare device_resident(pfsqltur_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqitur_block
	!$acc declare device_resident(pfsqitur_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfplsl_block
	!$acc declare device_resident(pfplsl_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfplsn_block
	!$acc declare device_resident(pfplsn_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfhpsl_block
	!$acc declare device_resident(pfhpsl_block)
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfhpsn_block
	!$acc declare device_resident(pfhpsn_block)
    
    
    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP

    ! double blocking variables
    INTEGER(KIND=JPIM) :: BLOCK_BUFFER_SIZE     ! block size for blocks in outer loop
    INTEGER(KIND=JPIM) :: BLOCK_COUNT           ! number of blocks
    INTEGER(KIND=JPIM) :: BLOCK_IDX             ! idx of current block in [1,BLOCK_COUNT]
    INTEGER(KIND=JPIM) :: BLOCK_START           ! start of current block in [1,NGPBLKS]
    INTEGER(KIND=JPIM) :: BLOCK_END             ! end of current block in [1,NGPBLKS]
    INTEGER(KIND=JPIM) :: IBLLOC                ! local loop idx inside inner block loop
 

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
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

    ! BLOCK SIZES
    BLOCK_BUFFER_SIZE=MIN(900, NGPBLKS)
    BLOCK_COUNT=(NGPBLKS+BLOCK_BUFFER_SIZE-1)/BLOCK_BUFFER_SIZE
    
    print *, 'BLOCK_BUFFER_SIZE=', BLOCK_BUFFER_SIZE
    print *, 'BLOCK_COUNT=', BLOCK_COUNT
   
    ! buffer allocations
   !copyin
    ALLOCATE(pt_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pq_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(buffer_tmp_block(NPROMA,NLEV,3+NCLV,BLOCK_BUFFER_SIZE))
    ALLOCATE(pvfa_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pvfl_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pvfi_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pdyna_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pdynl_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pdyni_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(phrsw_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(phrlw_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pvervel_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pap_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(paph_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(plsm_block(NPROMA, BLOCK_BUFFER_SIZE))
    ALLOCATE(ldcum_block(NPROMA, BLOCK_BUFFER_SIZE))
    ALLOCATE(ktype_block(NPROMA, BLOCK_BUFFER_SIZE))
    ALLOCATE(plu_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(psnde_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pmfu_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pmfd_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pa_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pclv_block(NPROMA, NLEV, NCLV, BLOCK_BUFFER_SIZE))
    ALLOCATE(psupsat_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(plcrit_aer_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(picrit_aer_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pre_ice_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pccn_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pnice_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ! copy
    ALLOCATE(buffer_loc_block(NPROMA,NLEV,3+NCLV,BLOCK_BUFFER_SIZE))
    ALLOCATE(plude_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(pcovptot_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    ALLOCATE(prainfrac_toprfz_block(NPROMA, BLOCK_BUFFER_SIZE))
    !copyout
    ALLOCATE(pfsqlf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfsqif_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfcqnng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfcqlng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfsqrf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfsqsf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfcqrng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfcqsng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfsqltur_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfsqitur_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfplsl_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfplsn_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfhpsl_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    ALLOCATE(pfhpsn_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    
    
    
    DO BLOCK_IDX=0, BLOCK_COUNT-1
      BLOCK_START=BLOCK_IDX*BLOCK_BUFFER_SIZE+1
      BLOCK_END=MIN((BLOCK_IDX+1)*BLOCK_BUFFER_SIZE, NGPBLKS)
    
      ! data to device
      !$acc host_data &
        !$acc use_device(pt_block, &
     	!$acc pq_block, &
     	!$acc buffer_tmp_block, &
    !  	!$acc BUFFER_CML_block, &
     	!$acc pvfa_block, &
     	!$acc pvfl_block, &
     	!$acc pvfi_block, &
     	!$acc pdyna_block, &
     	!$acc pdynl_block, &
     	!$acc pdyni_block, &
     	!$acc phrsw_block, &
     	!$acc phrlw_block, &
     	!$acc pvervel_block, &
     	!$acc pap_block, &
     	!$acc paph_block, &
     	!$acc plsm_block, &
     	!$acc ldcum_block, &
     	!$acc ktype_block, &
     	!$acc plu_block, &
     	!$acc psnde_block, &
     	!$acc pmfu_block, &
     	!$acc pmfd_block, &
     	!$acc pa_block, &
     	!$acc pclv_block, &
     	!$acc psupsat_block, &
     	!$acc plcrit_aer_block, &
     	!$acc picrit_aer_block, &
     	!$acc pre_ice_block, &
     	!$acc pccn_block, &
     	!$acc pnice_block, &
      ! copy
     	!$acc buffer_loc_block, &
     	!$acc plude_block, &
     	!$acc pcovptot_block, &
     	!$acc prainfrac_toprfz_block)

      !copyin
        call acc_memcpy_to_device(pt_block, pt(:,:, BLOCK_START:BLOCK_END), SIZEOF(pt(:,:, BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pq_block, pq(:,:,BLOCK_START:BLOCK_END), SIZEOF(pq(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(buffer_tmp_block, buffer_tmp(:,:,:,BLOCK_START:BLOCK_END), SIZEOF(buffer_tmp(:,:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pvfa_block, pvfa(:,:,BLOCK_START:BLOCK_END), SIZEOF(pvfa(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pvfl_block, pvfl(:,:,BLOCK_START:BLOCK_END), SIZEOF(pvfl(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pvfi_block, pvfi(:,:,BLOCK_START:BLOCK_END), SIZEOF(pvfi(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pdyna_block, pdyna(:,:,BLOCK_START:BLOCK_END), SIZEOF(pdyna(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pdynl_block, pdynl(:,:,BLOCK_START:BLOCK_END), SIZEOF(pdynl(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pdyni_block, pdyni(:,:,BLOCK_START:BLOCK_END), SIZEOF(pdyni(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(phrsw_block, phrsw(:,:,BLOCK_START:BLOCK_END), SIZEOF(phrsw(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(phrlw_block, phrlw(:,:,BLOCK_START:BLOCK_END), SIZEOF(phrlw(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pvervel_block, pvervel(:,:,BLOCK_START:BLOCK_END), SIZEOF(pvervel(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pap_block, pap(:,:,BLOCK_START:BLOCK_END), SIZEOF(pap(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(paph_block, paph(:,:,BLOCK_START:BLOCK_END), SIZEOF(paph(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(plsm_block, plsm(:,BLOCK_START:BLOCK_END), SIZEOF(plsm(:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(ldcum_block, ldcum(:,BLOCK_START:BLOCK_END), SIZEOF(ldcum(:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(ktype_block, ktype(:,BLOCK_START:BLOCK_END), SIZEOF(ktype(:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(plu_block, plu(:,:,BLOCK_START:BLOCK_END), SIZEOF(plu(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(psnde_block, psnde(:,:,BLOCK_START:BLOCK_END), SIZEOF(psnde(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pmfu_block, pmfu(:,:,BLOCK_START:BLOCK_END), SIZEOF(pmfu(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pmfd_block, pmfd(:,:,BLOCK_START:BLOCK_END), SIZEOF(pmfd(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pa_block, pa(:,:,BLOCK_START:BLOCK_END), SIZEOF(pa(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pclv_block, pclv(:,:,:,BLOCK_START:BLOCK_END), SIZEOF(pclv(:,:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(psupsat_block, psupsat(:,:,BLOCK_START:BLOCK_END), SIZEOF(psupsat(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(plcrit_aer_block, plcrit_aer(:,:,BLOCK_START:BLOCK_END), SIZEOF(plcrit_aer(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(picrit_aer_block, picrit_aer(:,:,BLOCK_START:BLOCK_END), SIZEOF(picrit_aer(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pre_ice_block, pre_ice(:,:,BLOCK_START:BLOCK_END), SIZEOF(pre_ice(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pccn_block, pccn(:,:,BLOCK_START:BLOCK_END), SIZEOF(pccn(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pnice_block, pnice(:,:,BLOCK_START:BLOCK_END), SIZEOF(pnice(:,:,BLOCK_START:BLOCK_END)))
    !   copy)
        call acc_memcpy_to_device(buffer_loc_block, buffer_loc(:,:,:,BLOCK_START:BLOCK_END), SIZEOF(buffer_loc(:,:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(plude_block, plude(:,:,BLOCK_START:BLOCK_END), SIZEOF(plude(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(pcovptot_block, pcovptot(:,:,BLOCK_START:BLOCK_END), SIZEOF(pcovptot(:,:,BLOCK_START:BLOCK_END)))
        call acc_memcpy_to_device(prainfrac_toprfz_block, prainfrac_toprfz(:,BLOCK_START:BLOCK_END), SIZEOF(prainfrac_toprfz(:,BLOCK_START:BLOCK_END)))
      !$acc end host_data
    
    
    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_ACC_START(TID)


    !$acc parallel loop gang vector_length(NPROMA)
    DO IBLLOC=1, BLOCK_BUFFER_SIZE ! just a way to loop over NGPBLKS
        IBL= BLOCK_BUFFER_SIZE*BLOCK_IDX +IBLLOC
        JKGLO=(IBL-1)*NPROMA+1
        ICEND=MIN(NPROMA, NGPTOT-JKGLO+1)

      !$acc loop vector
      DO JL=1,ICEND
       CALL CLOUDSC_SCC_K_CACHING &
        & (1, ICEND, NPROMA, NLEV, PTSPHY,&
        & PT_BLOCK(:,:,IBLLOC), PQ_BLOCK(:,:,IBLLOC), &
        & BUFFER_TMP_BLOCK(:,:,1,IBLLOC), BUFFER_TMP_BLOCK(:,:,3,IBLLOC), BUFFER_TMP_BLOCK(:,:,2,IBLLOC), BUFFER_TMP_BLOCK(:,:,4:8,IBLLOC), &
        & BUFFER_LOC_BLOCK(:,:,1,IBLLOC), BUFFER_LOC_BLOCK(:,:,3,IBLLOC), BUFFER_LOC_BLOCK(:,:,2,IBLLOC), BUFFER_LOC_BLOCK(:,:,4:8,IBLLOC), &
        & PVFA_BLOCK(:,:,IBLLOC), PVFL_BLOCK(:,:,IBLLOC), PVFI_BLOCK(:,:,IBLLOC), PDYNA_BLOCK(:,:,IBLLOC), PDYNL_BLOCK(:,:,IBLLOC), PDYNI_BLOCK(:,:,IBLLOC), &
        & PHRSW_BLOCK(:,:,IBLLOC),    PHRLW_BLOCK(:,:,IBLLOC),&
        & PVERVEL_BLOCK(:,:,IBLLOC),  PAP_BLOCK(:,:,IBLLOC),      PAPH_BLOCK(:,:,IBLLOC),&
        & PLSM_BLOCK(:,IBLLOC),       LDCUM_BLOCK(:,IBLLOC),      KTYPE_BLOCK(:,IBLLOC), &
        & PLU_BLOCK(:,:,IBLLOC),      PLUDE_BLOCK(:,:,IBLLOC),    PSNDE_BLOCK(:,:,IBLLOC),    PMFU_BLOCK(:,:,IBLLOC),     PMFD_BLOCK(:,:,IBLLOC),&
        !---prognostic fields
        & PA_BLOCK(:,:,IBLLOC),       PCLV_BLOCK(:,:,:,IBLLOC),   PSUPSAT_BLOCK(:,:,IBLLOC),&
        !-- arrays for aerosol-cloud interactions
        & PLCRIT_AER_BLOCK(:,:,IBLLOC),PICRIT_AER_BLOCK(:,:,IBLLOC),&
        & PRE_ICE_BLOCK(:,:,IBLLOC),&
        & PCCN_BLOCK(:,:,IBLLOC),     PNICE_BLOCK(:,:,IBLLOC),&
        !---diagnostic output
        & PCOVPTOT_BLOCK(:,:,IBLLOC), PRAINFRAC_TOPRFZ_BLOCK(:,IBLLOC),&
        !---resulting fluxes
        & PFSQLF_BLOCK(:,:,IBLLOC),   PFSQIF_BLOCK(:,:,IBLLOC),  PFCQNNG_BLOCK(:,:,IBLLOC),  PFCQLNG_BLOCK(:,:,IBLLOC),&
        & PFSQRF_BLOCK(:,:,IBLLOC),   PFSQSF_BLOCK(:,:,IBLLOC),  PFCQRNG_BLOCK(:,:,IBLLOC),  PFCQSNG_BLOCK(:,:,IBLLOC),&
        & PFSQLTUR_BLOCK(:,:,IBLLOC), PFSQITUR_BLOCK(:,:,IBLLOC), &
        & PFPLSL_BLOCK(:,:,IBLLOC),   PFPLSN_BLOCK(:,:,IBLLOC),   PFHPSL_BLOCK(:,:,IBLLOC),   PFHPSN_BLOCK(:,:,IBLLOC),&
        & YRECLDP=LOCAL_YRECLDP, JL=JL)
      ENDDO
    ENDDO
    !$acc end parallel loop

    CALL TIMER%THREAD_ACC_END(TID)

    ! data to host
   !$acc host_data &
      !$acc use_device( &
      !$acc buffer_loc_block, &
      !$acc plude_block, &
      !$acc pcovptot_block, &
      !$acc prainfrac_toprfz_block, &
      !$acc pfsqlf_block, &
      !$acc pfsqif_block, &
      !$acc pfcqnng_block, &
      !$acc pfcqlng_block, &
      !$acc pfsqrf_block, &
      !$acc pfsqsf_block, &
      !$acc pfcqrng_block, &
      !$acc pfcqsng_block, &
      !$acc pfsqltur_block, &
      !$acc pfsqitur_block, &
      !$acc pfplsl_block, &
      !$acc pfplsn_block, &
      !$acc pfhpsl_block, &
      !$acc pfhpsn_block)
    ! copy
      call acc_memcpy_from_device(buffer_loc(:,:,:, BLOCK_START:BLOCK_END), buffer_loc_block, SIZEOF(buffer_loc(:,:,:,BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(plude(:, :, BLOCK_START:BLOCK_END), plude_block, SIZEOF(plude(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pcovptot(:, :, BLOCK_START:BLOCK_END), pcovptot_block, SIZEOF(pcovptot(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(prainfrac_toprfz(:, BLOCK_START:BLOCK_END), prainfrac_toprfz_block, SIZEOF(prainfrac_toprfz(:, BLOCK_START:BLOCK_END)))
    !copyout
      call acc_memcpy_from_device(pfsqlf(:, :, BLOCK_START:BLOCK_END), pfsqlf_block, SIZEOF(pfsqlf(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfsqif(:, :, BLOCK_START:BLOCK_END), pfsqif_block, SIZEOF(pfsqif(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfcqnng(:, :, BLOCK_START:BLOCK_END), pfcqnng_block, SIZEOF(pfcqnng(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfcqlng(:, :, BLOCK_START:BLOCK_END), pfcqlng_block, SIZEOF(pfcqlng(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfsqrf(:, :, BLOCK_START:BLOCK_END), pfsqrf_block, SIZEOF(pfsqrf(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfsqsf(:, :, BLOCK_START:BLOCK_END), pfsqsf_block, SIZEOF(pfsqsf(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfcqrng(:, :, BLOCK_START:BLOCK_END), pfcqrng_block, SIZEOF(pfcqrng(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfcqsng(:, :, BLOCK_START:BLOCK_END), pfcqsng_block, SIZEOF(pfcqsng(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfsqltur(:, :, BLOCK_START:BLOCK_END), pfsqltur_block, SIZEOF(pfsqltur(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfsqitur(:, :, BLOCK_START:BLOCK_END), pfsqitur_block, SIZEOF(pfsqitur(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfplsl(:, :, BLOCK_START:BLOCK_END), pfplsl_block, SIZEOF(pfplsl(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfplsn(:, :, BLOCK_START:BLOCK_END), pfplsn_block, SIZEOF(pfplsn(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfhpsl(:, :, BLOCK_START:BLOCK_END), pfhpsl_block, SIZEOF(pfhpsl(:, :, BLOCK_START:BLOCK_END)))
      call acc_memcpy_from_device(pfhpsn(:, :, BLOCK_START:BLOCK_END), pfhpsn_block, SIZEOF(pfhpsn(:, :, BLOCK_START:BLOCK_END)))
    !$acc end host_data
!
  
  ENDDO ! end of outer block loop


   ! buffer allocations
   !copyin
    DEALLOCATE(pt_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pq_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(buffer_tmp_block(NPROMA,NLEV,3+NCLV,BLOCK_BUFFER_SIZE))
    ! DEALLOCATE(BUFFER_CML_block(NPROMA,NLEV,3+NCLV,BLOCK_BUFFER_SIZE))
    DEALLOCATE(pvfa_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pvfl_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pvfi_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pdyna_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pdynl_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pdyni_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(phrsw_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(phrlw_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pvervel_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pap_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(paph_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(plsm_block(NPROMA, BLOCK_BUFFER_SIZE))
    DEALLOCATE(ldcum_block(NPROMA, BLOCK_BUFFER_SIZE))
    DEALLOCATE(ktype_block(NPROMA, BLOCK_BUFFER_SIZE))
    DEALLOCATE(plu_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(psnde_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pmfu_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pmfd_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pa_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pclv_block(NPROMA, NLEV, NCLV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(psupsat_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(plcrit_aer_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(picrit_aer_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pre_ice_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pccn_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pnice_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))

    ! copy
    DEALLOCATE(buffer_loc_block(NPROMA,NLEV,3+NCLV,BLOCK_BUFFER_SIZE))
    DEALLOCATE(plude_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pcovptot_block(NPROMA, NLEV, BLOCK_BUFFER_SIZE))
    DEALLOCATE(prainfrac_toprfz_block(NPROMA, BLOCK_BUFFER_SIZE))

    !copyout
    DEALLOCATE(pfsqlf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfsqif_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfcqnng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfcqlng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfsqrf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfsqsf_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfcqrng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfcqsng_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfsqltur_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfsqitur_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfplsl_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfplsn_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfhpsl_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))
    DEALLOCATE(pfhpsn_block(NPROMA, NLEV+1, BLOCK_BUFFER_SIZE))


    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_DBLK_K_CACHING

END MODULE CLOUDSC_DRIVER_GPU_SCC_DBLK_K_CACHING_MOD
