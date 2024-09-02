! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_GPU_SCC_DBLK_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE CLOUDSC_GPU_SCC_MOD, ONLY: CLOUDSC_SCC
  USE FIELD_MODULE, ONLY: GET_DEVICE_DATA, ENSURE_HOST


  IMPLICIT NONE
  

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_DBLK( &
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

    ! Device pointers for kernel data for FIELD API device offload
    ! copyin
    CLASS(FIELD_3RB), POINTER :: pt_field_block
    CLASS(FIELD_3RB), POINTER :: pq_field_block
    CLASS(FIELD_4RB), POINTER :: buffer_tmp_field_block
    CLASS(FIELD_3RB), POINTER :: pvfa_field_block
    CLASS(FIELD_3RB), POINTER :: pvfl_field_block
    CLASS(FIELD_3RB), POINTER :: pvfi_field_block
    CLASS(FIELD_3RB), POINTER :: pdyna_field_block
    CLASS(FIELD_3RB), POINTER :: pdynl_field_block
    CLASS(FIELD_3RB), POINTER :: pdyni_field_block
    CLASS(FIELD_3RB), POINTER :: phrsw_field_block
    CLASS(FIELD_3RB), POINTER :: phrlw_field_block
    CLASS(FIELD_3RB), POINTER :: pvervel_field_block
    CLASS(FIELD_3RB), POINTER :: pap_field_block
    CLASS(FIELD_3RB), POINTER :: paph_field_block
    CLASS(FIELD_2RB), POINTER :: plsm_field_block
    CLASS(FIELD_2LM), POINTER :: ldcum_field_block
    CLASS(FIELD_2IM), POINTER :: ktype_field_block
    CLASS(FIELD_3RB), POINTER :: plu_field_block
    CLASS(FIELD_3RB), POINTER :: psnde_field_block
    CLASS(FIELD_3RB), POINTER :: pmfu_field_block
    CLASS(FIELD_3RB), POINTER :: pmfd_field_block
    CLASS(FIELD_3RB), POINTER :: pa_field_block
    CLASS(FIELD_4RB), POINTER :: pclv_field_block
    CLASS(FIELD_3RB), POINTER :: psupsat_field_block
    CLASS(FIELD_3RB), POINTER :: plcrit_aer_field_block
    CLASS(FIELD_3RB), POINTER :: picrit_aer_field_block
    CLASS(FIELD_3RB), POINTER :: pre_ice_field_block
    CLASS(FIELD_3RB), POINTER :: pccn_field_block
    CLASS(FIELD_3RB), POINTER :: pnice_field_block
    ! copy
    CLASS(FIELD_4RB), POINTER :: buffer_loc_field_block
    CLASS(FIELD_3RB), POINTER :: plude_field_block
    CLASS(FIELD_3RB), POINTER :: pcovptot_field_block
    CLASS(FIELD_2RB), POINTER :: prainfrac_toprfz_field_block
    ! copyout
    CLASS(FIELD_3RB), POINTER :: pfsqlf_field_block
    CLASS(FIELD_3RB), POINTER :: pfsqif_field_block
    CLASS(FIELD_3RB), POINTER :: pfcqnng_field_block
    CLASS(FIELD_3RB), POINTER :: pfcqlng_field_block
    CLASS(FIELD_3RB), POINTER :: pfsqrf_field_block
    CLASS(FIELD_3RB), POINTER :: pfsqsf_field_block
    CLASS(FIELD_3RB), POINTER :: pfcqrng_field_block
    CLASS(FIELD_3RB), POINTER :: pfcqsng_field_block
    CLASS(FIELD_3RB), POINTER :: pfsqltur_field_block
    CLASS(FIELD_3RB), POINTER :: pfsqitur_field_block
    CLASS(FIELD_3RB), POINTER :: pfplsl_field_block
    CLASS(FIELD_3RB), POINTER :: pfplsn_field_block
    CLASS(FIELD_3RB), POINTER :: pfhpsl_field_block
    CLASS(FIELD_3RB), POINTER :: pfhpsn_field_block

    ! Field API VIEWS for the DEVICE pr
    ! copyin
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pt_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pq_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: buffer_tmp_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pvfa_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pvfl_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pvfi_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pdyna_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pdynl_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pdyni_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: phrsw_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: phrlw_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pvervel_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pap_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: paph_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:) :: plsm_block
    LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:,:) :: ldcum_block
    INTEGER(KIND=JPIM), POINTER, CONTIGUOUS, DIMENSION(:,:) :: ktype_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: plu_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: psnde_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pmfu_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pmfd_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pa_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: pclv_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: psupsat_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: plcrit_aer_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: picrit_aer_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pre_ice_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pccn_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pnice_block
    ! copy
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: buffer_loc_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: plude_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pcovptot_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:) :: prainfrac_toprfz_block
    ! copyout
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqlf_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqif_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfcqnng_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfcqlng_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqrf_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqsf_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfcqrng_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfcqsng_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqltur_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfsqitur_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfplsl_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfplsn_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfhpsl_block
    REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: pfhpsn_block

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP
    
    ! double blocking variables
    INTEGER(KIND=JPIM) :: BLOCK_BUFFER_SIZE     ! block size for buffer of blocks in outer loop
    INTEGER(KIND=JPIM) :: BLOCK_SIZE            ! size of a block (BLOCK_END-BLOCK_START+1)
    INTEGER(KIND=JPIM) :: BLOCK_COUNT           ! number of blocks
    INTEGER(KIND=JPIM) :: BLOCK_IDX             ! idx of current block in [1,BLOCK_COUNT]
    INTEGER(KIND=JPIM) :: BLOCK_START           ! start of current block in [1,NGPBLKS]
    INTEGER(KIND=JPIM) :: BLOCK_END             ! end of current block in [1,NGPBLKS]
    INTEGER(KIND=JPIM) :: IBLLOC                ! local loop idx inside inner block loop
    INTEGER(KIND=JPIM) :: NUM_STREAMS           ! number of stream for async copy/execution
    INTEGER(KIND=JPIM) :: STREAM                ! number of stream for async copy/execution
    
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
    BLOCK_BUFFER_SIZE= MIN(200,NGPBLKS)
    BLOCK_COUNT=(NGPBLKS+BLOCK_BUFFER_SIZE-1)/BLOCK_BUFFER_SIZE
    ! Number of streams
    NUM_STREAMS = 5
    
    print *, 'BLOCK_BUFFER_SIZE=', BLOCK_BUFFER_SIZE
    print *, 'BLOCK_COUNT=', BLOCK_COUNT
    print *, 'NUM_STREAMS=', NUM_STREAMS
    


    DO BLOCK_IDX=0, BLOCK_COUNT-1
      BLOCK_START=BLOCK_IDX*BLOCK_BUFFER_SIZE+1
      BLOCK_END=MIN((BLOCK_IDX+1)*BLOCK_BUFFER_SIZE, NGPBLKS)
      BLOCK_SIZE = BLOCK_END-BLOCK_START+1
      
      CALL FIELD_NEW(F_PT, DATA=PT(:,:,lower:upper))

      CALL F_PT%GET_DEVICE_DATA_RDONLY(PT_BLOCK)
      
      STREAM = MODULO(BLOCK_IDX, NUM_STREAMS) + 1
    ! NEW FIELDS FOR EACH BLOCK
      !copyin
        CALL FIELD_NEW(pt_field_block, DATA=pt(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pq_field_block, DATA=pq(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(buffer_tmp_field_block, DATA=buffer_tmp(:,:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pvfa_field_block, DATA=pvfa(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pvfl_field_block, DATA=pvfl(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pvfi_field_block, DATA=pvfi(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pdyna_field_block, DATA=pdyna(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pdynl_field_block, DATA=pdynl(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pdyni_field_block, DATA=pdyni(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(phrsw_field_block, DATA=phrsw(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(phrlw_field_block, DATA=phrlw(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pvervel_field_block, DATA=pvervel(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pap_field_block, DATA=pap(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(paph_field_block, DATA=paph(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(plsm_field_block, DATA=plsm(:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(ldcum_field_block, DATA=ldcum(:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(ktype_field_block, DATA=ktype(:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(plu_field_block, DATA=plu(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(psnde_field_block, DATA=psnde(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pmfu_field_block, DATA=pmfu(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pmfd_field_block, DATA=pmfd(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pa_field_block, DATA=pa(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pclv_field_block, DATA=pclv(:,:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(psupsat_field_block, DATA=psupsat(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(plcrit_aer_field_block, DATA=plcrit_aer(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(picrit_aer_field_block, DATA=picrit_aer(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pre_ice_field_block, DATA=pre_ice(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pccn_field_block, DATA=pccn(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pnice_field_block, DATA=pnice(:,:,BLOCK_START:BLOCK_END))
        ! copy
        CALL FIELD_NEW(buffer_loc_field_block, DATA=buffer_loc(:,:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(plude_field_block, DATA=plude(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pcovptot_field_block, DATA=pcovptot(:,:,BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(prainfrac_toprfz_field_block, DATA=prainfrac_toprfz(:,BLOCK_START:BLOCK_END))
        ! copyout
        CALL FIELD_NEW(pfsqlf_field_block, DATA=pfsqlf(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfsqif_field_block, DATA=pfsqif(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfcqnng_field_block, DATA=pfcqnng(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfcqlng_field_block, DATA=pfcqlng(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfsqrf_field_block, DATA=pfsqrf(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfsqsf_field_block, DATA=pfsqsf(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfcqrng_field_block, DATA=pfcqrng(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfcqsng_field_block, DATA=pfcqsng(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfsqltur_field_block, DATA=pfsqltur(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfsqitur_field_block, DATA=pfsqitur(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfplsl_field_block, DATA=pfplsl(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfplsn_field_block, DATA=pfplsn(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfhpsl_field_block, DATA=pfhpsl(:, :, BLOCK_START:BLOCK_END))
        CALL FIELD_NEW(pfhpsn_field_block, DATA=pfhpsn(:, :, BLOCK_START:BLOCK_END))

    

    ! copy data to device
        !copyin
        CALL pt_field_block%GET_DEVICE_DATA_RDONLY(pt_block)
        CALL pq_field_block%GET_DEVICE_DATA_RDONLY(pq_block)
        CALL buffer_tmp_field_block%GET_DEVICE_DATA_RDONLY(buffer_tmp_block)
        CALL pvfa_field_block%GET_DEVICE_DATA_RDONLY(pvfa_block)
        CALL pvfl_field_block%GET_DEVICE_DATA_RDONLY(pvfl_block)
        CALL pvfi_field_block%GET_DEVICE_DATA_RDONLY(pvfi_block)
        CALL pdyna_field_block%GET_DEVICE_DATA_RDONLY(pdyna_block)
        CALL pdynl_field_block%GET_DEVICE_DATA_RDONLY(pdynl_block)
        CALL pdyni_field_block%GET_DEVICE_DATA_RDONLY(pdyni_block)
        CALL phrsw_field_block%GET_DEVICE_DATA_RDONLY(phrsw_block)
        CALL phrlw_field_block%GET_DEVICE_DATA_RDONLY(phrlw_block)
        CALL pvervel_field_block%GET_DEVICE_DATA_RDONLY(pvervel_block)
        CALL pap_field_block%GET_DEVICE_DATA_RDONLY(pap_block)
        CALL paph_field_block%GET_DEVICE_DATA_RDONLY(paph_block)
        CALL plsm_field_block%GET_DEVICE_DATA_RDONLY(plsm_block)
        CALL ldcum_field_block%GET_DEVICE_DATA_RDONLY(ldcum_block)
        CALL ktype_field_block%GET_DEVICE_DATA_RDONLY(ktype_block)
        CALL plu_field_block%GET_DEVICE_DATA_RDONLY(plu_block)
        CALL psnde_field_block%GET_DEVICE_DATA_RDONLY(psnde_block)
        CALL pmfu_field_block%GET_DEVICE_DATA_RDONLY(pmfu_block)
        CALL pmfd_field_block%GET_DEVICE_DATA_RDONLY(pmfd_block)
        CALL pa_field_block%GET_DEVICE_DATA_RDONLY(pa_block)
        CALL pclv_field_block%GET_DEVICE_DATA_RDONLY(pclv_block)
        CALL psupsat_field_block%GET_DEVICE_DATA_RDONLY(psupsat_block)
        CALL plcrit_aer_field_block%GET_DEVICE_DATA_RDONLY(plcrit_aer_block)
        CALL picrit_aer_field_block%GET_DEVICE_DATA_RDONLY(picrit_aer_block)
        CALL pre_ice_field_block%GET_DEVICE_DATA_RDONLY(pre_ice_block)
        CALL pccn_field_block%GET_DEVICE_DATA_RDONLY(pccn_block)
        CALL pnice_field_block%GET_DEVICE_DATA_RDONLY(pnice_block)
        ! copy
        CALL buffer_loc_field_block%GET_DEVICE_DATA_RDWR(buffer_loc_block)
        CALL plude_field_block%GET_DEVICE_DATA_RDWR(plude_block)
        CALL pcovptot_field_block%GET_DEVICE_DATA_RDWR(pcovptot_block)
        CALL prainfrac_toprfz_field_block%GET_DEVICE_DATA_RDWR(prainfrac_toprfz_block)
        ! copyout
        CALL pfsqlf_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqlf_block)
        CALL pfsqif_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqif_block)
        CALL pfcqnng_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfcqnng_block)
        CALL pfcqlng_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfcqlng_block)
        CALL pfsqrf_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqrf_block)
        CALL pfsqsf_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqsf_block)
        CALL pfcqrng_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfcqrng_block)
        CALL pfcqsng_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfcqsng_block)
        CALL pfsqltur_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqltur_block)
        CALL pfsqitur_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfsqitur_block)
        CALL pfplsl_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfplsl_block)
        CALL pfplsn_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfplsn_block)
        CALL pfhpsl_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfhpsl_block)
        CALL pfhpsn_field_block%GET_DEVICE_DATA_WRITE_ONLY(pfhpsn_block)



        !$acc data copyin(yrecldp), deviceptr(pt_block, &
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

        ! Local timer for each thread
        TID = GET_THREAD_NUM()
        CALL TIMER%THREAD_ACC_START(TID)
      
        
        !$acc parallel loop gang vector_length(NPROMA) copy(LOCAL_YRECLDP)
        DO IBLLOC=1, BLOCK_BUFFER_SIZE ! just a way to loop over NGPBLKS
            IBL= BLOCK_BUFFER_SIZE*BLOCK_IDX +IBLLOC
            JKGLO=(IBL-1)*NPROMA+1
            ICEND=MIN(NPROMA, NGPTOT-JKGLO+1)

             CALL CLOUDSC_SCC &
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
              & YRECLDP=LOCAL_YRECLDP)

          ENDDO
      !$acc end parallel loop
      !$acc and data
      
      CALL TIMER%THREAD_ACC_END(TID)
   
      CALL buffer_loc_field_block%SYNC_HOST_RWDR()
      CALL plude_field_block%SYNC_HOST_RWDR()
      CALL pcovptot_field_block%SYNC_HOST_RWDR()
      CALL prainfrac_toprfz_field_block%SYNC_HOST_RWDR()
      ! copyout
      CALL pfsqlf_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfsqif_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfcqnng_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfcqlng_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfsqrf_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfsqsf_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfcqrng_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfcqsng_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfsqltur_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfsqitur_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfplsl_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfplsn_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfhpsl_field_block%SYNC_HOST_WRITE_ONLY()
      CALL pfhpsn_field_block%SYNC_HOST_WRITE_ONLY()
  
    ENDDO ! end of outer block loop

    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)


  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC_DBLK

END MODULE CLOUDSC_DRIVER_GPU_SCC_DBLK_MOD
