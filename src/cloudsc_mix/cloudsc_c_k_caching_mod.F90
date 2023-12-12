module cloudsc_c_k_caching_mod 

  use iso_fortran_env
  use iso_c_binding
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOECLDP, ONLY: TECLDP 
  implicit none

  public :: cloudsc_c_wrapper

  interface 

    subroutine cloudsc_c_launch(numcols, nproma, kidia, kfdia, klon, ptsphy, pt, &
          & pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
          & tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
          & pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, & 
          & pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
          & pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
          & pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn, yrecldp, &
          & ngpblks, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, r3les, r3ies, r4les, &
          & r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, & 
          & rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2) bind(c, name='cloudsc_c_launch')
   
        use iso_c_binding
        USE YOECLDP, ONLY: TECLDP
        implicit none

        type(c_ptr), value :: pt, &
          & pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
          & tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
          & pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, &
          & pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
          & pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
          & pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn
        integer(c_int), value :: numcols, nproma, kidia, kfdia, klon, ngpblks
        real(c_double), value :: ptsphy, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, & 
          & r3les, r3ies, r4les, r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, & 
          & ralfdcp, rtwat, rtice, rticecu, rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2
        type(TECLDP) :: yrecldp
        ! type(c_ptr) :: yrecldp
   end subroutine cloudsc_c_launch
  
  end interface

  contains

    subroutine cloudsc_c_wrapper(numcols, nproma, kidia, kfdia, klon, klev, ptsphy, pt, &
          & pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
          & tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
          & pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, &
          & pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
          & pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
          & pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn, yrecldp, &
          & ngpblks, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, r3les, r3ies, r4les, &
          & r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
          & rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)
      
  use iso_c_binding

        REAL(KIND=JPRB), INTENT(IN) :: PLCRIT_AER(:, :, :)
        REAL(KIND=JPRB), INTENT(IN) :: PICRIT_AER(:, :, :)
        REAL(KIND=JPRB), INTENT(IN) :: PRE_ICE(:, :, :)
        REAL(KIND=JPRB), INTENT(IN) :: PCCN(:, :, :)    ! liquid cloud condensation nuclei
        REAL(KIND=JPRB), INTENT(IN) :: PNICE(:, :, :)    ! ice number concentration (cf. CCN)

        INTEGER(KIND=JPIM), VALUE, INTENT(IN) :: KLON    ! Number of grid points
        INTEGER(KIND=JPIM), VALUE, INTENT(IN) :: KLEV    ! Number of levels
        INTEGER(KIND=JPIM), VALUE, INTENT(IN) :: KIDIA
        INTEGER(KIND=JPIM), VALUE, INTENT(IN) :: KFDIA
        INTEGER(KIND=JPIM), VALUE, INTENT(IN) :: ngpblks, nproma, numcols
        REAL(KIND=JPRB), VALUE, INTENT(IN) :: PTSPHY    ! Physics timestep
        REAL(KIND=JPRB), INTENT(IN) :: PT(:, :, :)    ! T at start of callpar
        REAL(KIND=JPRB), INTENT(IN) :: PQ(:, :, :)    ! Q at start of callpar
        ! REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP(:, :, :, :) ! TODO: tendencies split up?!
        ! REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC(:, :, :, :) ! TODO: tendencies split up?!
        REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_T(:, :, :) 
        REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_Q(:, :, :)   
        REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_A(:, :, :)
        REAL(KIND=JPRB), INTENT(IN) :: TENDENCY_TMP_CLD(:, :, :, :)
        REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_T(:, :, :)
        REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_Q(:, :, :)
        REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_A(:, :, :)
        REAL(KIND=JPRB), INTENT(INOUT) :: TENDENCY_LOC_CLD(:, :, :, :)
        REAL(KIND=JPRB), INTENT(IN) :: PVFA(:, :, :)    ! CC from VDF scheme
        REAL(KIND=JPRB), INTENT(IN) :: PVFL(:, :, :)    ! Liq from VDF scheme
        REAL(KIND=JPRB), INTENT(IN) :: PVFI(:, :, :)    ! Ice from VDF scheme
        REAL(KIND=JPRB), INTENT(IN) :: PDYNA(:, :, :)    ! CC from Dynamics
        REAL(KIND=JPRB), INTENT(IN) :: PDYNL(:, :, :)    ! Liq from Dynamics
        REAL(KIND=JPRB), INTENT(IN) :: PDYNI(:, :, :)    ! Liq from Dynamics
        REAL(KIND=JPRB), INTENT(IN) :: PHRSW(:, :, :)    ! Short-wave heating rate
        REAL(KIND=JPRB), INTENT(IN) :: PHRLW(:, :, :)    ! Long-wave heating rate
        REAL(KIND=JPRB), INTENT(IN) :: PVERVEL(:, :, :)    !Vertical velocity
        REAL(KIND=JPRB), INTENT(IN) :: PAP(:, :, :)    ! Pressure on full levels
        REAL(KIND=JPRB), INTENT(IN) :: PAPH(:, :, :)    ! Pressure on half levels
        REAL(KIND=JPRB), INTENT(IN) :: PLSM(:, :)    ! Land fraction (0-1)
        ! LOGICAL, INTENT(IN) :: LDCUM(:, :)    ! Convection active
        INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(:, :)    ! Convection type 0,1,2
        REAL(KIND=JPRB), INTENT(IN) :: PLU(:, :, :)    ! Conv. condensate
        REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(:, :, :)    ! Conv. detrained water
        REAL(KIND=JPRB), INTENT(IN) :: PSNDE(:, :, :)    ! Conv. detrained snow
        REAL(KIND=JPRB), INTENT(IN) :: PMFU(:, :, :)    ! Conv. mass flux up
        REAL(KIND=JPRB), INTENT(IN) :: PMFD(:, :, :)    ! Conv. mass flux down
        REAL(KIND=JPRB), INTENT(IN) :: PA(:, :, :)    ! Original Cloud fraction (t)

        REAL(KIND=JPRB), INTENT(IN) :: PCLV(:, :, :, :)

        ! Supersat clipped at previous time level in SLTEND
        REAL(KIND=JPRB), INTENT(IN) :: PSUPSAT(:, :, :)
        REAL(KIND=JPRB), INTENT(OUT) :: PCOVPTOT(:, :, :)    ! Precip fraction
        REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(:, :)
        ! Flux diagnostics for DDH budget
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(:, :, :)    ! Flux of liquid
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(:, :, :)    ! Flux of ice
        REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(:, :, :)    ! -ve corr for liq
        REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(:, :, :)    ! -ve corr for ice
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(:, :, :)    ! Flux diagnostics
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(:, :, :)    !    for DDH, generic
        REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(:, :, :)    ! rain
        REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(:, :, :)    ! snow
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(:, :, :)    ! liquid flux due to VDF
        REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(:, :, :)    ! ice flux due to VDF
        REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(:, :, :)    ! liq+rain sedim flux
        REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(:, :, :)    ! ice+snow sedim flux
        REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(:, :, :)    ! Enthalpy flux for liq
        REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(:, :, :)    ! Enthalp flux for ice

        REAL(KIND=JPRB) :: rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, &
          & r3les, r3ies, r4les, r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, &
          & ralfdcp, rtwat, rtice, rticecu, rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2
        type(TECLDP) :: yrecldp
        ! type(c_ptr) :: yrecldp

#if GPU_OFFLOAD == OMP_OFFLOAD
! NOTE:
!  seems like, 
!  omp target data use_device_PTR works for both AMD and NVIDIA machines
!  but,
!  omp target data use_device_ADDR works for AMD but NOT for NVIDIA machines

        !$omp target data use_device_ptr(pt, &
         !$omp& pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
         !$omp& tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
         !$omp& pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, &
         !$omp& pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
         !$omp& pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
         !$omp& pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn)
#elif GPU_OFFLOAD == ACC_OFFLOAD
        !$acc host_data use_device(pt, &
         !$acc& pq, tendency_tmp_t, tendency_tmp_q, tendency_tmp_a, tendency_tmp_cld, tendency_loc_t, &
         !$acc& tendency_loc_q, tendency_loc_a, tendency_loc_cld, pvfa, pvfl, pvfi, pdyna, &
         !$acc& pdynl, pdyni, phrsw, phrlw, pvervel, pap, paph, plsm, ktype, plu, plude, psnde, &
         !$acc& pmfu, pmfd, pa, pclv, psupsat, plcrit_aer, picrit_aer, pre_ice, pccn, pnice, &
         !$acc& pcovptot, prainfrac_toprfz, pfsqlf, pfsqif, pfcqnng, pfcqlng, pfsqrf, pfsqsf, &
         !$acc& pfcqrng, pfcqsng, pfsqltur, pfsqitur, pfplsl, pfplsn, pfhpsl, pfhpsn)
#endif

        CALL cloudsc_c_launch(numcols, nproma, kidia, kfdia, klon, ptsphy, & 
              & c_loc(pt), c_loc(pq), c_loc(tendency_tmp_t), c_loc(tendency_tmp_q), & 
              & c_loc(tendency_tmp_a), c_loc(tendency_tmp_cld), c_loc(tendency_loc_t), &
              & c_loc(tendency_loc_q), c_loc(tendency_loc_a), c_loc(tendency_loc_cld), & 
              & c_loc(pvfa), c_loc(pvfl), c_loc(pvfi), c_loc(pdyna), c_loc(pdynl), & 
              & c_loc(pdyni), c_loc(phrsw), c_loc(phrlw), c_loc(pvervel), c_loc(pap), &
              & c_loc(paph), c_loc(plsm), c_loc(ktype), c_loc(plu), c_loc(plude), &
              & c_loc(psnde), c_loc(pmfu), c_loc(pmfd), c_loc(pa), c_loc(pclv), & 
              & c_loc(psupsat), c_loc(plcrit_aer), c_loc(picrit_aer), c_loc(pre_ice), & 
              & c_loc(pccn), c_loc(pnice), c_loc(pcovptot), c_loc(prainfrac_toprfz), & 
              & c_loc(pfsqlf), c_loc(pfsqif), c_loc(pfcqnng), c_loc(pfcqlng), c_loc(pfsqrf), & 
              & c_loc(pfsqsf), c_loc(pfcqrng), c_loc(pfcqsng), c_loc(pfsqltur), & 
              & c_loc(pfsqitur), c_loc(pfplsl), c_loc(pfplsn), c_loc(pfhpsl), c_loc(pfhpsn), yrecldp, & 
              & ngpblks, rg, rd, rcpd, retv, rlvtt, rlstt, rlmlt, rtt, rv, r2es, r3les, r3ies, r4les, &
              & r4ies, r5les, r5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, rtwat, rtice, rticecu, &
              & rtwat_rtice_r, rtwat_rticecu_r, rkoop1, rkoop2)
#if GPU_OFFLOAD == OMP_OFFLOAD
        !$omp end target data 
#elif GPU_OFFLOAD == ACC_OFFLOAD
        !$acc end host_data
#endif

    end subroutine cloudsc_c_wrapper

end module cloudsc_c_k_caching_mod
