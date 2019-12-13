MODULE CLOUDSC_DRIVER_MOD

  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV
  USE TIMER_MOD, ONLY : FTIMER
  USE EC_PMON_MOD, ONLY: EC_PMON

#ifdef _OPENMP
use omp_lib
#else
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

IMPLICIT NONE


CONTAINS
  SUBROUTINE CLOUDSC_DRIVER( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, KFLDX, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & LDSLPHY,  LDMAINCALL, PA, &
     & PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
     & PEXTRA    )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC kernel

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NLEV, NGPTOT
    INTEGER(KIND=JPIM)                                    :: KFLDX 
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(IN)    :: TENDENCY_TMP(:) ! cumulative tendency used as input
    TYPE(STATE_TYPE)  ,POINTER,             INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFA(:,:,:)  ! CC from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFL(:,:,:)  ! Liq from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVFI(:,:,:)  ! Ice from VDF scheme
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNA(:,:,:) ! CC from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNL(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PDYNI(:,:,:) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PHRSW(:,:,:) ! Short-wave heating rate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PHRLW(:,:,:) ! Long-wave heating rate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PVERVEL(:,:,:) !Vertical velocity
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLSM(:,:)    ! Land fraction (0-1) 
    LOGICAL           ,POINTER, CONTIGUOUS, INTENT(IN)    :: LDCUM(:,:)   ! Convection active
    INTEGER(KIND=JPIM),POINTER, CONTIGUOUS, INTENT(IN)    :: KTYPE(:,:)   ! Convection type 0,1,2
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PSNDE(:,:,:) ! Conv. detrained snow
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    LOGICAL                                               :: LDSLPHY 
    LOGICAL                                               :: LDMAINCALL   ! T if main call to cloudsc
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(INOUT) :: PEXTRA(:,:,:,:) ! extra fields
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PLCRIT_AER(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PICRIT_AER(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PRE_ICE(:,:,:) 
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PCCN(:,:,:)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(IN)    :: PNICE(:,:,:)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PRAINFRAC_TOPRFZ(:,:) 
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQLF(:,:,:)  ! Flux of liquid
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQIF(:,:,:)  ! Flux of ice
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQLNG(:,:,:) ! -ve corr for liq
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQNNG(:,:,:) ! -ve corr for ice
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQRF(:,:,:)  ! Flux diagnostics
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQSF(:,:,:)  !    for DDH, generic
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQRNG(:,:,:) ! rain
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFCQSNG(:,:,:) ! snow
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQLTUR(:,:,:) ! liquid flux due to VDF
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFSQITUR(:,:,:) ! ice flux due to VDF
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB)   ,POINTER, CONTIGUOUS, INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

    REAL(KIND=JPRD) :: t1, t2, tloc, tdiff
    REAL(KIND=JPRD), parameter :: zhpm = 12482329.0_JPRD ! IBM P7 HPM flop count for 100 points at L137
    REAL(KIND=JPRD) :: zmflops ! MFlops/s rate
    REAL(KIND=JPRD) :: zfrac ! fraction of gp columns handled by thread

    INTEGER(KIND=JPIM) :: tid ! thread id from 0 .. NUMOMP - 1
    INTEGER(KIND=JPIM) :: coreid ! core id thread belongs to
    INTEGER(KIND=JPIM) :: icalls ! number of calls to CLOUDSC == number of blocks handled by particular thread
    INTEGER(KIND=JPIM) :: igpc ! number of gp columns handled by particular thread
    REAL(KIND=JPRD) :: zinfo(4,0:NUMOMP - 1)

    INTEGER(KIND=JPIB) :: ENERGY, POWER, POWER_TOTAL, POWER_MAX, POWER_COUNT
    LOGICAL            :: LEC_PMON = .FALSE.
    CHARACTER(LEN=1)   :: CLEC_PMON

#include "mycpu.intfb.h"

    CALL GET_ENVIRONMENT_VARIABLE('EC_PMON', CLEC_PMON)
    IF (CLEC_PMON == '1') LEC_PMON = .TRUE.

    POWER_MAX = 0_JPIB
    POWER_TOTAL = 0_JPIB
    POWER_COUNT = 0_JPIB

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMOMP=',i0,', NGPTOT=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    write(0,1003) NUMOMP,NGPTOT,NPROMA,NGPBLKS

    t1 = ftimer()

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,tloc,tid,coreid,icalls,igpc,energy,power) &
    !$omp& num_threads(NUMOMP)
    tloc = ftimer()
    tid = omp_get_thread_num()
    coreid = mycpu()
    icalls = 0
    igpc = 0

    !$omp do schedule(runtime) reduction(+:power_total,power_count) reduction(max:power_max)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         CALL CLOUDSC &
              & (    1,    ICEND,    NPROMA,  NLEV,&
              & PTSPHY,&
              & PT(:,:,IBL), PQ(:,:,IBL), TENDENCY_CML(IBL), TENDENCY_TMP(IBL), TENDENCY_LOC(IBL), &
              & PVFA(:,:,IBL), PVFL(:,:,IBL), PVFI(:,:,IBL), PDYNA(:,:,IBL), PDYNL(:,:,IBL), PDYNI(:,:,IBL), &
              & PHRSW(:,:,IBL),    PHRLW(:,:,IBL),&
              & PVERVEL(:,:,IBL),  PAP(:,:,IBL),      PAPH(:,:,IBL),&
              & PLSM(:,IBL),       LDCUM(:,IBL),      KTYPE(:,IBL), &
              & PLU(:,:,IBL),      PLUDE(:,:,IBL),    PSNDE(:,:,IBL),    PMFU(:,:,IBL),     PMFD(:,:,IBL),&
              & LDSLPHY,  LDMAINCALL, &
              !---prognostic fields
              & PA(:,:,IBL),&
              & PCLV(:,:,:,IBL),  &
              & PSUPSAT(:,:,IBL),&
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
              & PEXTRA(:,:,:,IBL),   KFLDX)

         IF (LEC_PMON) THEN
           ! Sample power consuption
           IF (MOD(IBL, 100) == 0) THEN
             CALL EC_PMON(ENERGY, POWER)
             POWER_MAX = MAX(POWER_MAX, POWER)
             POWER_TOTAL = POWER_TOTAL + POWER
             POWER_COUNT = POWER_COUNT + 1
           END IF
         END IF

         icalls = icalls + 1
         igpc = igpc + ICEND
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      tloc = ftimer() - tloc
      zinfo(1,tid) = tloc
      zinfo(2,tid) = coreid
      zinfo(3,tid) = icalls
      zinfo(4,tid) = igpc
      !$omp end parallel

      t2 = ftimer()

1000  format(1x,5a10,1x,a4,' : ',2a10)
1001  format(1x,5i10,1x,i4,' : ',2i10,:,' @ core#',i0)
1002  format(1x,5i10,1x,i4,' : ',2i10,  ' : TOTAL')
      write(0,1000) 'NUMOMP','NGPTOT','#GP-cols','#BLKS','NPROMA','tid#','Time(msec)','MFlops/s'
      do tid=0,NUMOMP-1
         tloc = zinfo(1,tid)
         coreid = int(zinfo(2,tid))
         icalls = int(zinfo(3,tid))
         igpc = int(zinfo(4,tid))
         zfrac = real(igpc,JPRB)/real(NGPTOT,JPRB)
         if (tloc > 0.0_JPRB) then
            zmflops = 1.0e-06_JPRB * zfrac * zhpm * (real(NGPTOT,JPRB)/real(100,JPRB))/tloc
         else
            zmflops = 0.0_JPRB
         endif
         write(0,1001) numomp,ngptot,igpc,icalls,nproma,tid,&
              & int(tloc*1000.0_JPRB),int(zmflops),coreid
      enddo
      tdiff = t2-t1
      zfrac = 1.0_JPRB
      if (tdiff > 0.0_JPRB) then
         zmflops = 1.0e-06_JPRB * zfrac * zhpm * (real(NGPTOT,JPRB)/real(100,JPRB))/tdiff
      else
         zmflops = 0.0_JPRB
      endif
      write(0,1002) numomp,ngptot,int(sum(zinfo(4,:))),ngpblks,nproma,-1,&
           & int(tdiff*1000.0_JPRB),int(zmflops)

      IF (LEC_PMON) THEN
        print *, "Power usage (sampled):: max: ", POWER_MAX, "avg:", &
         & (REAL(POWER_TOTAL, KIND=JPRD) / REAL(POWER_COUNT, KIND=JPRD)), &
         & "count:", POWER_COUNT
      END IF
    
  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
