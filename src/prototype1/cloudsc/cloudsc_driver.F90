! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE CLOUDSC_DRIVER
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMPHYDER, ONLY : STATE_TYPE, DIMENSION_TYPE
USE YOECLDP, ONLY : NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
USE EXPAND_MOD, ONLY : loadjj, expand
USE DIAG_MOD, ONLY : NGPTOT, NPROMA, NGPBLKS, NPROMAS_IN, NUMOMP
USE timer_mod, ONLY : ftimer
USE diff_mod, ONLY : errhead, errcalc, saveref
USE serialize_mod, ONLY : query_dimensions, serialize, deserialize, serialize_reference

#ifdef _OPENMP
use omp_lib
#else
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

implicit none

INTEGER(KIND=JPIM), parameter :: iu = 123

TYPE (STATE_TYPE)  :: tendency_cml   ! cumulative tendency used for final output
TYPE (STATE_TYPE)  :: tendency_tmp   ! cumulative tendency used as input
TYPE (STATE_TYPE)  :: tendency_loc

TYPE (STATE_TYPE), ALLOCATABLE  :: ztendency_cml(:)   ! cumulative tendency used for final output
TYPE (STATE_TYPE), ALLOCATABLE  :: ztendency_tmp(:)   ! cumulative tendency used as input
TYPE (STATE_TYPE), ALLOCATABLE  :: ztendency_loc(:)

REAL(KIND=JPRB)   ,ALLOCATABLE    :: PLCRIT_AER(:,:), zPLCRIT_AER(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PICRIT_AER(:,:), zPICRIT_AER(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PRE_ICE(:,:), zPRE_ICE(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PCCN(:,:), zPCCN(:,:,:)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PNICE(:,:), zPNICE(:,:,:)    ! ice number concentration (cf. CCN)
REAL(KIND=JPRB)                   :: PTSPHY            ! Physics timestep
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PT(:,:), zPT(:,:,:)    ! T at start of callpar
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PQ(:,:), zPQ(:,:,:)    ! Q at start of callpar
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PVFA(:,:), zPVFA(:,:,:)  ! CC from VDF scheme
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PVFL(:,:), zPVFL(:,:,:)  ! Liq from VDF scheme
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PVFI(:,:), zPVFI(:,:,:)  ! Ice from VDF scheme
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PDYNA(:,:), zPDYNA(:,:,:) ! CC from Dynamics
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PDYNL(:,:), zPDYNL(:,:,:) ! Liq from Dynamics
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PDYNI(:,:), zPDYNI(:,:,:) ! Liq from Dynamics
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PHRSW(:,:), zPHRSW(:,:,:) ! Short-wave heating rate
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PHRLW(:,:), zPHRLW(:,:,:) ! Long-wave heating rate
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PVERVEL(:,:), zPVERVEL(:,:,:) !Vertical velocity
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PAP(:,:), zPAP(:,:,:)   ! Pressure on full levels
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PAPH(:,:), zPAPH(:,:,:)! Pressure on half levels
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PLSM(:), zPLSM(:,:)       ! Land fraction (0-1)
LOGICAL           ,ALLOCATABLE    :: LDCUM(:), LLCUM(:,:)      ! Convection active
INTEGER(KIND=JPIM),ALLOCATABLE    :: KTYPE(:), ITYPE(:,:)      ! Convection type 0,1,2
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PLU(:,:), zPLU(:,:,:)   ! Conv. condensate
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PLUDE(:,:), zPLUDE(:,:,:) ! Conv. detrained water
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PLUDE_tmp(:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PSNDE(:,:), zPSNDE(:,:,:) ! Conv. detrained snow
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PMFU(:,:), zPMFU(:,:,:)  ! Conv. mass flux up
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PMFD(:,:), zPMFD(:,:,:)  ! Conv. mass flux down
LOGICAL                           :: LDSLPHY
LOGICAL                           :: LDMAINCALL       ! T if main call to cloudsc
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PA(:,:), zPA(:,:,:)    ! Original Cloud fraction (t)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PEXTRA(:,:,:), zPEXTRA(:,:,:,:) ! extra fields
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PEXTRA_tmp(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PCLV(:,:,:), zPCLV(:,:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: PSUPSAT(:,:), zPSUPSAT(:,:,:)

REAL(KIND=JPRB)   ,ALLOCATABLE   :: PCOVPTOT(:,:), zPCOVPTOT(:,:,:) ! Precip fraction
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PRAINFRAC_TOPRFZ(:), zPRAINFRAC_TOPRFZ(:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQLF(:,:), zPFSQLF(:,:,:)  ! Flux of liquid
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQIF(:,:), zPFSQIF(:,:,:)  ! Flux of ice
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFCQLNG(:,:), zPFCQLNG(:,:,:) ! -ve corr for liq
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFCQNNG(:,:), zPFCQNNG(:,:,:) ! -ve corr for ice
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQRF(:,:), zPFSQRF(:,:,:)  ! Flux diagnostics
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQSF(:,:), zPFSQSF(:,:,:)  !    for DDH, generic
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFCQRNG(:,:), zPFCQRNG(:,:,:) ! rain
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFCQSNG(:,:), zPFCQSNG(:,:,:) ! snow
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQLTUR(:,:), zPFSQLTUR(:,:,:) ! liquid flux due to VDF
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFSQITUR(:,:), zPFSQITUR(:,:,:) ! ice flux due to VDF
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFPLSL(:,:), zPFPLSL(:,:,:) ! liq+rain sedim flux
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFPLSN(:,:), zPFPLSN(:,:,:) ! ice+snow sedim flux
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFHPSL(:,:), zPFHPSL(:,:,:) ! Enthalpy flux for liq
REAL(KIND=JPRB)   ,ALLOCATABLE   :: PFHPSN(:,:), zPFHPSN(:,:,:) ! Enthalp flux for ice

INTEGER(KIND=JPIM) :: KLON,KLEV,KFLDX
INTEGER(KIND=JPIM) :: KIDIA, KFDIA
INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND

INTEGER(KIND=JPIM) :: JPR, ICOUNT
REAL(KIND=JPRB)   ,ALLOCATABLE, TARGET   :: PSTATE_ARRL(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE, TARGET   :: PSTATE_CLD(:,:,:,:)

TYPE SPACE_T
   REAL(KIND=JPRB), POINTER :: arrl(:,:,:)
   REAL(KIND=JPRB), POINTER :: cld(:,:,:,:)
END TYPE SPACE_T

TYPE(SPACE_T) :: zcml, ztmp, zloc

REAL(KIND=JPRB) :: t1, t2, tloc, tdiff

character(len=80) :: cl_omp_schedule

CHARACTER(LEN=1) :: write_input, write_reference

REAL(KIND=JPRB), parameter :: zhpm = 12482329.0_JPRB ! IBM P7 HPM flop count for 100 points at L137
REAL(KIND=JPRB) :: zmflops ! MFlops/s rate
REAL(KIND=JPRB) :: zthrput ! Throughput col/s
REAL(KIND=JPRB) :: zfrac ! fraction of gp columns handled by thread

INTEGER(KIND=JPIM) :: tid ! thread id from 0 .. NUMOMP - 1
INTEGER(KIND=JPIM) :: coreid ! core id thread belongs to
INTEGER(KIND=JPIM) :: icalls ! number of calls to CLOUDSC == number of blocks handled by particular thread
INTEGER(KIND=JPIM) :: igpc ! number of gp columns handled by particular thread

REAL(KIND=JPRB) :: zinfo(4,0:NUMOMP - 1)
! zinfo(1,tid) == thread (wall clock) time
! zinfo(2,tid) == coreid
! zinfo(3,tid) == number of calls CLOUDSC == number of block processed by this tid
! zinfo(4,tid) == igpc

#include "cloudsc_in.intfb.h"
#include "cloudsc.intfb.h"
#include "mycpu.intfb.h"

CALL GET_ENVIRONMENT_VARIABLE('CLOUDSC_WRITE_INPUT', write_input)
CALL GET_ENVIRONMENT_VARIABLE('CLOUDSC_WRITE_REFERENCE', write_reference)

open(iu,file='cloudsc.bin',status='old',&
     & access='stream', form='unformatted', convert='BIG_ENDIAN')

read(iu) KLON,KLEV,KFLDX
write(0,*) 'KLON,KLEV,KFLDX,NCLV=',KLON,KLEV,KFLDX,NCLV

ALLOCATE(PLCRIT_AER(KLON,KLEV))
ALLOCATE(PICRIT_AER(KLON,KLEV))
ALLOCATE(PRE_ICE(KLON,KLEV))
ALLOCATE(PCCN(KLON,KLEV))
ALLOCATE(PNICE(KLON,KLEV))
ALLOCATE(PT(KLON,KLEV))
ALLOCATE(PQ(KLON,KLEV))
ALLOCATE(PVFA(KLON,KLEV))
ALLOCATE(PVFL(KLON,KLEV))
ALLOCATE(PVFI(KLON,KLEV))
ALLOCATE(PDYNA(KLON,KLEV))
ALLOCATE(PDYNL(KLON,KLEV))
ALLOCATE(PDYNI(KLON,KLEV))
ALLOCATE(PHRSW(KLON,KLEV))
ALLOCATE(PHRLW(KLON,KLEV))
ALLOCATE(PVERVEL(KLON,KLEV))
ALLOCATE(PAP(KLON,KLEV))
ALLOCATE(PAPH(KLON,KLEV+1))
ALLOCATE(PLSM(KLON))
ALLOCATE(LDCUM(KLON))
ALLOCATE(KTYPE(KLON))
ALLOCATE(PLU(KLON,KLEV))
ALLOCATE(PLUDE(KLON,KLEV))
ALLOCATE(PSNDE(KLON,KLEV))
ALLOCATE(PMFU(KLON,KLEV))
ALLOCATE(PMFD(KLON,KLEV))
ALLOCATE(PA(KLON,KLEV))
ALLOCATE(PEXTRA(KLON,KLEV,KFLDX))
ALLOCATE(PCLV(KLON,KLEV,NCLV))
ALLOCATE(PSUPSAT(KLON,KLEV))

ICOUNT = 0
ALLOCATE(PSTATE_ARRL(KLON,KLEV,18))
ICOUNT = ICOUNT + 1
tendency_cml%u  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_cml%v  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_cml%o3 => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_cml%a  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_cml%q  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_cml%T  => PSTATE_ARRL(:,:,ICOUNT)

ICOUNT = ICOUNT + 1
tendency_tmp%u  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%v  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%o3 => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%a  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%q  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%T  => PSTATE_ARRL(:,:,ICOUNT)

ICOUNT = ICOUNT + 1
tendency_loc%u  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%v  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%o3 => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%a  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%q  => PSTATE_ARRL(:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%T  => PSTATE_ARRL(:,:,ICOUNT)

PSTATE_ARRL = 0.0_JPRB

ICOUNT = 0
ALLOCATE(PSTATE_CLD(KLON,KLEV,NCLV,3))
ICOUNT = ICOUNT + 1
tendency_cml%cld  => PSTATE_CLD(:,:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_tmp%cld  => PSTATE_CLD(:,:,:,ICOUNT)
ICOUNT = ICOUNT + 1
tendency_loc%cld  => PSTATE_CLD(:,:,:,ICOUNT)

PSTATE_CLD = 0.0_JPRB

CALL CLOUDSC_IN &
     !---input
     & (iu,    KLON,    KLEV,&
     & PTSPHY,&
     & PT, PQ, tendency_cml,tendency_tmp, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW,&
     & PVERVEL,  PAP,      PAPH,&
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD,&
     & LDSLPHY,  LDMAINCALL, &
     !---prognostic fields
     & PA,&
     & PCLV,  &
     & PSUPSAT, &
     & PLCRIT_AER,PICRIT_AER,&
     & PRE_ICE,&
     & PCCN,     PNICE,&
     & PEXTRA,   KFLDX)

close(iu)

! Serialize the reference input (for debug and migration purposes)
if (write_input == '1') then
  call serialize( &
   & KLON, KLEV, PTSPHY,&
   & PT, PQ, TENDENCY_CML, TENDENCY_TMP, &
   & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
   & PHRSW, PHRLW, PVERVEL, PAP, PAPH, &
   & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, &
   & LDSLPHY,  LDMAINCALL, &
   & PA, PCLV, PSUPSAT, PLCRIT_AER,PICRIT_AER, &
   & PRE_ICE, PCCN, PNICE, PEXTRA, KFLDX)
end if

call query_dimensions(KLON, KLEV, KFLDX, name='input')
write(0,*) 'KLON,KLEV,KFLDX,NCLV=',KLON,KLEV,KFLDX,NCLV

call deserialize( &
 & KLON, KLEV, PTSPHY,&
 & PT, PQ, TENDENCY_CML, TENDENCY_TMP, &
 & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
 & PHRSW, PHRLW, PVERVEL, PAP, PAPH, &
 & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, &
 & LDSLPHY, LDMAINCALL, &
 & PA, PCLV, PSUPSAT, PLCRIT_AER,PICRIT_AER, &
 & PRE_ICE, PCCN, PNICE, PEXTRA, KFLDX)


! Create reference results dataset based on input data & one sweep over CLOUDSC

! These guys are INTENT(INOUT) -- we do not want to destroy the original values
ALLOCATE(PLUDE_tmp(KLON,KLEV))
PLUDE_tmp = PLUDE
ALLOCATE(PEXTRA_tmp(KLON,KLEV,KFLDX))
PEXTRA_tmp = PEXTRA

! These are INTENT(OUT)
ALLOCATE(PCOVPTOT(KLON,KLEV))
PCOVPTOT = 0.0_JPRB ! this is not fully initialized in CLOUDSC
ALLOCATE(PRAINFRAC_TOPRFZ(KLON))
ALLOCATE(PFSQLF(KLON,KLEV+1))
ALLOCATE(PFSQIF(KLON,KLEV+1))
ALLOCATE(PFCQLNG(KLON,KLEV+1))
ALLOCATE(PFCQNNG(KLON,KLEV+1))
ALLOCATE(PFSQRF(KLON,KLEV+1))
ALLOCATE(PFSQSF(KLON,KLEV+1))
ALLOCATE(PFCQRNG(KLON,KLEV+1))
ALLOCATE(PFCQSNG(KLON,KLEV+1))
ALLOCATE(PFSQLTUR(KLON,KLEV+1))
ALLOCATE(PFSQITUR(KLON,KLEV+1))
ALLOCATE(PFPLSL(KLON,KLEV+1))
ALLOCATE(PFPLSN(KLON,KLEV+1))
ALLOCATE(PFHPSL(KLON,KLEV+1))
ALLOCATE(PFHPSN(KLON,KLEV+1))

KIDIA = 1
KFDIA = KLON

CALL CLOUDSC &
     & (KIDIA,    KFDIA,    KLON,    KLEV,&
     & PTSPHY,&
     & PT, PQ, tendency_cml,tendency_tmp,tendency_loc, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW,&
     & PVERVEL,  PAP,      PAPH,&
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE_tmp,    PSNDE,    PMFU,     PMFD,&
     & LDSLPHY,  LDMAINCALL, &
     !---prognostic fields
     & PA,&
     & PCLV,  &
     & PSUPSAT,&
     !-- arrays for aerosol-cloud interactions
     & PLCRIT_AER,PICRIT_AER,&
     & PRE_ICE,&
     & PCCN,     PNICE,&
     !---diagnostic output
     & PCOVPTOT, PRAINFRAC_TOPRFZ,&
     !---resulting fluxes
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,&
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG,&
     & PFSQLTUR, PFSQITUR , &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN,&
     & PEXTRA_tmp,   KFLDX)

! Generate reference data if flag is set
if (write_reference == '1') then
  call serialize_reference( KLON, KLEV, KFLDX, &
   & PLUDE_tmp,    PCOVPTOT, PRAINFRAC_TOPRFZ,&
   & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,&
   & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG,&
   & PFSQLTUR, PFSQITUR , &
   & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
   & TENDENCY_LOC)
end if

CALL saveref('PLUDE',PLUDE_tmp)
DEALLOCATE(PLUDE_tmp)
DEALLOCATE(PEXTRA_tmp)

CALL saveref('PCOVPTOT',PCOVPTOT)
DEALLOCATE(PCOVPTOT)
CALL saveref('PRAINFRAC_TOPRFZ',PRAINFRAC_TOPRFZ)
DEALLOCATE(PRAINFRAC_TOPRFZ)
CALL saveref('PFSQLF',PFSQLF)
DEALLOCATE(PFSQLF)
CALL saveref('PFSQIF',PFSQIF)
DEALLOCATE(PFSQIF)
CALL saveref('PFCQLNG',PFCQLNG)
DEALLOCATE(PFCQLNG)
CALL saveref('PFCQNNG',PFCQNNG)
DEALLOCATE(PFCQNNG)
CALL saveref('PFSQRF',PFSQRF)
DEALLOCATE(PFSQRF)
CALL saveref('PFSQSF',PFSQSF)
DEALLOCATE(PFSQSF)
CALL saveref('PFCQRNG',PFCQRNG)
DEALLOCATE(PFCQRNG)
CALL saveref('PFCQSNG',PFCQSNG)
DEALLOCATE(PFCQSNG)
CALL saveref('PFSQLTUR',PFSQLTUR)
DEALLOCATE(PFSQLTUR)
CALL saveref('PFSQITUR',PFSQITUR)
DEALLOCATE(PFSQITUR)
CALL saveref('PFPLSL',PFPLSL)
DEALLOCATE(PFPLSL)
CALL saveref('PFPLSN',PFPLSN)
DEALLOCATE(PFPLSN)
CALL saveref('PFHPSL',PFHPSL)
DEALLOCATE(PFHPSL)
CALL saveref('PFHPSN',PFHPSN)
DEALLOCATE(PFHPSN)

!CALL saveref('tendency_loc%u',tendency_loc%u)
!CALL saveref('tendency_loc%v',tendency_loc%v)
!CALL saveref('tendency_loc%o3',tendency_loc%o3)
CALL saveref('tendency_loc%a',tendency_loc%a)
CALL saveref('tendency_loc%q',tendency_loc%q)
CALL saveref('tendency_loc%T',tendency_loc%T)
CALL saveref('tendency_loc%cld',tendency_loc%cld)


if (NGPTOT >= 1 .and. allocated(NPROMAS_IN)) then
   CALL get_environment_variable('OMP_SCHEDULE', cl_omp_schedule)

   ! Sweeps over varios NPROMAs
   do JPR=1,size(NPROMAS_IN)
      NPROMA = NPROMAS_IN(JPR)
      if (NPROMA < 1) cycle

      call loadjj(KLON) ! sets NGPBLKS etc.

1003  format(5x,'=> CASE#',i0,' : NUMOMP=',i0,', NGPTOT=',i0,', NPROMA=',i0,', NGPBLKS=',i0,', OMP_SCHEDULE=',a)
      write(0,1003) JPR,NUMOMP,NGPTOT,NPROMA,NGPBLKS,trim(adjustl(cl_omp_schedule))

      ! Expand KLON grid point column input as if we had NGPTOT of them
      ! Create blocks for multi-threading
      ! e.g. from array PA(KLON,KLEV) we create zPA(NPROMA,KLEV,NGPBLKS)

      ALLOCATE(zPLCRIT_AER(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PLCRIT_AER,zPLCRIT_AER)
      ALLOCATE(zPICRIT_AER(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PICRIT_AER,zPICRIT_AER)
      ALLOCATE(zPRE_ICE(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PRE_ICE,zPRE_ICE)
      ALLOCATE(zPCCN(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PCCN,zPCCN)
      ALLOCATE(zPNICE(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PNICE,zPNICE)
      ALLOCATE(zPT(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PT,zPT)
      ALLOCATE(zPQ(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PQ,zPQ)
      ALLOCATE(zPVFA(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PVFA,zPVFA)
      ALLOCATE(zPVFL(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PVFL,zPVFL)
      ALLOCATE(zPVFI(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PVFI,zPVFI)
      ALLOCATE(zPDYNA(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PDYNA,zPDYNA)
      ALLOCATE(zPDYNL(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PDYNL,zPDYNL)
      ALLOCATE(zPDYNI(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PDYNI,zPDYNI)
      ALLOCATE(zPHRSW(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PHRSW,zPHRSW)
      ALLOCATE(zPHRLW(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PHRLW,zPHRLW)
      ALLOCATE(zPVERVEL(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PVERVEL,zPVERVEL)
      ALLOCATE(zPAP(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PAP,zPAP)
      ALLOCATE(zPAPH(NPROMA,KLEV+1,NGPBLKS))
      call expand(KLON,KLEV+1,PAPH,zPAPH)
      ALLOCATE(zPLSM(NPROMA,NGPBLKS))
      call expand(KLON,PLSM,zPLSM)
      ALLOCATE(LLCUM(NPROMA,NGPBLKS))
      call expand(KLON,LDCUM,LLCUM)
      ALLOCATE(ITYPE(NPROMA,NGPBLKS))
      call expand(KLON,KTYPE,ITYPE)
      ALLOCATE(zPLU(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PLU,zPLU)
      ALLOCATE(zPLUDE(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PLUDE,zPLUDE)
      ALLOCATE(zPSNDE(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PSNDE,zPSNDE)
      ALLOCATE(zPMFU(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PMFU,zPMFU)
      ALLOCATE(zPMFD(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PMFD,zPMFD)
      ALLOCATE(zPA(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PA,zPA)
      ALLOCATE(zPEXTRA(NPROMA,KLEV,max(1,KFLDX),NGPBLKS))
      call expand(KLON,KLEV,KFLDX,PEXTRA,zPEXTRA)
      ALLOCATE(zPCLV(NPROMA,KLEV,NCLV,NGPBLKS))
      call expand(KLON,KLEV,NCLV,PCLV,zPCLV)
      ALLOCATE(zPSUPSAT(NPROMA,KLEV,NGPBLKS))
      call expand(KLON,KLEV,PSUPSAT,zPSUPSAT)
      CALL ALLOCATE_STATE_ARRAY(zcml, ztendency_cml,NPROMA,KLEV,NCLV,NGPBLKS)
      call expand(KLON,KLEV,NCLV,tendency_cml,ztendency_cml)
      CALL ALLOCATE_STATE_ARRAY(ztmp, ztendency_tmp,NPROMA,KLEV,NCLV,NGPBLKS)
      call expand(KLON,KLEV,NCLV,tendency_tmp,ztendency_tmp)

      ! Intent OUT
      ALLOCATE(zPCOVPTOT(NPROMA,KLEV,NGPBLKS))
      ALLOCATE(zPRAINFRAC_TOPRFZ(NPROMA,NGPBLKS))
      ALLOCATE(zPFSQLF(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFSQIF(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFCQLNG(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFCQNNG(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFSQRF(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFSQSF(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFCQRNG(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFCQSNG(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFSQLTUR(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFSQITUR(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFPLSL(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFPLSN(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFHPSL(NPROMA,KLEV+1,NGPBLKS))
      ALLOCATE(zPFHPSN(NPROMA,KLEV+1,NGPBLKS))
      CALL ALLOCATE_STATE_ARRAY(zloc,ztendency_loc,NPROMA,KLEV,NCLV,NGPBLKS)

      t1 = ftimer()

      !$omp parallel default(shared) private(JKGLO,IBL,ICEND,tloc,tid,coreid,icalls,igpc) &
      !$omp& num_threads(NUMOMP)
      tloc = ftimer()
      tid = omp_get_thread_num()
      coreid = mycpu()
      icalls = 0
      igpc = 0
      !$omp do schedule(runtime)
      DO JKGLO=1,NGPTOT,NPROMA
         IBL=(JKGLO-1)/NPROMA+1
         ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         zPCOVPTOT(:,:,IBL) = 0.0_JPRB
         ztendency_loc(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         CALL CLOUDSC &
              & (    1,    ICEND,    NPROMA,  KLEV,&
              & PTSPHY,&
              & zPT(:,:,IBL), zPQ(:,:,IBL), ztendency_cml(IBL),ztendency_tmp(IBL),ztendency_loc(IBL), &
              & zPVFA(:,:,IBL), zPVFL(:,:,IBL), zPVFI(:,:,IBL), zPDYNA(:,:,IBL), zPDYNL(:,:,IBL), zPDYNI(:,:,IBL), &
              & zPHRSW(:,:,IBL),    zPHRLW(:,:,IBL),&
              & zPVERVEL(:,:,IBL),  zPAP(:,:,IBL),      zPAPH(:,:,IBL),&
              & zPLSM(:,IBL),     LLCUM(:,IBL),    ITYPE(:,IBL), &
              & zPLU(:,:,IBL),      zPLUDE(:,:,IBL),    zPSNDE(:,:,IBL),    zPMFU(:,:,IBL),     zPMFD(:,:,IBL),&
              & LDSLPHY,  LDMAINCALL, &
              !---prognostic fields
              & zPA(:,:,IBL),&
              & zPCLV(:,:,:,IBL),  &
              & zPSUPSAT(:,:,IBL),&
              !-- arrays for aerosol-cloud interactions
              & zPLCRIT_AER(:,:,IBL),zPICRIT_AER(:,:,IBL),&
              & zPRE_ICE(:,:,IBL),&
              & zPCCN(:,:,IBL),     zPNICE(:,:,IBL),&
              !---diagnostic output
              & zPCOVPTOT(:,:,IBL), zPRAINFRAC_TOPRFZ(:,IBL),&
              !---resulting fluxes
              & zPFSQLF(:,:,IBL),   zPFSQIF (:,:,IBL),  zPFCQNNG(:,:,IBL),  zPFCQLNG(:,:,IBL),&
              & zPFSQRF(:,:,IBL),   zPFSQSF (:,:,IBL),  zPFCQRNG(:,:,IBL),  zPFCQSNG(:,:,IBL),&
              & zPFSQLTUR(:,:,IBL), zPFSQITUR (:,:,IBL), &
              & zPFPLSL(:,:,IBL),   zPFPLSN(:,:,IBL),   zPFHPSL(:,:,IBL),   zPFHPSN(:,:,IBL),&
              & zPEXTRA(:,:,:,IBL),   KFLDX)

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

1000  format(1x,5a10,1x,a4,' : ',3a10)
1001  format(1x,5i10,1x,i4,' : ',3i10,:,' @ core#',i0)
1002  format(1x,5i10,1x,i4,' : ',3i10,  ' : TOTAL')
1005  format(1x,'Reference MFLOP count for 100 columns :',1x,f12.8)
      write(0,1005) 1.0E-06_JPRB * ZHPM
      write(0,1000) 'NUMOMP','NGPTOT','#GP-cols','#BLKS','NPROMA','tid#','Time(msec)','MFlops/s','col/s'
      do tid=0,NUMOMP-1
         tloc = zinfo(1,tid)
         coreid = int(zinfo(2,tid))
         icalls = int(zinfo(3,tid))
         igpc = int(zinfo(4,tid))
         zfrac = real(igpc,JPRB)/real(NGPTOT,JPRB)
         if (tloc > 0.0_JPRB) then
            zmflops = 1.0e-06_JPRB * zfrac * zhpm * (real(NGPTOT,JPRB)/real(100,JPRB))/tloc
            zthrput = real(NGPTOT,JPRB)/tloc
         else
            zmflops = 0.0_JPRB
            zthrput = 0.0_JPRB
         endif
         write(0,1001) numomp,ngptot,igpc,icalls,nproma,tid,&
              & int(tloc*1000.0_JPRB),int(zmflops),int(zthrput),coreid
      enddo
      tdiff = t2-t1
      zfrac = 1.0_JPRB
      if (tdiff > 0.0_JPRB) then
         zmflops = 1.0e-06_JPRB * zfrac * zhpm * (real(NGPTOT,JPRB)/real(100,JPRB))/tdiff
         zthrput = real(NGPTOT,JPRB)/tdiff
      else
         zmflops = 0.0_JPRB
         zthrput = 0.0_JPRB
      endif
      write(0,1002) numomp,ngptot,int(sum(zinfo(4,:))),ngpblks,nproma,-1,&
           & int(tdiff*1000.0_JPRB),int(zmflops),int(zthrput)

      CALL errhead(0)

      DEALLOCATE(zPLCRIT_AER)
      DEALLOCATE(zPICRIT_AER)
      DEALLOCATE(zPRE_ICE)
      DEALLOCATE(zPCCN)
      DEALLOCATE(zPNICE)
      DEALLOCATE(zPT)
      DEALLOCATE(zPQ)
      DEALLOCATE(zPVFA)
      DEALLOCATE(zPVFL)
      DEALLOCATE(zPVFI)
      DEALLOCATE(zPDYNA)
      DEALLOCATE(zPDYNL)
      DEALLOCATE(zPDYNI)
      DEALLOCATE(zPHRSW)
      DEALLOCATE(zPHRLW)
      DEALLOCATE(zPVERVEL)
      DEALLOCATE(zPAP)
      DEALLOCATE(zPAPH)
      DEALLOCATE(zPLSM)
      DEALLOCATE(LLCUM)
      DEALLOCATE(ITYPE)
      DEALLOCATE(zPLU)
      CALL errcalc(0,'PLUDE',zPLUDE)
      DEALLOCATE(zPLUDE)
      DEALLOCATE(zPSNDE)
      DEALLOCATE(zPMFU)
      DEALLOCATE(zPMFD)
      DEALLOCATE(zPA)
      DEALLOCATE(zPEXTRA)
      DEALLOCATE(zPCLV)
      DEALLOCATE(zPSUPSAT)
      CALL DEALLOCATE_STATE_ARRAY(zcml, ztendency_cml)
      CALL DEALLOCATE_STATE_ARRAY(ztmp, ztendency_tmp)

      ! Intent OUT
      CALL errcalc(0,'PCOVPTOT',zPCOVPTOT)
      DEALLOCATE(zPCOVPTOT)
      CALL errcalc(0,'PRAINFRAC_TOPRFZ',zPRAINFRAC_TOPRFZ)
      DEALLOCATE(zPRAINFRAC_TOPRFZ)
      CALL errcalc(0,'PFSQLF',zPFSQLF)
      DEALLOCATE(zPFSQLF)
      CALL errcalc(0,'PFSQIF',zPFSQIF)
      DEALLOCATE(zPFSQIF)
      CALL errcalc(0,'PFCQLNG',zPFCQLNG)
      DEALLOCATE(zPFCQLNG)
      CALL errcalc(0,'PFCQNNG',zPFCQNNG)
      DEALLOCATE(zPFCQNNG)
      CALL errcalc(0,'PFSQRF',zPFSQRF)
      DEALLOCATE(zPFSQRF)
      CALL errcalc(0,'PFSQSF',zPFSQSF)
      DEALLOCATE(zPFSQSF)
      CALL errcalc(0,'PFCQRNG',zPFCQRNG)
      DEALLOCATE(zPFCQRNG)
      CALL errcalc(0,'PFCQSNG',zPFCQSNG)
      DEALLOCATE(zPFCQSNG)
      CALL errcalc(0,'PFSQLTUR',zPFSQLTUR)
      DEALLOCATE(zPFSQLTUR)
      CALL errcalc(0,'PFSQITUR',zPFSQITUR)
      DEALLOCATE(zPFSQITUR)
      CALL errcalc(0,'PFPLSL',zPFPLSL)
      DEALLOCATE(zPFPLSL)
      CALL errcalc(0,'PFPLSN',zPFPLSN)
      DEALLOCATE(zPFPLSN)
      CALL errcalc(0,'PFHPSL',zPFHPSL)
      DEALLOCATE(zPFHPSL)
      CALL errcalc(0,'PFHPSN',zPFHPSN)
      DEALLOCATE(zPFHPSN)

!      CALL errcalc(0,'tendency_loc%u',ztendency_loc)
!      CALL errcalc(0,'tendency_loc%v',ztendency_loc)
!      CALL errcalc(0,'tendency_loc%o3',ztendency_loc)
      CALL errcalc(0,'tendency_loc%a',ztendency_loc)
      CALL errcalc(0,'tendency_loc%q',ztendency_loc)
      CALL errcalc(0,'tendency_loc%T',ztendency_loc)
      CALL errcalc(0,'tendency_loc%cld',ztendency_loc)

      CALL DEALLOCATE_STATE_ARRAY(zloc, ztendency_loc)

1004  format(5x,'=> END CASE#',i0)
      write(0,1004) JPR

      call loadjj(0) ! reset
   enddo
endif ! if (allocate(NPROMAS_IN)) then

CONTAINS

SUBROUTINE ALLOCATE_STATE_ARRAY(zbase,zstate,KPROMA,KLEV,KCLV,KGPBLKS)
TYPE(SPACE_T), INTENT(OUT) :: zbase
TYPE(STATE_TYPE), allocatable, INTENT(OUT) :: zstate(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA,KLEV,KCLV,KGPBLKS
INTEGER(KIND=JPIM) :: ICOUNT, jbl
allocate(zstate(KGPBLKS))
ICOUNT = 0
ALLOCATE(zbase%arrl(KPROMA,KLEV,6*KGPBLKS))
do jbl=1,KGPBLKS
   ICOUNT = ICOUNT + 1
   zstate(jbl)%u  => zbase%arrl(:,:,ICOUNT)
   ICOUNT = ICOUNT + 1
   zstate(jbl)%v  => zbase%arrl(:,:,ICOUNT)
   ICOUNT = ICOUNT + 1
   zstate(jbl)%o3 => zbase%arrl(:,:,ICOUNT)
   ICOUNT = ICOUNT + 1
   zstate(jbl)%a  => zbase%arrl(:,:,ICOUNT)
   ICOUNT = ICOUNT + 1
   zstate(jbl)%q  => zbase%arrl(:,:,ICOUNT)
   ICOUNT = ICOUNT + 1
   zstate(jbl)%T  => zbase%arrl(:,:,ICOUNT)
enddo
ICOUNT = 0
ALLOCATE(zbase%cld(KPROMA,KLEV,KCLV,1*KGPBLKS))
do jbl=1,KGPBLKS
   ICOUNT = ICOUNT + 1
   zstate(jbl)%cld => zbase%cld(:,:,:,ICOUNT)
enddo
END SUBROUTINE ALLOCATE_STATE_ARRAY

SUBROUTINE DEALLOCATE_STATE_ARRAY(zbase,zstate)
TYPE(SPACE_T), INTENT(INOUT) :: zbase
TYPE(STATE_TYPE), allocatable, INTENT(INOUT) :: zstate(:)
deallocate(zbase%arrl)
deallocate(zbase%cld)
deallocate(zstate)
END SUBROUTINE DEALLOCATE_STATE_ARRAY

END SUBROUTINE CLOUDSC_DRIVER
