! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module expand_mod

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPIB
USE diag_mod, only : NGPTOT, NPROMA, NGPBLKS, NUMOMP
USE YOMPHYDER ,ONLY : STATE_TYPE

implicit none

save
private

INTEGER(KIND=JPIM), allocatable :: JJ(:)
INTEGER(KIND=JPIM) :: ICALL = 0
INTEGER(KIND=JPIB) :: TOTBYTES = 0

interface expand
   module procedure &
        & expand_L1, expand_I1, &
        & expand_R1, expand_R2, expand_r3, &
        & expand_state
end interface

public :: expand
public :: loadjj
public :: jj

contains

subroutine loadjj(KLON)
INTEGER(KIND=JPIM), intent(in) :: KLON
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL
INTEGER(KIND=JPIB) :: toGB, toMB, ibytes
toGB(ibytes) = (ibytes + 1073741823_JPIB) / 1073741824_JPIB
toMB(ibytes) = (ibytes +    1048575_JPIB) / 1048576_JPIB
if (allocated(JJ)) deallocate(JJ)
NGPBLKS = 0
IF (KLON > 0) THEN
   ICALL = 0
   TOTBYTES = 0
   NGPBLKS = 0
   DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      NGPBLKS = MAX(NGPBLKS,IBL)
   ENDDO
ENDIF
!write(0,'(1x,a,i3,i15,a)') 'loadjj   ',ICALL,toMB(TOTBYTES),' (in MB)'
IF (KLON > 0) THEN
!   write(0,'(1x,a,i0)')       'loadjj : NGPBLKS=',NGPBLKS
   allocate(JJ(NGPBLKS))
   JJ(:) = -1
   J=1
   DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      JJ(IBL) = J
      DO JL=1,ICEND
         J=MOD(J,KLON)+1
      ENDDO
   ENDDO
ENDIF
end subroutine loadjj

subroutine expand_L1(KLON, LDIN, LDOUT)
INTEGER(KIND=JPIM), intent(in) :: KLON
LOGICAL, intent(in) :: LDIN(KLON)
LOGICAL, intent(out):: LDOUT(NPROMA,NGPBLKS)
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL
INTEGER(KIND=JPIB) :: IBYTES
IBYTES = KIND(LDOUT) * SIZE(LDOUT)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_L1',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      LDOUT(JL,IBL) = LDIN(J)
      J=MOD(J,KLON)+1
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JL=ICEND+1,NPROMA
         LDOUT(JL,IBL) = .FALSE.
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_L1

subroutine expand_I1(KLON, KIN, KOUT)
INTEGER(KIND=JPIM), intent(in) :: KLON
INTEGER(KIND=JPIM), intent(in) :: KIN(KLON)
INTEGER(KIND=JPIM), intent(out):: KOUT(NPROMA,NGPBLKS)
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL
INTEGER(KIND=JPIB) :: IBYTES
IBYTES = KIND(KOUT) * SIZE(KOUT)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_I1',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      KOUT(JL,IBL) = KIN(J)
      J=MOD(J,KLON)+1
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JL=ICEND+1,NPROMA
         KOUT(JL,IBL) = 0
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_I1

subroutine expand_R1(KLON, PIN, POUT)
INTEGER(KIND=JPIM), intent(in) :: KLON
REAL(KIND=JPRB), intent(in)    :: PIN(KLON)
REAL(KIND=JPRB), intent(out)   :: POUT(NPROMA,NGPBLKS)
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL
INTEGER(KIND=JPIB) :: IBYTES
IBYTES = KIND(POUT) * SIZE(POUT)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_R1',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      POUT(JL,IBL) = PIN(J)
      J=MOD(J,KLON)+1
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JL=ICEND+1,NPROMA
         POUT(JL,IBL) = 0.0_JPRB
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_R1

subroutine expand_R2(KLON, KLEV, PIN, POUT)
INTEGER(KIND=JPIM), intent(in) :: KLON, KLEV
REAL(KIND=JPRB), intent(in)    :: PIN(KLON,KLEV)
REAL(KIND=JPRB), intent(out)   :: POUT(NPROMA,KLEV,NGPBLKS)
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL, JK
INTEGER(KIND=JPIM) :: JCACHE(NPROMA)
INTEGER(KIND=JPIB) :: IBYTES
IBYTES = KIND(POUT) * SIZE(POUT)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_R2',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL, JK, JCACHE) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,KLON)+1
   ENDDO
   DO JK=1,KLEV
      DO JL=1,ICEND
         J = JCACHE(JL)
         POUT(JL,JK,IBL) = PIN(J,JK)
      ENDDO
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JK=1,KLEV
         DO JL=ICEND+1,NPROMA
            POUT(JL,JK,IBL) = 0.0_JPRB
         ENDDO
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_R2

subroutine expand_R3(KLON, KLEV, K3RD, PIN, POUT)
INTEGER(KIND=JPIM), intent(in) :: KLON, KLEV, K3RD
REAL(KIND=JPRB), intent(in)    :: PIN(KLON,KLEV,K3RD)
REAL(KIND=JPRB), intent(out)   :: POUT(NPROMA,KLEV,K3RD,NGPBLKS)
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL, JK, JM
INTEGER(KIND=JPIM) :: JCACHE(NPROMA)
INTEGER(KIND=JPIB) :: IBYTES
if (K3RD <= 0) return
IBYTES = KIND(POUT) * SIZE(POUT)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_R3',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL, JK, JM, JCACHE) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,KLON)+1
   ENDDO
   DO JM=1,K3RD
      DO JK=1,KLEV
         DO JL=1,ICEND
            J = JCACHE(JL)
            POUT(JL,JK,JM,IBL) = PIN(J,JK,JM)
         ENDDO
      ENDDO
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JM=1,K3RD
         DO JK=1,KLEV
            DO JL=ICEND+1,NPROMA
               POUT(JL,JK,JM,IBL) = 0.0_JPRB
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_R3

subroutine expand_state(KLON, KLEV, K3RD, TIN, TOUT)
INTEGER(KIND=JPIM), intent(in) :: KLON, KLEV, K3RD
TYPE(STATE_TYPE), intent(in)    :: TIN
TYPE(STATE_TYPE), intent(out)   :: TOUT(:) ! length NGPBLKS with elems 6xNPROMAxKLEV + 1xNPROMAxKLEV*NCLV = (6 + K3RD) * NPROMA * KLEV
INTEGER(KIND=JPIM) :: J, JKGLO, IBL, ICEND, JL, JK, JM
INTEGER(KIND=JPIM) :: JCACHE(NPROMA)
INTEGER(KIND=JPIB) :: IBYTES
if (K3RD <= 0) return
IBYTES = JPRB * SIZE(TOUT) * ((6 + K3RD) * NPROMA * KLEV)
TOTBYTES = TOTBYTES + IBYTES
ICALL = ICALL + 1
!write(0,'(1x,a,i3,2i15)') 'expand_state',ICALL,IBYTES
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J, JKGLO, IBL, ICEND, JL, JK, JM, JCACHE) &
!$OMP& num_threads(NUMOMP)
!$OMP DO SCHEDULE(RUNTIME)
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,KLON)+1
   ENDDO
   DO JK=1,KLEV
      DO JL=1,ICEND
         J = JCACHE(JL)
         TOUT(IBL)%u(JL,JK)  = TIN%u(J,JK)
         TOUT(IBL)%v(JL,JK)  = TIN%v(J,JK)
         TOUT(IBL)%o3(JL,JK) = TIN%o3(J,JK)
         TOUT(IBL)%a(JL,JK)  = TIN%a(J,JK)
         TOUT(IBL)%q(JL,JK)  = TIN%q(J,JK)
         TOUT(IBL)%T(JL,JK)  = TIN%T(J,JK)
      ENDDO
   ENDDO
   DO JM=1,K3RD
      DO JK=1,KLEV
         DO JL=1,ICEND
            J = JCACHE(JL)
            TOUT(IBL)%cld(JL,JK,JM) = TIN%cld(J,JK,JM)
         ENDDO
      ENDDO
   ENDDO
   IF (ICEND < NPROMA) THEN
      DO JK=1,KLEV
         DO JL=ICEND+1,NPROMA
            TOUT(IBL)%u(JL,JK)  = 0.0_JPRB
            TOUT(IBL)%v(JL,JK)  = 0.0_JPRB
            TOUT(IBL)%o3(JL,JK) = 0.0_JPRB
            TOUT(IBL)%a(JL,JK)  = 0.0_JPRB
            TOUT(IBL)%q(JL,JK)  = 0.0_JPRB
            TOUT(IBL)%T(JL,JK)  = 0.0_JPRB
         ENDDO
      ENDDO
      DO JM=1,K3RD
         DO JK=1,KLEV
            DO JL=ICEND+1,NPROMA
               TOUT(IBL)%cld(JL,JK,JM) = 0.0_JPRB
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
!$OMP END DO
!$OMP END PARALLEL
end subroutine expand_state

end module expand_mod
