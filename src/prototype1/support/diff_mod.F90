! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE diff_mod
USE PARKIND1  ,ONLY : JPIM, JPRB
USE diag_mod, ONLY  : NPROMA, NGPTOT, NGPBLKS, NUMOMP
USE expand_mod, ONLY : JJ
USE YOMPHYDER ,ONLY : STATE_TYPE

implicit none

save
private

interface saveref
   module procedure saveref_1d, saveref_2d, saveref_3d
end interface

interface errcalc
   module procedure errcalc_1d, errcalc_2d, errcalc_3d, errcalc_st
end interface

public :: saveref, errcalc
public :: errhead

INTEGER(KIND=JPIM), parameter :: NDIFFS = 5
INTEGER(KIND=JPIM), parameter :: NDIFF_MINVAL = 1
INTEGER(KIND=JPIM), parameter :: NDIFF_MAXVAL = 2
INTEGER(KIND=JPIM), parameter :: NDIFF_MAXERR = 3
INTEGER(KIND=JPIM), parameter :: NDIFF_ERRSUM = 4
INTEGER(KIND=JPIM), parameter :: NDIFF_SUM = 5

REAL(KIND=JPRB), parameter :: zeps = epsilon(1.0_JPRB)

TYPE refdata_t
   character(len=20) :: label
   INTEGER(KIND=JPIM) :: idim
   REAL(KIND=JPRB), allocatable :: d(:,:,:) ! typically klon x {klev or klev+1} x {1 or nclv}
   TYPE (refdata_t), pointer :: next => null()
END type refdata_t

TYPE (refdata_t), pointer :: beg => null()
TYPE (refdata_t), pointer :: cur => null()

CONTAINS

SUBROUTINE nextcur
if (.not.associated(cur)) then
   allocate(cur)
   beg => cur
else
   allocate(cur%next)
   cur => cur%next
endif
END SUBROUTINE nextcur

SUBROUTINE find(cdlabel,this)
character(len=*), intent(in) :: cdlabel
TYPE (refdata_t), pointer :: this
TYPE (refdata_t), pointer :: x
this => null()
x => beg
do while (associated(x))
   if (x%label == cdlabel) then
      this => x
      return
   endif
   x => x%next
end do
END SUBROUTINE find

SUBROUTINE saveref_1d(cdlabel, P)
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB) , intent(in) :: P(:)
INTEGER(KIND=JPIM) :: KLON
KLON = size(P,dim=1)
CALL nextcur
cur%label = cdlabel
cur%idim = 1
allocate(cur%d(KLON,1,1))
cur%d(:,1,1) = P(:)
END SUBROUTINE saveref_1d

SUBROUTINE saveref_2d(cdlabel, P)
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB) , intent(in) :: P(:,:)
INTEGER(KIND=JPIM) :: KLON, KLEV
KLON = size(P,dim=1)
KLEV = size(P,dim=2)
CALL nextcur
cur%label = cdlabel
cur%idim = 2
allocate(cur%d(KLON,KLEV,1))
cur%d(:,:,1) = P(:,:)
END SUBROUTINE saveref_2d

SUBROUTINE saveref_3d(cdlabel, P)
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB) , intent(in) :: P(:,:,:)
INTEGER(KIND=JPIM) :: KLON, KLEV, KCLV
KLON = size(P,dim=1)
KLEV = size(P,dim=2)
KCLV = size(P,dim=3)
CALL nextcur
cur%label = cdlabel
cur%idim = 3
allocate(cur%d(KLON,KLEV,KCLV))
cur%d(:,:,:) = P(:,:,:)
END SUBROUTINE saveref_3d

subroutine errcalc_st(iu, cdlabel, P)
INTEGER(KIND=JPIM), intent(in) :: iu
character(len=*), intent(in) :: cdlabel
TYPE(STATE_TYPE), intent(in) :: P(:) ! state_type with nproma_blocks
REAL(KIND=JPRB), allocatable :: Z(:,:,:,:) ! kproma x klev x kclv x nproma_blocks
INTEGER(KIND=JPIM) :: idx, ILEV, ICLV, JB
if (size(P) /= NGPBLKS) return
if (size(P) <= 0) return
idx = scan(cdlabel,'%')
if (idx <= 0) return ! not a state type struct
select case (cdlabel(idx+1:))
   case ('u')
      ILEV = size(P(1)%u,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%u ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('v')
      ILEV = size(P(1)%v,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%v ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('o3')
      ILEV = size(P(1)%o3,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%o3 ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('a')
      ILEV = size(P(1)%a,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%a ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('q')
      ILEV = size(P(1)%q,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%q ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('T')
      ILEV = size(P(1)%T,dim=2)
      allocate(Z(NPROMA,ILEV,NGPBLKS,1))
      DO JB=1,NGPBLKS ; Z(:,:,JB,1) = P(JB)%T ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,1))
   case ('cld')
      ILEV = size(P(1)%cld,dim=2)
      ICLV = size(P(1)%cld,dim=3)
      allocate(Z(NPROMA,ILEV,ICLV,NGPBLKS))
      DO JB=1,NGPBLKS ; Z(:,:,:,JB) = P(JB)%cld ; ENDDO
      CALL errcalc(iu, cdlabel, Z(:,:,:,:))
end select
if (allocated(z)) deallocate(z)
end subroutine errcalc_st

subroutine errcalc_1d(iu, cdlabel, P)
INTEGER(KIND=JPIM), intent(in) :: iu
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB), intent(in) :: P(:,:) ! nproma x nproma_blocks
INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, JL, JK, JM, ILON
INTEGER(KIND=JPIM),parameter :: ILEV = 1
INTEGER(KIND=JPIM),parameter :: ICLV = 1
INTEGER(KIND=JPIM) :: JCACHE(NPROMA), J, idim, iopt
REAL(KIND=JPRB) :: z,zerrsum,zmaxerr,zsum,zdiff,zminval,zmaxval,zrelerr,zavgpgp
REAL(KIND=JPRB) :: diff(NDIFFS)
TYPE (refdata_t), pointer :: this => null()
if (iu < 0) return
CALL find(cdlabel,this)
if (.not.associated(this)) then
   write(0,*) 'ERRCALC_1D: Could not associate this with label='//cdlabel
   return
endif
diff(NDIFF_MINVAL) = +huge(z)
diff(NDIFF_MAXVAL) = -huge(z)
diff(NDIFF_MAXERR) = 0.0_JPRB
diff(NDIFF_ERRSUM) = 0.0_JPRB
diff(NDIFF_SUM) = 0.0_JPRB
ILON = size(this%d,dim=1)
JM = 1
!$omp parallel default(shared) private(JKGLO, IBL, ICEND, JL, JK) &
!$omp& private(J,JCACHE,zerrsum,zsum,zmaxerr,zdiff,zminval,zmaxval)
zminval = +huge(z)
zmaxval = -huge(z)
zmaxerr = 0.0_JPRB
zerrsum = 0.0_JPRB
zsum = 0.0_JPRB
!$omp do
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,ILON)+1
   ENDDO
   DO JK=1,ILEV
      DO JL=1,ICEND
         J = JCACHE(JL)
         ! Difference against previously stored original-run result on Host (in one-norm sense)
         zminval = min(zminval,P(JL,IBL))
         zmaxval = max(zmaxval,P(JL,IBL))
         zdiff = abs(P(JL,IBL) - this%d(J,JK,JM))
         zmaxerr = MAX(zmaxerr,zdiff)
         zerrsum = zerrsum + zdiff
         zsum = zsum + abs(this%d(J,JK,JM))
      ENDDO
   ENDDO
ENDDO
!$omp end do
!$omp critical (oned)
diff(NDIFF_MINVAL) = MIN(diff(NDIFF_MINVAL),zminval)
diff(NDIFF_MAXVAL) = MAX(diff(NDIFF_MAXVAL),zmaxval)
diff(NDIFF_MAXERR) = MAX(diff(NDIFF_MAXERR),zmaxerr)
diff(NDIFF_ERRSUM) = diff(NDIFF_ERRSUM) + zerrsum
diff(NDIFF_SUM) = diff(NDIFF_SUM) + zsum
!$omp end critical (oned)
!$omp end parallel
CALL errprt(iu,cdlabel,this%idim,diff)
end subroutine errcalc_1d

subroutine errcalc_2d(iu, cdlabel, P)
INTEGER(KIND=JPIM), intent(in) :: iu
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB), intent(in) :: P(:,:,:) ! nproma x levels x nproma_blocks
INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, JL, JK, JM, ILON, ILEV
INTEGER(KIND=JPIM),parameter :: ICLV = 1
INTEGER(KIND=JPIM) :: JCACHE(NPROMA), J, idim, iopt
REAL(KIND=JPRB) :: z,zerrsum,zmaxerr,zsum,zdiff,zminval,zmaxval,zrelerr,zavgpgp
REAL(KIND=JPRB) :: diff(NDIFFS)
TYPE (refdata_t), pointer :: this => null()
if (iu < 0) return
CALL find(cdlabel,this)
if (.not.associated(this)) then
   write(0,*) 'ERRCALC_2D: Could not associate this with label='//cdlabel
   return
endif
diff(NDIFF_MINVAL) = +huge(z)
diff(NDIFF_MAXVAL) = -huge(z)
diff(NDIFF_MAXERR) = 0.0_JPRB
diff(NDIFF_ERRSUM) = 0.0_JPRB
diff(NDIFF_SUM) = 0.0_JPRB
ILEV = size(P,dim=2)
ILON = size(this%d,dim=1)
JM = 1
!$omp parallel default(shared) private(JKGLO, IBL, ICEND, JL, JK) &
!$omp& private(J,JCACHE,zerrsum,zsum,zmaxerr,zdiff,zminval,zmaxval)
zminval = +huge(z)
zmaxval = -huge(z)
zmaxerr = 0.0_JPRB
zerrsum = 0.0_JPRB
zsum = 0.0_JPRB
!$omp do
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,ILON)+1
   ENDDO
   DO JK=1,ILEV
      DO JL=1,ICEND
         J = JCACHE(JL)
         ! Difference against previously stored original-run result on Host (in one-norm sense)
         zminval = min(zminval,P(JL,JK,IBL))
         zmaxval = max(zmaxval,P(JL,JK,IBL))
         zdiff = abs(P(JL,JK,IBL) - this%d(J,JK,JM))
         zmaxerr = MAX(zmaxerr,zdiff)
         zerrsum = zerrsum + zdiff
         zsum = zsum + abs(this%d(J,JK,JM))
      ENDDO
   ENDDO
ENDDO
!$omp end do
!$omp critical (twod)
diff(NDIFF_MINVAL) = MIN(diff(NDIFF_MINVAL),zminval)
diff(NDIFF_MAXVAL) = MAX(diff(NDIFF_MAXVAL),zmaxval)
diff(NDIFF_MAXERR) = MAX(diff(NDIFF_MAXERR),zmaxerr)
diff(NDIFF_ERRSUM) = diff(NDIFF_ERRSUM) + zerrsum
diff(NDIFF_SUM) = diff(NDIFF_SUM) + zsum
!$omp end critical (twod)
!$omp end parallel
CALL errprt(iu,cdlabel,this%idim,diff)
end subroutine errcalc_2d

subroutine errcalc_3d(iu, cdlabel, P)
INTEGER(KIND=JPIM), intent(in) :: iu
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB), intent(in) :: P(:,:,:,:) ! nproma x levels x nclv x nproma_blocks
INTEGER(KIND=JPIM) :: JKGLO, IBL, ICEND, JL, JK, JM, ILON, ILEV, ICLV
INTEGER(KIND=JPIM) :: JCACHE(NPROMA), J, idim, iopt
REAL(KIND=JPRB) :: z,zerrsum,zmaxerr,zsum,zdiff,zminval,zmaxval,zrelerr,zavgpgp
REAL(KIND=JPRB) :: diff(NDIFFS)
TYPE (refdata_t), pointer :: this => null()
if (iu < 0) return
CALL find(cdlabel,this)
if (.not.associated(this)) then
   write(0,*) 'ERRCALC_3D: Could not associate this with label='//cdlabel
   return
endif
diff(NDIFF_MINVAL) = +huge(z)
diff(NDIFF_MAXVAL) = -huge(z)
diff(NDIFF_MAXERR) = 0.0_JPRB
diff(NDIFF_ERRSUM) = 0.0_JPRB
diff(NDIFF_SUM) = 0.0_JPRB
ILEV = size(P,dim=2)
ICLV = size(P,dim=3)
ILON = size(this%d,dim=1)
!$omp parallel default(shared) private(JKGLO, IBL, ICEND, JL, JK) &
!$omp& private(J,JCACHE,zerrsum,zsum,zmaxerr,zdiff,zminval,zmaxval)
zminval = +huge(z)
zmaxval = -huge(z)
zmaxerr = 0.0_JPRB
zerrsum = 0.0_JPRB
zsum = 0.0_JPRB
!$omp do
DO JKGLO=1,NGPTOT,NPROMA
   IBL=(JKGLO-1)/NPROMA+1
   ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
   J = JJ(IBL)
   DO JL=1,ICEND
      JCACHE(JL) = J
      J=MOD(J,ILON)+1
   ENDDO
   DO JM=1,ICLV
      DO JK=1,ILEV
         DO JL=1,ICEND
            J = JCACHE(JL)
            ! Difference against previously stored original-run result on Host (in one-norm sense)
            zminval = min(zminval,P(JL,JK,JM,IBL))
            zmaxval = max(zmaxval,P(JL,JK,JM,IBL))
            zdiff = abs(P(JL,JK,JM,IBL) - this%d(J,JK,JM))
            zmaxerr = MAX(zmaxerr,zdiff)
            zerrsum = zerrsum + zdiff
            zsum = zsum + abs(this%d(J,JK,JM))
         ENDDO
      ENDDO
   ENDDO
ENDDO
!$omp end do
!$omp critical (threed)
diff(NDIFF_MINVAL) = MIN(diff(NDIFF_MINVAL),zminval)
diff(NDIFF_MAXVAL) = MAX(diff(NDIFF_MAXVAL),zmaxval)
diff(NDIFF_MAXERR) = MAX(diff(NDIFF_MAXERR),zmaxerr)
diff(NDIFF_ERRSUM) = diff(NDIFF_ERRSUM) + zerrsum
diff(NDIFF_SUM) = diff(NDIFF_SUM) + zsum
!$omp end critical (threed)
!$omp end parallel
CALL errprt(iu,cdlabel,this%idim,diff)
end subroutine errcalc_3d

SUBROUTINE errhead(iu)
INTEGER(KIND=JPIM), intent(in) :: iu
if (iu < 0) return
write(iu,1000) 'Variable','Dim',&
     & 'MinValue','MaxValue','AbsMaxErr','AvgAbsErr/GP','MaxRelErr-%'
1000 format(1X,A20,1X,A3,5(1X,A20))
END SUBROUTINE errhead

SUBROUTINE errprt(iu,cdlabel,idim,diff)
INTEGER(KIND=JPIM), intent(in) :: iu, idim
character(len=*), intent(in) :: cdlabel
REAL(KIND=JPRB), intent(in) :: diff(:)
REAL(KIND=JPRB) :: zminval,zmaxval,zmaxerr,zavgpgp,zsum,zrelerr
INTEGER(KIND=JPIM) iopt
character(len=5) clwarn
if (iu < 0) return

zminval = diff(NDIFF_MINVAL)
zmaxval = diff(NDIFF_MAXVAL)
zmaxerr = diff(NDIFF_MAXERR)
zavgpgp = diff(NDIFF_ERRSUM)/REAL(NGPTOT,JPRB)
zsum = diff(NDIFF_SUM)

iopt = 0
if (diff(NDIFF_ERRSUM) < zeps) then
   zrelerr = 0.0_JPRB
   iopt = 1
elseif (zsum < zeps) then
   zrelerr = diff(NDIFF_ERRSUM)/(1.0_JPRB + zsum)
   iopt = 2
else
   zrelerr = diff(NDIFF_ERRSUM)/zsum
   iopt = 3
endif
 
!-- If you get 4 exclamation marks next to your error output,
!   then it is likely that some uninitialized variables exists or
!   some other screw-up -- watch out this !!!!
 
clwarn = ' '
if (zrelerr > 10.0_JPRB * zeps) clwarn = ' !!!!'

zrelerr = 100.0_JPRB * zrelerr

write(iu,1000) cdlabel,idim,iopt, &
     & zminval,zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn
1000 format(1X,A20,1X,I1,'D',I1,5(1X,E20.13),A)

END SUBROUTINE errprt

END MODULE diff_mod
