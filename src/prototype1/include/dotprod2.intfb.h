! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTERFACE
SUBROUTINE DOTPROD2(PGP1,PGP2,PPROD)
USE PARKIND1 ,ONLY : JPIM, JPRB
REAL(KIND=JPRB),INTENT(IN) :: PGP1(:,:,:)
REAL(KIND=JPRB),INTENT(IN) :: PGP2(:,:,:)
REAL(KIND=JPRB),INTENT(OUT) :: PPROD
END SUBROUTINE DOTPROD2
END INTERFACE
