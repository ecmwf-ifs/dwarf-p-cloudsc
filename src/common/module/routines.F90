! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE routines
CHARACTER(LEN=*), PARAMETER :: CLPREFIX = 'ERROR: YOU SHOULD HAVE NEVER BE CALLING ROUTINE '
#include "abor1.intfb.h"
ENTRY DOTPROD2
CALL ABOR1(CLPREFIX//'DOTPROD2')
RETURN
ENTRY DOTPROD3
CALL ABOR1(CLPREFIX//'DOTPROD3')
RETURN
ENTRY VDIV
CALL ABOR1(CLPREFIX//'VDIV')
RETURN
ENTRY VEXP
CALL ABOR1(CLPREFIX//'VEXP')
RETURN
ENTRY VREC
CALL ABOR1(CLPREFIX//'VREC')
RETURN
ENTRY VPOW
CALL ABOR1(CLPREFIX//'VPOW')
RETURN
END SUBROUTINE routines
