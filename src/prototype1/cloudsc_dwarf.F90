! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM CLOUDSC_DWARF

USE PARKIND1  ,ONLY : JPIM
use diag_mod  ,ONLY : NUMOMP, NGPTOT, NPROMAS_IN

implicit none

#include "cloudsc_driver.intfb.h"

CHARACTER(LEN=20) :: clarg
INTEGER(KIND=JPIM) :: iargs, lenarg, jarg

iargs = COMMAND_ARGUMENT_COUNT()
NUMOMP = 1
NGPTOT = 0

if (iargs >= 1) then
   CALL GET_COMMAND_ARGUMENT(1, clarg, lenarg)
   read(clarg(1:lenarg),*) NUMOMP
   if (iargs >= 2) then
      CALL GET_COMMAND_ARGUMENT(2, clarg, lenarg)
      read(clarg(1:lenarg),*) NGPTOT
      if (iargs >= 3) then
         allocate(NPROMAS_IN(iargs-2))
         do jarg=1,size(NPROMAS_IN)
            CALL GET_COMMAND_ARGUMENT(jarg+2, clarg, lenarg)
            read(clarg(1:lenarg),*) NPROMAS_IN(jarg)
         enddo
      endif
   endif
endif

CALL CLOUDSC_DRIVER

END PROGRAM CLOUDSC_DWARF
