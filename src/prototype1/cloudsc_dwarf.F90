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
