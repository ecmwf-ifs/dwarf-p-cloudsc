! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM DWARF_CLOUDSC

USE PARKIND1, ONLY: JPIM
USE CLOUDSC_MPI_MOD, ONLY: CLOUDSC_MPI_INIT, CLOUDSC_MPI_END, NUMPROC, IRANK
USE CLOUDSC_GLOBAL_ATLAS_STATE_MOD, ONLY: CLOUDSC_GLOBAL_ATLAS_STATE
USE CLOUDSC_DRIVER_SCC_MOD, ONLY: CLOUDSC_DRIVER

USE ATLAS_MODULE
USE, INTRINSIC :: ISO_C_BINDING
USE ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS_MODULE
!USE ATLAS_MULTIFIELD_MODULE

#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

CHARACTER(LEN=20) :: CLARG
INTEGER(KIND=JPIM) :: IARGS, LENARG, JARG, I

INTEGER(KIND=JPIM) :: NUMOMP   = 1     ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOTG  = 16384 ! Number of grid points (as read from command line)
INTEGER(KIND=JPIM) :: NPROMA   = 32    ! NPROMA blocking factor (currently active)
INTEGER(KIND=JPIM) :: NITER    = 1     ! Call CLOUDSC driver multiple times, but deactivates validation when > 1
INTEGER(KIND=JPIM) :: JITER

TYPE(ATLAS_FIELDSET) :: FSET
TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
TYPE(CLOUDSC_GLOBAL_ATLAS_STATE) :: GLOBAL_ATLAS_STATE
TYPE(ATLAS_TRACE) :: TRACE

#include "abor1.intfb.h"

IARGS = COMMAND_ARGUMENT_COUNT()

! Get the number of OpenMP threads to use for the benchmark
if (IARGS >= 1) then
  CALL GET_COMMAND_ARGUMENT(1, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NUMOMP
  if (NUMOMP <= 0) then
#ifdef _OPENMP
    NUMOMP = OMP_GET_MAX_THREADS()
#else
    ! if arg is 0 or negative, and OpenMP disabled; defaults to 1
    NUMOMP = 1
#endif
  end if
end if

! Initialize MPI environment
CALL CLOUDSC_MPI_INIT(NUMOMP)
CALL ATLAS_LIBRARY%INITIALISE()
TRACE = ATLAS_TRACE("dwarf_cloudsc_atlas.F90",__LINE__,"program")

! Get total number of grid points (NGPTOTG) with which to run the benchmark
IF (IARGS >= 2) THEN
  CALL GET_COMMAND_ARGUMENT(2, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NGPTOTG
END IF

! Get the block size (NPROMA) for which to run the benchmark  
IF (IARGS >= 3) THEN
  CALL GET_COMMAND_ARGUMENT(3, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NPROMA
ENDIF

! timing of memory pool allocations
FSET = ATLAS_FIELDSET()
CALL GLOBAL_ATLAS_STATE%LOAD(FSET, FSPACE, NPROMA, NGPTOTG)
write(0,*) " ### Ignore the above timer, it is timing the one-time memory pool allocation"

FSET = ATLAS_FIELDSET()
CALL GLOBAL_ATLAS_STATE%LOAD(FSET, FSPACE, NPROMA, NGPTOTG)
write(0,*) " ### Above timer is host allocation"

CALL GET_ENV_INT("NITER",NITER)

! Call the driver to perform the parallel loop over our kernel
DO  JITER = 1, NITER
  write(0,'(A,I0,A,I0)') "### ITERATION ", JITER, '/', NITER
  CALL CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOTG, GLOBAL_ATLAS_STATE%KFLDX, GLOBAL_ATLAS_STATE%PTSPHY)

  ! Validate the output against serialized reference data
  IF (JITER == NITER) CALL GLOBAL_ATLAS_STATE%VALIDATE(FSET, FSPACE, NGPTOTG)
ENDDO

CALL FSET%FINAL()
CALL FSPACE%FINAL()
CALL TRACE%FINAL()

! Tear down MPI environment
CALL ATLAS_LIBRARY%FINALISE()
CALL CLOUDSC_MPI_END()

CONTAINS

SUBROUTINE GET_ENV_INT(name,value)
    character(len=*), intent(in):: name
    INTEGER, intent(inout) :: value
    character(len=32) :: value_str
    integer :: value_str_length
    CALL GET_ENVIRONMENT_VARIABLE(NAME=trim(name), VALUE=value_str, LENGTH=value_str_length)
    IF (value_str_length > 0 ) THEN
      READ(value_str(1:value_str_length), *) value
    ENDIF
END SUBROUTINE

END PROGRAM DWARF_CLOUDSC
