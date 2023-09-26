! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM DWARF_CLOUDSC

USE PARKIND1, ONLY: JPIM, JPIB
USE CLOUDSC_MPI_MOD, ONLY: CLOUDSC_MPI_INIT, CLOUDSC_MPI_END, NUMPROC, IRANK
USE CLOUDSC_GLOBAL_ATLAS_STATE_MOD, ONLY: CLOUDSC_GLOBAL_ATLAS_STATE
USE CLOUDSC_DRIVER_MOD, ONLY: CLOUDSC_DRIVER
USE EC_PMON_MOD, ONLY: EC_PMON

USE ATLAS_MODULE
USE, INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

CHARACTER(LEN=20) :: CLARG
INTEGER(KIND=JPIM) :: IARGS, LENARG, JARG, I

INTEGER(KIND=JPIM) :: NUMOMP   = 1     ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOTG  = 16384 ! Number of grid points (as read from command line)
INTEGER(KIND=JPIM) :: NPROMA   = 32    ! NPROMA blocking factor (currently active)

type(atlas_fieldset) :: fset
type(atlas_field) :: field

TYPE(CLOUDSC_GLOBAL_ATLAS_STATE) :: GLOBAL_ATLAS_STATE

INTEGER(KIND=JPIB) :: ENERGY, POWER
CHARACTER(LEN=1)   :: CLEC_PMON

CALL GET_ENVIRONMENT_VARIABLE('EC_PMON', CLEC_PMON)
IF (CLEC_PMON == '1') THEN
  CALL EC_PMON(ENERGY, POWER)
  print *, "EC_PMON:: Initial (idle) power: ", POWER
END IF

IARGS = COMMAND_ARGUMENT_COUNT()

! Get the number of OpenMP threads to use for the benchmark
if (IARGS >= 1) then
   CALL GET_COMMAND_ARGUMENT(1, CLARG, LENARG)
   READ(CLARG(1:LENARG),*) NUMOMP
end if

! Initialize MPI environment
CALL CLOUDSC_MPI_INIT(NUMOMP)
CALL ATLAS_LIBRARY%INITIALISE()

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

FSET = ATLAS_FIELDSET()

! TODO: Create a global memory state from serialized input data
CALL GLOBAL_ATLAS_STATE%LOAD(FSET, NPROMA, NGPTOTG)

! Call the driver to perform the parallel loop over our kernel
CALL CLOUDSC_DRIVER(FSET, NUMOMP, NGPTOTG, GLOBAL_ATLAS_STATE%KFLDX, GLOBAL_ATLAS_STATE%PTSPHY)

! Validate the output against serialized reference data
CALL GLOBAL_ATLAS_STATE%VALIDATE(FSET, NGPTOTG)

! Tear down MPI environment
CALL ATLAS_LIBRARY%FINALISE()
CALL CLOUDSC_MPI_END()

END PROGRAM DWARF_CLOUDSC
