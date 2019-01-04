PROGRAM DWARF_CLOUDSC

USE PARKIND1, ONLY : JPIM
USE CLOUDSC_DRIVER, ONLY: CLOUDSC_GLOBAL_STATE, CLOUDSC_GLOBAL_STATE_LOAD, CLOUDSC_GLOBAL_STATE_VALIDATE

IMPLICIT NONE

CHARACTER(LEN=20) :: CLARG
INTEGER(KIND=JPIM) :: IARGS, LENARG, JARG, I

INTEGER(KIND=JPIM) :: NUMOMP  = 1     ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOT  = 16384 ! Number of grid points (as read from command line)
INTEGER(KIND=JPIM) :: NPROMA  = 32    ! NPROMA blocking factor (currently active)

TYPE(CLOUDSC_GLOBAL_STATE) :: GLOBAL_STATE

IARGS = COMMAND_ARGUMENT_COUNT()

! Get the number of OpenMP threads to use for the benchmark
if (IARGS >= 1) then
   CALL GET_COMMAND_ARGUMENT(1, CLARG, LENARG)
   READ(CLARG(1:LENARG),*) NUMOMP
end if

! Get total number of grid points (NGPTOT) with which to run the benchmark
IF (IARGS >= 2) THEN
  CALL GET_COMMAND_ARGUMENT(2, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NGPTOT
END IF

! Get the block size (NPROMA) for which to run the benchmark  
IF (IARGS >= 3) THEN
  CALL GET_COMMAND_ARGUMENT(3, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NPROMA
ENDIF

! TODO: Create a global global memory state from serialized input data
CALL CLOUDSC_GLOBAL_STATE_LOAD()


! TODO: Call the parallel execution loop with the CLOUDSC kernel
! For now, this should be a plain outer-parallel OpenMP loop,
! but later it will use a generic "parallel_loop" abstraction.


! TODO: Validate the output against serialized reference data
CALL CLOUDSC_GLOBAL_STATE_VALIDATE()

END PROGRAM DWARF_CLOUDSC
