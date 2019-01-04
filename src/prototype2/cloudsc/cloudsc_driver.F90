MODULE CLOUDSC_DRIVER
  ! Driver module to manage the setup and teardown of the dwarf memory state
  USE PARKIND1, ONLY: JPIM, JPRB

  IMPLICIT NONE

  TYPE CLOUDSC_GLOBAL_STATE
    ! Memory state containing raw field variables for CLOUDSC dwarf
    INTEGER(KIND=JPIM) :: NGPTOT, KLEV  ! Number of grid points and vertical levels
  END TYPE CLOUDSC_GLOBAL_STATE

CONTAINS

  SUBROUTINE CLOUDSC_EXECUTE_KERNEL()
    ! Execute the benchmark by calling the CLOUDSC kernel
    ! with an outer-parallel block loop (NPROMA).

  END SUBROUTINE CLOUDSC_EXECUTE_KERNEL


  SUBROUTINE CLOUDSC_GLOBAL_STATE_LOAD()
    ! Load reference input data via serialbox

  END SUBROUTINE CLOUDSC_GLOBAL_STATE_LOAD


  SUBROUTINE CLOUDSC_GLOBAL_STATE_VALIDATE()
    ! Validate the correctness of output against reference data

  END SUBROUTINE CLOUDSC_GLOBAL_STATE_VALIDATE

END MODULE CLOUDSC_DRIVER
