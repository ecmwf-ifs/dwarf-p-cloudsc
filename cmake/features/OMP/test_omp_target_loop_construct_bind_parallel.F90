PROGRAM OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL
USE iso_fortran_env
IMPLICIT NONE

INTEGER, PARAMETER :: NB = 10
INTEGER, PARAMETER :: N = 10
INTEGER :: I, J
REAL(KIND=REAL32) :: TMP(N, NB)

!$omp target data map(tofrom: TMP)

!$omp target teams loop
DO I=1,NB
!$omp loop bind(parallel)
    DO J=1,N
        TMP(J, I) = REAL(J * I, KIND=REAL32)
    ENDDO
ENDDO

!$omp end target data

END PROGRAM OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL
