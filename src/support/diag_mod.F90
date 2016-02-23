module diag_mod
USE PARKIND1  ,ONLY : JPIM
implicit none
save
public

INTEGER(KIND=JPIM) :: NUMOMP  = 1 ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOT  = 0 ! Number of grid points (as read from command line)

INTEGER(KIND=JPIM) :: NPROMA  = 0 ! NPROMA blocking factor (currently active)
INTEGER(KIND=JPIM) :: NGPBLKS = 0 ! Number of NPROMA-blocks (currently active)

INTEGER(KIND=JPIM), ALLOCATABLE :: NPROMAS_IN(:) ! NPROMAs as read from command line

end module diag_mod
