! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

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
