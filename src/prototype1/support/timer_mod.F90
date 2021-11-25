! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module timer_mod
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPIB
implicit none
CONTAINS
function ftimer()
implicit none
REAL(KIND=JPRB) :: ftimer
INTEGER(KIND=JPIB) :: t, rate
call system_clock(t,count_rate=rate)
ftimer = real(t,kind(ftimer))/real(rate,kind(ftimer))
end function ftimer
end module timer_mod
