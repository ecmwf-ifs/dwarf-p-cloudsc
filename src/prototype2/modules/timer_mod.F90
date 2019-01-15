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
