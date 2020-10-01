MODULE VALIDATE_MOD
  USE PARKIND1, ONLY: JPIM, JPRB

  IMPLICIT NONE

  INTERFACE VALIDATE
     PROCEDURE VALIDATE_R1, VALIDATE_R2, VALIDATE_R3
  END INTERFACE VALIDATE

CONTAINS

  SUBROUTINE VALIDATE_R1(NAME, REF, FIELD, NLON, NGPTOT, NBLOCKS)
    ! Computes and prints errors "in the L1 norm sense"
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(INOUT) :: REF(:,:), FIELD(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON, NBLOCKS, NGPTOT

    INTEGER :: B, BSIZE, JK
    REAL(KIND=JPRB) :: ZMINVAL, ZMAXVAL, ZDIFF, ZMAXERR, ZERRSUM, ZSUM, ZRELERR, ZAVGPGP

    ZMINVAL = +HUGE(ZMINVAL)
    ZMAXVAL = -HUGE(ZMAXVAL)
    ZMAXERR = 0.0_JPRB
    ZERRSUM = 0.0_JPRB
    ZSUM = 0.0_JPRB

    !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
    !& REDUCTION(MIN:ZMINVAL, MAX:ZMAXVAL, MAX:ZMAXERR, +:ZERRSUM, +:ZSUM)
    DO B=1, NBLOCKS
      BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
      ZMINVAL = MIN(ZMINVAL,MINVAL(FIELD(:,B)))
      ZMAXVAL = MAX(ZMAXVAL,MAXVAL(FIELD(:,B)))
      DO JK=1, bsize
        ! Difference against reference result in one-norm sense
        ZDIFF = ABS(FIELD(JK,B) - REF(JK,B))
        ZMAXERR = MAX(ZMAXERR,ZDIFF)
        ZERRSUM = ZERRSUM + ZDIFF
        ZSUM = ZSUM + ABS(REF(JK,B))
      END DO
    END DO
    ZAVGPGP = ZERRSUM / REAL(NGPTOT,JPRB)

    CALL ERROR_PRINT(NAME, ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP, NDIM=1)
  END SUBROUTINE VALIDATE_R1

  SUBROUTINE VALIDATE_R2(NAME, REF, FIELD, NLON, NLEV, NGPTOT, NBLOCKS)
    ! Computes and prints errors "in the L1 norm sense"
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(INOUT) :: REF(:,:,:), FIELD(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON, NLEV, NBLOCKS, NGPTOT

    INTEGER :: B, BSIZE, JL, JK
    REAL(KIND=JPRB) :: ZMINVAL, ZMAXVAL, ZDIFF, ZMAXERR, ZERRSUM, ZSUM, ZRELERR, ZAVGPGP

    ZMINVAL = +HUGE(ZMINVAL)
    ZMAXVAL = -HUGE(ZMAXVAL)
    ZMAXERR = 0.0_JPRB
    ZERRSUM = 0.0_JPRB
    ZSUM = 0.0_JPRB

    !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
    !& REDUCTION(MIN:ZMINVAL, MAX:ZMAXVAL, MAX:ZMAXERR, +:ZERRSUM, +:ZSUM)
    DO B=1, NBLOCKS
      BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
      ZMINVAL = MIN(ZMINVAL,MINVAL(FIELD(:,:,B)))
      ZMAXVAL = MAX(ZMAXVAL,MAXVAL(FIELD(:,:,B)))
      DO JL=1, NLEV
        DO JK=1, bsize
          ! Difference against reference result in one-norm sense
          ZDIFF = ABS(FIELD(JK,JL,B) - REF(JK,JL,B))
          ZMAXERR = MAX(ZMAXERR,ZDIFF)
          ZERRSUM = ZERRSUM + ZDIFF
          ZSUM = ZSUM + ABS(REF(JK,JL,B))
        ENDDO
      END DO
    END DO
    ZAVGPGP = ZERRSUM/REAL(NGPTOT,JPRB)

    CALL ERROR_PRINT(NAME, ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP, NDIM=2)
  END SUBROUTINE VALIDATE_R2

  SUBROUTINE VALIDATE_R3(NAME, REF, FIELD, NLON, NLEV, NDIM, NGPTOT, NBLOCKS)
    ! Computes and prints errors "in the L1 norm sense"
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(INOUT) :: REF(:,:,:,:), FIELD(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON, NLEV, NDIM, NBLOCKS, NGPTOT

    INTEGER :: B, BSIZE, JL, JK, JM
    REAL(KIND=JPRB) :: ZMINVAL, ZMAXVAL, ZDIFF, ZMAXERR, ZERRSUM, ZSUM, ZRELERR, ZAVGPGP

    ZMINVAL = +HUGE(ZMINVAL)
    ZMAXVAL = -HUGE(ZMAXVAL)
    ZMAXERR = 0.0_JPRB
    ZERRSUM = 0.0_JPRB
    ZSUM = 0.0_JPRB

    !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
    !& REDUCTION(MIN:ZMINVAL, MAX:ZMAXVAL, MAX:ZMAXERR, +:ZERRSUM, +:ZSUM)
    DO B=1, NBLOCKS
      BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
      ZMINVAL = MIN(ZMINVAL,MINVAL(FIELD(:,:,:,B)))
      ZMAXVAL = MAX(ZMAXVAL,MAXVAL(FIELD(:,:,:,B)))
      DO JM=1, NDIM
        DO JL=1, NLEV
          DO JK=1, bsize
            ! Difference against reference result in one-norm sense
            ZDIFF = ABS(FIELD(JK,JL,JM,B) - REF(JK,JL,JM,B))
            ZMAXERR = MAX(ZMAXERR,ZDIFF)
            ZERRSUM = ZERRSUM + ZDIFF
            ZSUM = ZSUM + ABS(REF(JK,JL,JM,B))
          END DO
        END DO
      END DO
    END DO
    ZAVGPGP = ZERRSUM/REAL(NGPTOT,JPRB)

    CALL ERROR_PRINT(NAME, ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP, NDIM=3)
  END SUBROUTINE VALIDATE_R3

  SUBROUTINE ERROR_PRINT(NAME, ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP, NDIM)
    ! Print error statistic for a single variable (adapted from diff_mod.F90)
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(IN) :: ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP
    INTEGER(KIND=JPIM), INTENT(IN) :: NDIM
    REAL(KIND=JPRB) :: zrelerr
    REAL(KIND=JPRB), parameter :: zeps = epsilon(1.0_JPRB)
    INTEGER :: IOPT
    character(len=5) clwarn

    iopt = 0
    if (zerrsum < zeps) then
      zrelerr = 0.0_JPRB
      iopt = 1
    elseif (zsum < zeps) then
      zrelerr = zerrsum/(1.0_JPRB + zsum)
      iopt = 2
    else
      zrelerr = zerrsum/zsum
      iopt = 3
    endif

    !-- If you get 4 exclamation marks next to your error output,
    !   then it is likely that some uninitialized variables exists or
    !   some other screw-up -- watch out this !!!!
    clwarn = ' '
    if (zrelerr > 10.0_JPRB * zeps) clwarn = ' !!!!'
    zrelerr = 100.0_JPRB * zrelerr

    write(*,1000) name,ndim,iopt, &
     & zminval,zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn
1000 format(1X,A20,1X,I1,'D',I1,5(1X,E20.13),A)

  END SUBROUTINE ERROR_PRINT

END MODULE VALIDATE_MOD
