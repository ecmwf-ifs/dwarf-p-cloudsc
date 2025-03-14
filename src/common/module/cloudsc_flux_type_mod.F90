! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_FLUX_TYPE_MOD
  USE PARKIND1,  ONLY : JPIM, JPRB

  USE FIELD_MODULE, ONLY: FIELD_3RB, FIELD_4RB, FIELD_3RB_PTR
  USE FIELD_FACTORY_MODULE, ONLY: FIELD_NEW, FIELD_DELETE

  IMPLICIT NONE

  TYPE CLOUDSC_FLUX_TYPE
  
    INTEGER(KIND=JPIM) :: NLEV
    LOGICAL :: PACKED
    
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQLF(:,:)   ! Flux of liquid
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQIF(:,:)   ! Flux of ice
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFCQLNG(:,:)  ! -ve corr for liq
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFCQNNG(:,:)  ! -ve corr for ice
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQRF(:,:)   ! Flux diagnostics
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQSF(:,:)   ! for DDH, generic
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFCQRNG(:,:)  ! rain
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFCQSNG(:,:)  ! snow
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQLTUR(:,:) ! liquid flux due to VDF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFSQITUR(:,:) ! ice flux due to VDF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFPLSL(:,:)   ! liq+rain sedim flux
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFPLSN(:,:)   ! ice+snow sedim flux
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFHPSL(:,:)   ! Enthalpy flux for liq
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PFHPSN(:,:)   ! Enthalp flux for ice
    
    CLASS(FIELD_3RB), POINTER :: F_PFSQLF, F_PFSQIF, F_PFCQLNG, F_PFCQNNG, F_PFSQRF, F_PFSQSF, &
      & F_PFCQRNG, F_PFCQSNG, F_PFSQLTUR, F_PFSQITUR, F_PFPLSL, F_PFPLSN, F_PFHPSL, F_PFHPSN
  
    CLASS(FIELD_4RB), POINTER :: DATA_WRONLY ! ACTUAL FIELD storing data
    TYPE(FIELD_3RB_PTR), PRIVATE, ALLOCATABLE :: FIELDS_WRONLY(:) ! Array of field pointers
  
  CONTAINS
    PROCEDURE :: INIT => FLUX_TYPE_INIT
    PROCEDURE :: UPDATE_VIEW => FLUX_TYPE_UPDATE_VIEW
    PROCEDURE :: SYNC_HOST => FLUX_TYPE_SYNC_HOST
    PROCEDURE :: FINAL => FLUX_TYPE_FINAL
  
  END TYPE CLOUDSC_FLUX_TYPE

CONTAINS
  
  SUBROUTINE FLUX_TYPE_INIT(SELF,NPROMA, NGPTOT, KLON, KLEV, KFLDX, NBLOCKS, NGPTOTG, USE_PACKED)
    CLASS(CLOUDSC_FLUX_TYPE) :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: NPROMA, NGPTOT, KLON, KLEV, KFLDX, NBLOCKS 
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG
    LOGICAL, INTENT(IN), OPTIONAL :: USE_PACKED
    
    INTEGER(KIND=JPIM), PARAMETER :: NFIELDS = 14

    SELF%PACKED = .FALSE.
    IF (PRESENT(USE_PACKED)) SELF%PACKED = USE_PACKED
  
    
    IF (SELF%PACKED) THEN
      CALL FIELD_NEW(SELF%DATA_WRONLY, SELF%FIELDS_WRONLY, UBOUNDS=[NPROMA, KLEV+1, NFIELDS, NBLOCKS], &
      &                          PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      SELF%F_PFSQLF   => SELF%FIELDS_WRONLY(1)%PTR
      SELF%F_PFSQIF   => SELF%FIELDS_WRONLY(2)%PTR
      SELF%F_PFCQLNG  => SELF%FIELDS_WRONLY(3)%PTR
      SELF%F_PFCQNNG  => SELF%FIELDS_WRONLY(4)%PTR
      SELF%F_PFSQRF   => SELF%FIELDS_WRONLY(5)%PTR
      SELF%F_PFSQSF   => SELF%FIELDS_WRONLY(6)%PTR
      SELF%F_PFCQRNG  => SELF%FIELDS_WRONLY(7)%PTR
      SELF%F_PFCQSNG  => SELF%FIELDS_WRONLY(8)%PTR
      SELF%F_PFSQLTUR => SELF%FIELDS_WRONLY(9)%PTR
      SELF%F_PFSQITUR => SELF%FIELDS_WRONLY(10)%PTR
      SELF%F_PFPLSL   => SELF%FIELDS_WRONLY(11)%PTR
      SELF%F_PFPLSN   => SELF%FIELDS_WRONLY(12)%PTR
      SELF%F_PFHPSL   => SELF%FIELDS_WRONLY(13)%PTR
      SELF%F_PFHPSN   => SELF%FIELDS_WRONLY(14)%PTR
    ELSE
      CALL FIELD_NEW(SELF%F_PFSQLF, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFSQIF, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFCQLNG, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFCQNNG, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFSQRF, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFSQSF, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFCQRNG, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFCQSNG, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFSQLTUR, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFSQITUR, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFPLSL, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFPLSN, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFHPSL, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
      CALL FIELD_NEW(SELF%F_PFHPSN, UBOUNDS=[NPROMA, KLEV+1, NBLOCKS], PERSISTENT=.TRUE., INIT_VALUE=0._JPRB)
    END IF

  END SUBROUTINE FLUX_TYPE_INIT

  SUBROUTINE FLUX_TYPE_UPDATE_VIEW(SELF, BLOCK_INDEX)
    CLASS(CLOUDSC_FLUX_TYPE) :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    
   IF(ASSOCIATED(SELF%F_PFSQLF))    SELF%PFSQLF   => SELF%F_PFSQLF%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFSQIF))    SELF%PFSQIF   => SELF%F_PFSQIF%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFCQLNG))   SELF%PFCQLNG  => SELF%F_PFCQLNG%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFCQNNG))   SELF%PFCQNNG  => SELF%F_PFCQNNG%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFSQRF))    SELF%PFSQRF   => SELF%F_PFSQRF%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFSQSF))    SELF%PFSQSF   => SELF%F_PFSQSF%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFCQRNG))   SELF%PFCQRNG  => SELF%F_PFCQRNG%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFCQSNG))   SELF%PFCQSNG  => SELF%F_PFCQSNG%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFSQLTUR))  SELF%PFSQLTUR => SELF%F_PFSQLTUR%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFSQITUR))  SELF%PFSQITUR => SELF%F_PFSQITUR%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFPLSL))    SELF%PFPLSL   => SELF%F_PFPLSL%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFPLSN))    SELF%PFPLSN   => SELF%F_PFPLSN%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFHPSL))    SELF%PFHPSL   => SELF%F_PFHPSL%GET_VIEW(BLOCK_INDEX)
   IF(ASSOCIATED(SELF%F_PFHPSN))    SELF%PFHPSN   => SELF%F_PFHPSN%GET_VIEW(BLOCK_INDEX)

  END SUBROUTINE FLUX_TYPE_UPDATE_VIEW

  SUBROUTINE FLUX_TYPE_SYNC_HOST(SELF)
    CLASS(CLOUDSC_FLUX_TYPE) :: SELF
    
    IF (SELF%PACKED) THEN
      CALL SELF%DATA_WRONLY%SYNC_HOST_RDWR()
    ELSE 
      CALL SELF%F_PFSQLF%SYNC_HOST_RDWR()
      CALL SELF%F_PFSQIF%SYNC_HOST_RDWR()
      CALL SELF%F_PFCQLNG%SYNC_HOST_RDWR()
      CALL SELF%F_PFCQNNG%SYNC_HOST_RDWR()
      CALL SELF%F_PFSQRF%SYNC_HOST_RDWR()
      CALL SELF%F_PFSQSF%SYNC_HOST_RDWR()
      CALL SELF%F_PFCQRNG%SYNC_HOST_RDWR()
      CALL SELF%F_PFCQSNG%SYNC_HOST_RDWR()
      CALL SELF%F_PFSQLTUR%SYNC_HOST_RDWR()
      CALL SELF%F_PFSQITUR%SYNC_HOST_RDWR()
      CALL SELF%F_PFPLSL%SYNC_HOST_RDWR()
      CALL SELF%F_PFPLSN%SYNC_HOST_RDWR()
      CALL SELF%F_PFHPSL%SYNC_HOST_RDWR()
      CALL SELF%F_PFHPSN%SYNC_HOST_RDWR()
    ENDIF
  END SUBROUTINE FLUX_TYPE_SYNC_HOST

  SUBROUTINE FLUX_TYPE_FINAL(SELF)
    CLASS(CLOUDSC_FLUX_TYPE) :: SELF
    
    IF (SELF%PACKED) THEN
      CALL FIELD_DELETE(SELF%DATA_WRONLY)
      DEALLOCATE(SELF%FIELDS_WRONLY)
    ELSE
      CALL FIELD_DELETE(SELF%F_PFSQLF)
      CALL FIELD_DELETE(SELF%F_PFSQIF)
      CALL FIELD_DELETE(SELF%F_PFCQLNG)
      CALL FIELD_DELETE(SELF%F_PFCQNNG)
      CALL FIELD_DELETE(SELF%F_PFSQRF)
      CALL FIELD_DELETE(SELF%F_PFSQSF)
      CALL FIELD_DELETE(SELF%F_PFCQRNG)
      CALL FIELD_DELETE(SELF%F_PFCQSNG)
      CALL FIELD_DELETE(SELF%F_PFSQLTUR)
      CALL FIELD_DELETE(SELF%F_PFSQITUR)
      CALL FIELD_DELETE(SELF%F_PFPLSL)
      CALL FIELD_DELETE(SELF%F_PFPLSN)
      CALL FIELD_DELETE(SELF%F_PFHPSL)
      CALL FIELD_DELETE(SELF%F_PFHPSN)
    ENDIF
  END SUBROUTINE FLUX_TYPE_FINAL

END MODULE CLOUDSC_FLUX_TYPE_MOD

