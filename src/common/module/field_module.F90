! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE FIELD_MODULE
  ! The FIELD types provided by this module provide data abstractions that
  ! decouple data storage in memory from the data views used in thread-parallel
  ! sections of the code. They are intended to thinly wrap ATLAS_FIELD
  ! objects and provide additional features that may later be
  ! incorporated into Atlas. They can also provide backward-compatibility
  ! for non-Atlas execution modes.

USE PARKIND1, ONLY: JPIM, JPRB
USE OML_MOD, ONLY: OML_MAX_THREADS, OML_MY_THREAD
USE IEEE_ARITHMETIC, ONLY: IEEE_SIGNALING_NAN

USE CUDAFOR

use openacc

use iso_c_binding

IMPLICIT NONE

TYPE FIELD_2D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  REAL(KIND=JPRB), POINTER :: VIEW(:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  REAL(KIND=JPRB), POINTER :: PTR(:,:) => NULL()
  ! REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DATA(:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:) => NULL()

  ! Number of blocks used in the data layout
  INTEGER :: NBLOCKS

  ! Flag indicating whether this field stores real data
  LOGICAL :: ACTIVE = .FALSE.
  ! Flag indicating the use a single block-buffer per thread
  LOGICAL :: THREAD_BUFFER = .FALSE.
  ! Flag indicating whether we own the allocated base array
  LOGICAL :: OWNED = .TRUE.
  ! Flag indicating whether latest data currently resides on device
  LOGICAL :: ON_DEVICE = .FALSE.

CONTAINS

  PROCEDURE :: CLONE => FIELD_2D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_2D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_2D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_2D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_2D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_2D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_2D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_2D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_2D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_2D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_2D_DELETE_DEVICE
END TYPE FIELD_2D

TYPE FIELD_3D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  REAL(KIND=JPRB), POINTER :: VIEW(:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  REAL(KIND=JPRB), POINTER :: PTR(:,:,:) => NULL()
  ! REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:,:)
  ! REAL(KIND=JPRB), ALLOCATABLE, PINNED :: DATA(:,:,:)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DATA(:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! ! A separate data pointer that can be used to create
  ! ! a contiguous chunk of host memory to cleanly map to
  ! ! device, should the %DATA pointer be discontiguous.
  ! REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:) => NULL()

  REAL(KIND=JPRB), DEVICE, ALLOCATABLE :: DEVDATA(:,:,:)

  ! Number of blocks used in the data layout
  INTEGER :: NBLOCKS

  ! Flag indicating whether this field stores real data
  LOGICAL :: ACTIVE = .FALSE.
  ! Flag indicating the use a single block-buffer per thread
  LOGICAL :: THREAD_BUFFER = .FALSE.
  ! Flag indicating whether we own the allocated base array
  LOGICAL :: OWNED = .TRUE.
  ! Flag indicating whether latest data currently resides on device
  LOGICAL :: ON_DEVICE = .FALSE.

CONTAINS

  PROCEDURE :: CLONE => FIELD_3D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_3D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_3D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_3D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_3D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_3D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_3D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_3D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_3D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_3D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_3D_DELETE_DEVICE
END TYPE FIELD_3D

TYPE FIELD_4D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  REAL(KIND=JPRB), POINTER :: VIEW(:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  REAL(KIND=JPRB), POINTER :: PTR(:,:,:,:) => NULL()
  ! REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:,:,:)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DATA(:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:) => NULL()

  ! Number of blocks used in the data layout
  INTEGER :: NBLOCKS

  ! Flag indicating whether this field stores real data
  LOGICAL :: ACTIVE = .FALSE.
  ! Flag indicating the use a single block-buffer per thread
  LOGICAL :: THREAD_BUFFER = .FALSE.
  ! Flag indicating whether we own the allocated base array
  LOGICAL :: OWNED = .TRUE.
  ! Flag indicating whether latest data currently resides on device
  LOGICAL :: ON_DEVICE = .FALSE.

CONTAINS

  PROCEDURE :: CLONE => FIELD_4D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_4D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_4D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_4D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_4D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_4D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_4D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_4D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_4D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_4D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_4D_DELETE_DEVICE
END TYPE FIELD_4D


TYPE FIELD_INT2D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  INTEGER(KIND=JPIM), POINTER :: VIEW(:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  INTEGER(KIND=JPIM), POINTER :: PTR(:,:) => NULL()
  ! INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:,:)
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DATA(:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:) => NULL()

  ! Number of blocks used in the data layout
  INTEGER :: NBLOCKS

  ! Flag indicating whether this field stores real data
  LOGICAL :: ACTIVE = .FALSE.
  ! Flag indicating the use a single block-buffer per thread
  LOGICAL :: THREAD_BUFFER = .FALSE.
  ! Flag indicating whether we own the allocated base array
  LOGICAL :: OWNED = .TRUE.
  ! Flag indicating whether latest data currently resides on device
  LOGICAL :: ON_DEVICE = .FALSE.

CONTAINS

  PROCEDURE :: CLONE => FIELD_INT2D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_INT2D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_INT2D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_INT2D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_INT2D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_INT2D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_INT2D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_INT2D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_INT2D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_INT2D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_INT2D_DELETE_DEVICE
END TYPE FIELD_INT2D


TYPE FIELD_LOG2D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  LOGICAL, POINTER :: VIEW(:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  LOGICAL, POINTER :: PTR(:,:) => NULL()
  ! LOGICAL, ALLOCATABLE :: DATA(:,:)
  LOGICAL, POINTER, CONTIGUOUS :: DATA(:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  LOGICAL, POINTER, CONTIGUOUS :: BASE_PTR(:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:) => NULL()

  ! Number of blocks used in the data layout
  INTEGER :: NBLOCKS

  ! Flag indicating whether this field stores real data
  LOGICAL :: ACTIVE = .FALSE.
  ! Flag indicating the use a single block-buffer per thread
  LOGICAL :: THREAD_BUFFER = .FALSE.
  ! Flag indicating whether we own the allocated base array
  LOGICAL :: OWNED = .TRUE.
  ! Flag indicating whether latest data currently resides on device
  LOGICAL :: ON_DEVICE = .FALSE.

CONTAINS

  PROCEDURE :: CLONE => FIELD_LOG2D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_LOG2D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_LOG2D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_LOG2D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_LOG2D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_LOG2D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_LOG2D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_LOG2D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_LOG2D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_LOG2D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_LOG2D_DELETE_DEVICE
END TYPE FIELD_LOG2D


TYPE FIELD_2D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_2D), POINTER :: PTR => NULL()
END TYPE FIELD_2D_PTR

TYPE FIELD_2D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:) => NULL()
END TYPE FIELD_2D_VIEW
TYPE FIELD_3D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_3D), POINTER :: PTR => NULL()
END TYPE FIELD_3D_PTR

TYPE FIELD_3D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:,:) => NULL()
END TYPE FIELD_3D_VIEW
TYPE FIELD_4D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_4D), POINTER :: PTR => NULL()
END TYPE FIELD_4D_PTR

TYPE FIELD_4D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:,:,:) => NULL()
END TYPE FIELD_4D_VIEW


INTERFACE FIELD_2D
  MODULE PROCEDURE :: FIELD_2D_WRAP
  MODULE PROCEDURE :: FIELD_2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_2D_EMPTY
  MODULE PROCEDURE :: FIELD_2D_ALLOCATE
END INTERFACE

INTERFACE FIELD_3D
  MODULE PROCEDURE :: FIELD_3D_WRAP
  MODULE PROCEDURE :: FIELD_3D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_3D_EMPTY
  MODULE PROCEDURE :: FIELD_3D_ALLOCATE
END INTERFACE

INTERFACE FIELD_4D
  MODULE PROCEDURE :: FIELD_4D_WRAP
  MODULE PROCEDURE :: FIELD_4D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_4D_EMPTY
  MODULE PROCEDURE :: FIELD_4D_ALLOCATE
END INTERFACE


INTERFACE FIELD_INT2D
  MODULE PROCEDURE :: FIELD_INT2D_WRAP
  MODULE PROCEDURE :: FIELD_INT2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_INT2D_EMPTY
  MODULE PROCEDURE :: FIELD_INT2D_ALLOCATE
END INTERFACE


INTERFACE FIELD_LOG2D
  MODULE PROCEDURE :: FIELD_LOG2D_WRAP
  MODULE PROCEDURE :: FIELD_LOG2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_LOG2D_EMPTY
  MODULE PROCEDURE :: FIELD_LOG2D_ALLOCATE
END INTERFACE


INTERFACE FILL_BUFFER
  MODULE PROCEDURE :: FILL_BUFFER_2D, FILL_BUFFER_3D, FILL_BUFFER_4D
  MODULE PROCEDURE :: FILL_BUFFER_INT2D, FILL_BUFFER_LOG2D
END INTERFACE FILL_BUFFER

INTERFACE FIELD_CREATE_DEVICE
  MODULE PROCEDURE :: FIELD_2D_CREATE_DEVICE
  MODULE PROCEDURE :: FIELD_3D_CREATE_DEVICE
  MODULE PROCEDURE :: FIELD_4D_CREATE_DEVICE
  MODULE PROCEDURE :: FIELD_INT2D_CREATE_DEVICE
  MODULE PROCEDURE :: FIELD_LOG2D_CREATE_DEVICE
END INTERFACE FIELD_CREATE_DEVICE

INTERFACE FIELD_UPDATE_DEVICE
  MODULE PROCEDURE :: FIELD_2D_UPDATE_DEVICE
  MODULE PROCEDURE :: FIELD_3D_UPDATE_DEVICE
  MODULE PROCEDURE :: FIELD_4D_UPDATE_DEVICE
  MODULE PROCEDURE :: FIELD_INT2D_UPDATE_DEVICE
  MODULE PROCEDURE :: FIELD_LOG2D_UPDATE_DEVICE
END INTERFACE FIELD_UPDATE_DEVICE

INTERFACE FIELD_UPDATE_HOST
  MODULE PROCEDURE :: FIELD_2D_UPDATE_HOST
  MODULE PROCEDURE :: FIELD_3D_UPDATE_HOST
  MODULE PROCEDURE :: FIELD_4D_UPDATE_HOST
  MODULE PROCEDURE :: FIELD_INT2D_UPDATE_HOST
  MODULE PROCEDURE :: FIELD_LOG2D_UPDATE_HOST
END INTERFACE FIELD_UPDATE_HOST

INTERFACE FIELD_DELETE_DEVICE
  MODULE PROCEDURE :: FIELD_2D_DELETE_DEVICE
  MODULE PROCEDURE :: FIELD_3D_DELETE_DEVICE
  MODULE PROCEDURE :: FIELD_4D_DELETE_DEVICE
  MODULE PROCEDURE :: FIELD_INT2D_DELETE_DEVICE
  MODULE PROCEDURE :: FIELD_LOG2D_DELETE_DEVICE
END INTERFACE FIELD_DELETE_DEVICE

INTERFACE GET_DEVICE_DATA
  MODULE PROCEDURE :: FIELD_2D_GET_DEVICE_DATA
  MODULE PROCEDURE :: FIELD_3D_GET_DEVICE_DATA
  MODULE PROCEDURE :: FIELD_4D_GET_DEVICE_DATA
  MODULE PROCEDURE :: FIELD_INT2D_GET_DEVICE_DATA
  MODULE PROCEDURE :: FIELD_LOG2D_GET_DEVICE_DATA
END INTERFACE GET_DEVICE_DATA

INTERFACE FIELD_ENSURE_DEVICE
  MODULE PROCEDURE :: FIELD_2D_ENSURE_DEVICE
  MODULE PROCEDURE :: FIELD_3D_ENSURE_DEVICE
  MODULE PROCEDURE :: FIELD_4D_ENSURE_DEVICE
  MODULE PROCEDURE :: FIELD_INT2D_ENSURE_DEVICE
  MODULE PROCEDURE :: FIELD_LOG2D_ENSURE_DEVICE
END INTERFACE FIELD_ENSURE_DEVICE

INTERFACE FIELD_ENSURE_HOST
  MODULE PROCEDURE :: FIELD_2D_ENSURE_HOST
  MODULE PROCEDURE :: FIELD_3D_ENSURE_HOST
  MODULE PROCEDURE :: FIELD_4D_ENSURE_HOST
  MODULE PROCEDURE :: FIELD_INT2D_ENSURE_HOST
  MODULE PROCEDURE :: FIELD_LOG2D_ENSURE_HOST
END INTERFACE FIELD_ENSURE_HOST

CONTAINS
  SUBROUTINE FILL_BUFFER_2D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: BUFFER(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: IDX

    IDX = INDEX+1
    BUFFER(IDX:) = BUFFER(INDEX)
  END SUBROUTINE FILL_BUFFER_2D

  SUBROUTINE FILL_BUFFER_3D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: BUFFER(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      BUFFER(IDX:,I) = BUFFER(INDEX,I)
    END DO
  END SUBROUTINE FILL_BUFFER_3D

  SUBROUTINE FILL_BUFFER_4D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: BUFFER(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        BUFFER(IDX:,I,J) = BUFFER(INDEX,I,J)
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_4D

  SUBROUTINE FILL_BUFFER_INT2D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: BUFFER(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: IDX

    IDX = INDEX+1
    BUFFER(IDX:) = BUFFER(INDEX)
  END SUBROUTINE FILL_BUFFER_INT2D

  SUBROUTINE FILL_BUFFER_LOG2D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    LOGICAL, POINTER, INTENT(INOUT) :: BUFFER(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: IDX

    IDX = INDEX+1
    BUFFER(IDX:) = BUFFER(INDEX)
  END SUBROUTINE FILL_BUFFER_LOG2D

  FUNCTION FIELD_2D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_2D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(1)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_2D_EMPTY

  FUNCTION FIELD_3D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_3D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(2)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1),SHAPE(2)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_3D_EMPTY

  FUNCTION FIELD_4D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_4D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(3)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1),SHAPE(2),SHAPE(3)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_4D_EMPTY

  FUNCTION FIELD_INT2D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_INT2D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(1)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT2D_EMPTY

  FUNCTION FIELD_LOG2D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_LOG2D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(1)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG2D_EMPTY

  FUNCTION FIELD_2D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_2D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_2D_WRAP

  FUNCTION FIELD_3D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_3D_WRAP

  FUNCTION FIELD_4D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_4D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_4D_WRAP

  FUNCTION FIELD_INT2D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT2D_WRAP

  FUNCTION FIELD_LOG2D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_LOG2D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG2D_WRAP

  FUNCTION FIELD_2D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_2D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_2D_WRAP_PACKED

  FUNCTION FIELD_3D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_3D_WRAP_PACKED

  FUNCTION FIELD_4D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_4D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_4D_WRAP_PACKED

  FUNCTION FIELD_INT2D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_INT2D_WRAP_PACKED

  FUNCTION FIELD_LOG2D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_LOG2D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 2)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_LOG2D_WRAP_PACKED

  FUNCTION FIELD_2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    INTEGER(KIND=JPIM) :: istat, arrsize
    type(c_ptr) :: hptr

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_2D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ! ALLOCATE(SELF%DATA(SHAPE(1),NBLK))

    arrsize = SHAPE(1) * NBLK * sizeof(1.0_JPRB)
    istat = cudaSetDeviceFlags(cudadevicemaphost)
    istat = cudaHostAlloc(hptr, arrsize, cudaHostAllocMapped)
    call c_f_pointer(hptr, self%data, [SHAPE(1), NBLK] )

    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_2D_ALLOCATE

  FUNCTION FIELD_3D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(2)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    INTEGER(KIND=JPIM) :: istat, arrsize
    type(c_ptr) :: hptr

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_3D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ! ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),NBLK))

    arrsize = SHAPE(1) * SHAPE(2) * NBLK * sizeof(1.0_JPRB)
    istat = cudaSetDeviceFlags(cudadevicemaphost)
    istat = cudaHostAlloc(hptr, arrsize, cudaHostAllocMapped)
    call c_f_pointer(hptr, self%data, [SHAPE(1), SHAPE(2), NBLK] )

    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_3D_ALLOCATE

  FUNCTION FIELD_4D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_4D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(3)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    INTEGER(KIND=JPIM) :: istat, arrsize
    type(c_ptr) :: hptr

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_4D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ! ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),NBLK))

    arrsize = SHAPE(1) * SHAPE(2) * SHAPE(3) * NBLK * sizeof(1.0_JPRB)
    istat = cudaSetDeviceFlags(cudadevicemaphost)
    istat = cudaHostAlloc(hptr, arrsize, cudaHostAllocMapped)
    call c_f_pointer(hptr, self%data, [SHAPE(1), SHAPE(2), SHAPE(3), NBLK] )

    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_4D_ALLOCATE

  FUNCTION FIELD_INT2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    INTEGER(KIND=JPIM) :: istat, arrsize
    type(c_ptr) :: hptr

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_INT2D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ! ALLOCATE(SELF%DATA(SHAPE(1),NBLK))

    arrsize = SHAPE(1) * NBLK * sizeof(1.0_JPRB)
    istat = cudaSetDeviceFlags(cudadevicemaphost)
    istat = cudaHostAlloc(hptr, arrsize, cudaHostAllocMapped)
    call c_f_pointer(hptr, self%data, [SHAPE(1), NBLK] )

    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT2D_ALLOCATE

  FUNCTION FIELD_LOG2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_LOG2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    INTEGER(KIND=JPIM) :: istat, arrsize
    type(c_ptr) :: hptr

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_LOG2D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ! ALLOCATE(SELF%DATA(SHAPE(1),NBLK))

    arrsize = SHAPE(1) * NBLK * sizeof(1.0_JPRB)
    istat = cudaSetDeviceFlags(cudadevicemaphost)
    istat = cudaHostAlloc(hptr, arrsize, cudaHostAllocMapped)
    call c_f_pointer(hptr, self%data, [SHAPE(1), NBLK] )

    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG2D_ALLOCATE

  FUNCTION FIELD_2D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_2D) :: SELF
    TYPE(FIELD_2D), POINTER :: NEWOBJ

    ALLOCATE(NEWOBJ)
    ! For owned storage data, re-allocate but do not copy data over
    IF (SELF%OWNED) THEN
      ALLOCATE(NEWOBJ%DATA, MOLD=SELF%DATA)
      NEWOBJ%PTR => NEWOBJ%DATA
    ELSE
      NEWOBJ%PTR => SELF%PTR
    END IF
    NEWOBJ%VIEW => NULL()
    NEWOBJ%NBLOCKS = SELF%NBLOCKS
    NEWOBJ%THREAD_BUFFER = SELF%THREAD_BUFFER
    NEWOBJ%OWNED = .FALSE.
  END FUNCTION FIELD_2D_CLONE

  FUNCTION FIELD_3D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_3D) :: SELF
    TYPE(FIELD_3D), POINTER :: NEWOBJ

    ALLOCATE(NEWOBJ)
    ! For owned storage data, re-allocate but do not copy data over
    IF (SELF%OWNED) THEN
      ALLOCATE(NEWOBJ%DATA, MOLD=SELF%DATA)
      NEWOBJ%PTR => NEWOBJ%DATA
    ELSE
      NEWOBJ%PTR => SELF%PTR
    END IF
    NEWOBJ%VIEW => NULL()
    NEWOBJ%NBLOCKS = SELF%NBLOCKS
    NEWOBJ%THREAD_BUFFER = SELF%THREAD_BUFFER
    NEWOBJ%OWNED = .FALSE.
  END FUNCTION FIELD_3D_CLONE

  FUNCTION FIELD_4D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_4D) :: SELF
    TYPE(FIELD_4D), POINTER :: NEWOBJ

    ALLOCATE(NEWOBJ)
    ! For owned storage data, re-allocate but do not copy data over
    IF (SELF%OWNED) THEN
      ALLOCATE(NEWOBJ%DATA, MOLD=SELF%DATA)
      NEWOBJ%PTR => NEWOBJ%DATA
    ELSE
      NEWOBJ%PTR => SELF%PTR
    END IF
    NEWOBJ%VIEW => NULL()
    NEWOBJ%NBLOCKS = SELF%NBLOCKS
    NEWOBJ%THREAD_BUFFER = SELF%THREAD_BUFFER
    NEWOBJ%OWNED = .FALSE.
  END FUNCTION FIELD_4D_CLONE

  FUNCTION FIELD_INT2D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_INT2D) :: SELF
    TYPE(FIELD_INT2D), POINTER :: NEWOBJ

    ALLOCATE(NEWOBJ)
    ! For owned storage data, re-allocate but do not copy data over
    IF (SELF%OWNED) THEN
      ALLOCATE(NEWOBJ%DATA, MOLD=SELF%DATA)
      NEWOBJ%PTR => NEWOBJ%DATA
    ELSE
      NEWOBJ%PTR => SELF%PTR
    END IF
    NEWOBJ%VIEW => NULL()
    NEWOBJ%NBLOCKS = SELF%NBLOCKS
    NEWOBJ%THREAD_BUFFER = SELF%THREAD_BUFFER
    NEWOBJ%OWNED = .FALSE.
  END FUNCTION FIELD_INT2D_CLONE

  FUNCTION FIELD_LOG2D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_LOG2D) :: SELF
    TYPE(FIELD_LOG2D), POINTER :: NEWOBJ

    ALLOCATE(NEWOBJ)
    ! For owned storage data, re-allocate but do not copy data over
    IF (SELF%OWNED) THEN
      ALLOCATE(NEWOBJ%DATA, MOLD=SELF%DATA)
      NEWOBJ%PTR => NEWOBJ%DATA
    ELSE
      NEWOBJ%PTR => SELF%PTR
    END IF
    NEWOBJ%VIEW => NULL()
    NEWOBJ%NBLOCKS = SELF%NBLOCKS
    NEWOBJ%THREAD_BUFFER = SELF%THREAD_BUFFER
    NEWOBJ%OWNED = .FALSE.
  END FUNCTION FIELD_LOG2D_CLONE


  SUBROUTINE FIELD_2D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_2D_UPDATE_VIEW

  SUBROUTINE FIELD_3D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_3D_UPDATE_VIEW

  SUBROUTINE FIELD_4D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_4D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,:,:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:,:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_4D_UPDATE_VIEW

  SUBROUTINE FIELD_INT2D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT2D_UPDATE_VIEW


  SUBROUTINE FIELD_LOG2D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_LOG2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG2D_UPDATE_VIEW

  SUBROUTINE FIELD_2D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_2D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_2D_EXTRACT_VIEW

  SUBROUTINE FIELD_3D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_3D_EXTRACT_VIEW

  SUBROUTINE FIELD_4D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_4D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_4D_EXTRACT_VIEW

  SUBROUTINE FIELD_INT2D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT2D_EXTRACT_VIEW

  SUBROUTINE FIELD_LOG2D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG2D), TARGET :: SELF
    LOGICAL, POINTER, INTENT(INOUT) :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG2D_EXTRACT_VIEW

  FUNCTION FIELD_2D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_2D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE)) THEN
      IF (BLOCK_INDEX == SELF%NBLOCKS) THEN
        ! Fill the the buffer by replicating the last entry
        CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
      END IF
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = 0.0_JPRB
    END IF
  END FUNCTION FIELD_2D_GET_VIEW

  FUNCTION FIELD_3D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER :: VIEW_PTR(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE)) THEN
      IF (BLOCK_INDEX == SELF%NBLOCKS) THEN
        ! Fill the the buffer by replicating the last entry
        CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
      END IF
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:) = 0.0_JPRB
    END IF
  END FUNCTION FIELD_3D_GET_VIEW

  FUNCTION FIELD_4D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_4D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER :: VIEW_PTR(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE)) THEN
      IF (BLOCK_INDEX == SELF%NBLOCKS) THEN
        ! Fill the the buffer by replicating the last entry
        CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
      END IF
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:,:) = 0.0_JPRB
    END IF
  END FUNCTION FIELD_4D_GET_VIEW

  FUNCTION FIELD_INT2D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE)) THEN
      IF (BLOCK_INDEX == SELF%NBLOCKS) THEN
        ! Fill the the buffer by replicating the last entry
        CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
      END IF
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = 0.0_JPIM
    END IF
  END FUNCTION FIELD_INT2D_GET_VIEW


  FUNCTION FIELD_LOG2D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG2D), TARGET :: SELF
    LOGICAL, POINTER :: VIEW_PTR(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE)) THEN
      IF (BLOCK_INDEX == SELF%NBLOCKS) THEN
        ! Fill the the buffer by replicating the last entry
        CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
      END IF
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:) = .FALSE.
    END IF
  END FUNCTION FIELD_LOG2D_GET_VIEW


  SUBROUTINE FIELD_2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_2D_CREATE_DEVICE

  SUBROUTINE FIELD_3D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: ARRSIZE

    ARRSIZE = SIZE(SELF%PTR) * SIZEOF(1.0_JPRB)
    ALLOCATE(SELF%DEVDATA, MOLD=SELF%PTR)
    CALL ACC_MAP_DATA(SELF%PTR, SELF%DEVDATA, ARRSIZE)

    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_3D_CREATE_DEVICE

  SUBROUTINE FIELD_4D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_4D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_4D_CREATE_DEVICE

  SUBROUTINE FIELD_INT2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT2D_CREATE_DEVICE

  SUBROUTINE FIELD_LOG2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG2D_CREATE_DEVICE

  FUNCTION FIELD_2D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_2D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_2D_GET_DEVICE_DATA

  FUNCTION FIELD_3D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:)

    type(c_ptr) :: hptr

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      hptr = acc_hostptr(self%devdata)
      call c_f_pointer(hptr, devptr, shape(self%devdata))
    END IF
  END FUNCTION FIELD_3D_GET_DEVICE_DATA

  FUNCTION FIELD_4D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_4D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_4D_GET_DEVICE_DATA

  FUNCTION FIELD_INT2D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_INT2D_GET_DEVICE_DATA

  FUNCTION FIELD_LOG2D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG2D), TARGET :: SELF
    LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_LOG2D_GET_DEVICE_DATA

  SUBROUTINE FIELD_2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      !$acc update device(SELF%DATA(:,:))
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_2D_UPDATE_DEVICE

  SUBROUTINE FIELD_3D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    INTEGER(KIND=JPIM) :: istat, arrsize, blksize
    type(c_ptr) :: hptr
    integer(kind=jpim) :: shape(3)

    logical :: pres

    arrsize = size(self%ptr) * sizeof(1.0_JPRB)
    blksize = arrsize / self%nblocks
    ALLOCATE(SELF%DEVDATA, mold=SELF%PTR)

    IF (SELF%OWNED) THEN
      call acc_map_data(self%data, self%devdata, arrsize)
      call acc_memcpy_to_device(self%devdata(:,:,:), self%data(:,:,:), arrsize)

    ELSE
      ! TODO: This is a dirty trick to fool the OpenACC runtime!
      ! We allocate the associated data array (full size), so that we can
      ! add it to the OpenACC host-device map (it's contiguous!)
      ! Then, we copy the data in a strided fashio from the discontiguous pointer.
      ALLOCATE(SELF%DATA, MOLD=SELF%PTR)
      call acc_map_data(self%data, self%devdata, arrsize)
      DO IBL=1, SELF%NBLOCKS
        call acc_memcpy_to_device(self%devdata(:,:,ibl), self%base_ptr(:,:,self%fidx,ibl), blksize)
      END DO
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_3D_UPDATE_DEVICE

  SUBROUTINE FIELD_4D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_4D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      !$acc update device(SELF%DATA(:,:,:,:))
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_4D_UPDATE_DEVICE

  SUBROUTINE FIELD_INT2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      !$acc update device(SELF%DATA(:,:))
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT2D_UPDATE_DEVICE

  SUBROUTINE FIELD_LOG2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_LOG2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      !$acc update device(SELF%DATA(:,:))
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG2D_UPDATE_DEVICE

  SUBROUTINE FIELD_2D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_2D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc update host(SELF%DATA(:,:))
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:) = SELF%DEVPTR(:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_2D_UPDATE_HOST

  SUBROUTINE FIELD_3D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_3D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    INTEGER(KIND=JPIM) :: istat, arrsize, blksize
    type(c_ptr) :: hptr

    arrsize = size(self%ptr) * sizeof(1.0_JPRB)
    blksize = arrsize / self%nblocks

    IF (SELF%OWNED) THEN
      call acc_memcpy_from_device(self%data(:,:,:), self%devdata(:,:,:), arrsize)
      call acc_unmap_data(self%data)

    ELSE
      DO IBL=1, SELF%NBLOCKS
        call acc_memcpy_from_device(self%ptr(:,:,ibl), self%devdata(:,:,ibl), blksize)
      END DO
      call acc_unmap_data(self%data)
      DEALLOCATE(SELF%DATA)
    END IF

    DEALLOCATE(SELF%DEVDATA)
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_3D_UPDATE_HOST

  SUBROUTINE FIELD_4D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_4D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:,:) = SELF%DEVPTR(:,:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_4D_UPDATE_HOST

  SUBROUTINE FIELD_INT2D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_INT2D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:) = SELF%DEVPTR(:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT2D_UPDATE_HOST

  SUBROUTINE FIELD_LOG2D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_LOG2D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:) = SELF%DEVPTR(:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG2D_UPDATE_HOST
 
  SUBROUTINE FIELD_2D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_2D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_2D_DELETE_DEVICE

  SUBROUTINE FIELD_3D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF

    IF (SELF%OWNED) THEN
      CALL ACC_UNMAP_DATA(SELF%DATA)
    ELSE
      CALL ACC_UNMAP_DATA(SELF%PTR)
    END IF
    DEALLOCATE(SELF%DEVDATA)
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_3D_DELETE_DEVICE

  SUBROUTINE FIELD_4D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_4D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_4D_DELETE_DEVICE

  SUBROUTINE FIELD_INT2D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT2D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT2D_DELETE_DEVICE

  SUBROUTINE FIELD_LOG2D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG2D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG2D_DELETE_DEVICE

  SUBROUTINE FIELD_2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_2D_ENSURE_HOST

  SUBROUTINE FIELD_3D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_3D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_3D_ENSURE_HOST

  SUBROUTINE FIELD_4D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_4D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_4D_ENSURE_HOST

  SUBROUTINE FIELD_INT2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_INT2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_INT2D_ENSURE_HOST

  SUBROUTINE FIELD_LOG2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_LOG2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_LOG2D_ENSURE_HOST

  SUBROUTINE FIELD_2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_2D_ENSURE_DEVICE

  SUBROUTINE FIELD_3D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_3D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_3D_ENSURE_DEVICE

  SUBROUTINE FIELD_4D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_4D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_4D_ENSURE_DEVICE

  SUBROUTINE FIELD_INT2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_INT2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_INT2D_ENSURE_DEVICE

  SUBROUTINE FIELD_LOG2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_LOG2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_LOG2D_ENSURE_DEVICE

  SUBROUTINE FIELD_2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_2D_FINAL

  SUBROUTINE FIELD_3D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_3D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_3D_FINAL

  SUBROUTINE FIELD_4D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_4D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_4D_FINAL

  SUBROUTINE FIELD_INT2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_INT2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_INT2D_FINAL

  SUBROUTINE FIELD_LOG2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_LOG2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_LOG2D_FINAL

END MODULE FIELD_MODULE
