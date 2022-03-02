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

IMPLICIT NONE

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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
  REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:)

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

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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
  REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:) => NULL()

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

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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
  REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:,:,:)

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

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_5D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  REAL(KIND=JPRB), POINTER :: VIEW(:,:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  REAL(KIND=JPRB), POINTER :: PTR(:,:,:,:,:) => NULL()
  REAL(KIND=JPRB), ALLOCATABLE :: DATA(:,:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_5D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_5D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_5D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_5D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_5D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_5D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_5D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_5D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_5D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_5D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_5D_DELETE_DEVICE
END TYPE FIELD_5D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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
  INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:,:)

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

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_INT3D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  INTEGER(KIND=JPIM), POINTER :: VIEW(:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  INTEGER(KIND=JPIM), POINTER :: PTR(:,:,:) => NULL()
  INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_INT3D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_INT3D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_INT3D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_INT3D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_INT3D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_INT3D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_INT3D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_INT3D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_INT3D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_INT3D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_INT3D_DELETE_DEVICE
END TYPE FIELD_INT3D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_INT4D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  INTEGER(KIND=JPIM), POINTER :: VIEW(:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  INTEGER(KIND=JPIM), POINTER :: PTR(:,:,:,:) => NULL()
  INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_INT4D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_INT4D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_INT4D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_INT4D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_INT4D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_INT4D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_INT4D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_INT4D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_INT4D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_INT4D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_INT4D_DELETE_DEVICE
END TYPE FIELD_INT4D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_INT5D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  INTEGER(KIND=JPIM), POINTER :: VIEW(:,:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  INTEGER(KIND=JPIM), POINTER :: PTR(:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM), ALLOCATABLE :: DATA(:,:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_INT5D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_INT5D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_INT5D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_INT5D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_INT5D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_INT5D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_INT5D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_INT5D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_INT5D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_INT5D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_INT5D_DELETE_DEVICE
END TYPE FIELD_INT5D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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
  LOGICAL, ALLOCATABLE :: DATA(:,:)

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

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_LOG3D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  LOGICAL, POINTER :: VIEW(:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  LOGICAL, POINTER :: PTR(:,:,:) => NULL()
  LOGICAL, ALLOCATABLE :: DATA(:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  LOGICAL, POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_LOG3D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_LOG3D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_LOG3D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_LOG3D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_LOG3D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_LOG3D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_LOG3D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_LOG3D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_LOG3D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_LOG3D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_LOG3D_DELETE_DEVICE
END TYPE FIELD_LOG3D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_LOG4D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  LOGICAL, POINTER :: VIEW(:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  LOGICAL, POINTER :: PTR(:,:,:,:) => NULL()
  LOGICAL, ALLOCATABLE :: DATA(:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  LOGICAL, POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_LOG4D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_LOG4D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_LOG4D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_LOG4D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_LOG4D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_LOG4D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_LOG4D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_LOG4D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_LOG4D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_LOG4D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_LOG4D_DELETE_DEVICE
END TYPE FIELD_LOG4D

# 29 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 30 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_LOG5D
  ! A FIELD encapsulates a single multi-dimensional array and can
  ! provide block-indexed "views" of the data for automating the
  ! allocation and parallel iterration of NPROMA blocks.

  ! The data view to be used in thread-parallel sections
  !
  ! The underlying view pointer is of rank-1, since we always
  ! the horizontal component as a single dimension.
  LOGICAL, POINTER :: VIEW(:,:,:,:) => NULL()

  ! TODO: Atlas-based field data storage field
  ! TODO: Do we still need to use pointers here?
  ! TYPE(ATLAS_FIELD), POINTER :: DATA

  ! Storage pointer for non-Atlas backward-compatibility mode
  !
  ! The underlying storage pointer has the rank as the dimension,
  ! where the innermost dimension represents the horizontal and
  ! the outermost one is the block index.
  LOGICAL, POINTER :: PTR(:,:,:,:,:) => NULL()
  LOGICAL, ALLOCATABLE :: DATA(:,:,:,:,:)

  ! For wrapping discontiguous fields in co-allocated storage
  ! arrays (eg. GFL/GMV) also store a CONTIGUOUS base pointer
  ! and integer index, to allow block pointer extraction that
  ! conforms with CUDA device pointers in PGI.
  LOGICAL, POINTER, CONTIGUOUS :: BASE_PTR(:,:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: FIDX

  ! A separate data pointer that can be used to create
  ! a contiguous chunk of host memory to cleanly map to
  ! device, should the %DATA pointer be discontiguous.
  LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:) => NULL()

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

  PROCEDURE :: CLONE => FIELD_LOG5D_CLONE
  PROCEDURE :: UPDATE_VIEW => FIELD_LOG5D_UPDATE_VIEW
  PROCEDURE :: EXTRACT_VIEW => FIELD_LOG5D_EXTRACT_VIEW
  PROCEDURE :: GET_VIEW => FIELD_LOG5D_GET_VIEW
  PROCEDURE :: FINAL => FIELD_LOG5D_FINAL

  ! GPU-specific device data transfer API
  PROCEDURE :: CREATE_DEVICE => FIELD_LOG5D_CREATE_DEVICE
  PROCEDURE :: UPDATE_DEVICE => FIELD_LOG5D_UPDATE_DEVICE
  PROCEDURE :: UPDATE_HOST => FIELD_LOG5D_UPDATE_HOST
  PROCEDURE :: ENSURE_DEVICE => FIELD_LOG5D_ENSURE_DEVICE
  PROCEDURE :: ENSURE_HOST => FIELD_LOG5D_ENSURE_HOST
  PROCEDURE :: DELETE_DEVICE => FIELD_LOG5D_DELETE_DEVICE
END TYPE FIELD_LOG5D

# 95 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 97 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_2D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_2D), POINTER :: PTR => NULL()
END TYPE FIELD_2D_PTR

TYPE FIELD_2D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:) => NULL()
END TYPE FIELD_2D_VIEW

# 97 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_3D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_3D), POINTER :: PTR => NULL()
END TYPE FIELD_3D_PTR

TYPE FIELD_3D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:,:) => NULL()
END TYPE FIELD_3D_VIEW

# 97 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_4D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_4D), POINTER :: PTR => NULL()
END TYPE FIELD_4D_PTR

TYPE FIELD_4D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:,:,:) => NULL()
END TYPE FIELD_4D_VIEW

# 97 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
TYPE FIELD_5D_PTR
  ! Struct to hold references to field objects
  TYPE(FIELD_5D), POINTER :: PTR => NULL()
END TYPE FIELD_5D_PTR

TYPE FIELD_5D_VIEW
  ! Struct to hold array views, so we can make arrays of them
  REAL(KIND=JPRB), POINTER :: P(:,:,:,:) => NULL()
END TYPE FIELD_5D_VIEW

# 108 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_2D
  MODULE PROCEDURE :: FIELD_2D_WRAP
  MODULE PROCEDURE :: FIELD_2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_2D_EMPTY
  MODULE PROCEDURE :: FIELD_2D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_3D
  MODULE PROCEDURE :: FIELD_3D_WRAP
  MODULE PROCEDURE :: FIELD_3D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_3D_EMPTY
  MODULE PROCEDURE :: FIELD_3D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_4D
  MODULE PROCEDURE :: FIELD_4D_WRAP
  MODULE PROCEDURE :: FIELD_4D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_4D_EMPTY
  MODULE PROCEDURE :: FIELD_4D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_5D
  MODULE PROCEDURE :: FIELD_5D_WRAP
  MODULE PROCEDURE :: FIELD_5D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_5D_EMPTY
  MODULE PROCEDURE :: FIELD_5D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_INT2D
  MODULE PROCEDURE :: FIELD_INT2D_WRAP
  MODULE PROCEDURE :: FIELD_INT2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_INT2D_EMPTY
  MODULE PROCEDURE :: FIELD_INT2D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_INT3D
  MODULE PROCEDURE :: FIELD_INT3D_WRAP
  MODULE PROCEDURE :: FIELD_INT3D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_INT3D_EMPTY
  MODULE PROCEDURE :: FIELD_INT3D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_INT4D
  MODULE PROCEDURE :: FIELD_INT4D_WRAP
  MODULE PROCEDURE :: FIELD_INT4D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_INT4D_EMPTY
  MODULE PROCEDURE :: FIELD_INT4D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_INT5D
  MODULE PROCEDURE :: FIELD_INT5D_WRAP
  MODULE PROCEDURE :: FIELD_INT5D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_INT5D_EMPTY
  MODULE PROCEDURE :: FIELD_INT5D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_LOG2D
  MODULE PROCEDURE :: FIELD_LOG2D_WRAP
  MODULE PROCEDURE :: FIELD_LOG2D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_LOG2D_EMPTY
  MODULE PROCEDURE :: FIELD_LOG2D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_LOG3D
  MODULE PROCEDURE :: FIELD_LOG3D_WRAP
  MODULE PROCEDURE :: FIELD_LOG3D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_LOG3D_EMPTY
  MODULE PROCEDURE :: FIELD_LOG3D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_LOG4D
  MODULE PROCEDURE :: FIELD_LOG4D_WRAP
  MODULE PROCEDURE :: FIELD_LOG4D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_LOG4D_EMPTY
  MODULE PROCEDURE :: FIELD_LOG4D_ALLOCATE
END INTERFACE

# 110 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 111 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
INTERFACE FIELD_LOG5D
  MODULE PROCEDURE :: FIELD_LOG5D_WRAP
  MODULE PROCEDURE :: FIELD_LOG5D_WRAP_PACKED
  ! MODULE PROCEDURE :: FIELD_LOG5D_EMPTY
  MODULE PROCEDURE :: FIELD_LOG5D_ALLOCATE
END INTERFACE

# 119 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

INTERFACE FILL_BUFFER
# 122 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FILL_BUFFER_2D, FILL_BUFFER_3D
  MODULE PROCEDURE :: FILL_BUFFER_4D, FILL_BUFFER_5D
# 122 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FILL_BUFFER_INT2D, FILL_BUFFER_INT3D
  MODULE PROCEDURE :: FILL_BUFFER_INT4D, FILL_BUFFER_INT5D
# 122 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FILL_BUFFER_LOG2D, FILL_BUFFER_LOG3D
  MODULE PROCEDURE :: FILL_BUFFER_LOG4D, FILL_BUFFER_LOG5D
# 125 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FILL_BUFFER


INTERFACE FIELD_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_CREATE_DEVICE
# 130 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 131 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_CREATE_DEVICE
# 133 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_CREATE_DEVICE

INTERFACE FIELD_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_UPDATE_DEVICE
# 137 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 138 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_UPDATE_DEVICE
# 140 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_UPDATE_DEVICE

INTERFACE FIELD_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_UPDATE_HOST
# 144 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 145 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_UPDATE_HOST
# 147 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_UPDATE_HOST

INTERFACE FIELD_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_DELETE_DEVICE
# 151 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 152 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_DELETE_DEVICE
# 154 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_DELETE_DEVICE

INTERFACE GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_GET_DEVICE_DATA
# 158 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 159 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_GET_DEVICE_DATA
# 161 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE GET_DEVICE_DATA

INTERFACE FIELD_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_ENSURE_DEVICE
# 165 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 166 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_ENSURE_DEVICE
# 168 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_ENSURE_DEVICE

INTERFACE FIELD_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_2D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_3D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_4D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_5D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT2D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT3D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT4D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_INT5D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG2D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG3D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG4D_ENSURE_HOST
# 172 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 173 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  MODULE PROCEDURE :: FIELD_LOG5D_ENSURE_HOST
# 175 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
END INTERFACE FIELD_ENSURE_HOST

CONTAINS

# 180 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

  SUBROUTINE FILL_BUFFER_5D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: BUFFER(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, K, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        DO K=1, SIZE(BUFFER, 4)
          BUFFER(IDX:,I,J,K) = BUFFER(INDEX,I,J,K)
        END DO
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_5D

# 180 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FILL_BUFFER_INT2D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: BUFFER(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: IDX

    IDX = INDEX+1
    BUFFER(IDX:) = BUFFER(INDEX)
  END SUBROUTINE FILL_BUFFER_INT2D

  SUBROUTINE FILL_BUFFER_INT3D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: BUFFER(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      BUFFER(IDX:,I) = BUFFER(INDEX,I)
    END DO
  END SUBROUTINE FILL_BUFFER_INT3D

  SUBROUTINE FILL_BUFFER_INT4D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: BUFFER(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        BUFFER(IDX:,I,J) = BUFFER(INDEX,I,J)
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_INT4D

  SUBROUTINE FILL_BUFFER_INT5D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: BUFFER(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, K, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        DO K=1, SIZE(BUFFER, 4)
          BUFFER(IDX:,I,J,K) = BUFFER(INDEX,I,J,K)
        END DO
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_INT5D

# 180 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FILL_BUFFER_LOG2D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    LOGICAL, POINTER, INTENT(INOUT) :: BUFFER(:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: IDX

    IDX = INDEX+1
    BUFFER(IDX:) = BUFFER(INDEX)
  END SUBROUTINE FILL_BUFFER_LOG2D

  SUBROUTINE FILL_BUFFER_LOG3D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    LOGICAL, POINTER, INTENT(INOUT) :: BUFFER(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      BUFFER(IDX:,I) = BUFFER(INDEX,I)
    END DO
  END SUBROUTINE FILL_BUFFER_LOG3D

  SUBROUTINE FILL_BUFFER_LOG4D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    LOGICAL, POINTER, INTENT(INOUT) :: BUFFER(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        BUFFER(IDX:,I,J) = BUFFER(INDEX,I,J)
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_LOG4D

  SUBROUTINE FILL_BUFFER_LOG5D(BUFFER, INDEX)
    ! Utility routine to fill data buffers (views)
    LOGICAL, POINTER, INTENT(INOUT) :: BUFFER(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: INDEX
    INTEGER(KIND=JPIM) :: I, J, K, IDX

    IDX = INDEX+1
    DO I=1, SIZE(BUFFER, 2)
      DO J=1, SIZE(BUFFER, 3)
        DO K=1, SIZE(BUFFER, 4)
          BUFFER(IDX:,I,J,K) = BUFFER(INDEX,I,J,K)
        END DO
      END DO
    END DO
  END SUBROUTINE FILL_BUFFER_LOG5D

# 233 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_5D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(4)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_5D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_INT3D) :: SELF
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
  END FUNCTION FIELD_INT3D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_INT4D) :: SELF
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
  END FUNCTION FIELD_INT4D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_INT5D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(4)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT5D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_LOG3D) :: SELF
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
  END FUNCTION FIELD_LOG3D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_LOG4D) :: SELF
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
  END FUNCTION FIELD_LOG4D_EMPTY

# 235 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 236 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_EMPTY(SHAPE) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    !
    ! If a SHAPE is provided, a single empty buffer block-sized buffer
    ! will be allocated under %VIEW and used by all threads in a
    ! thread-parallel region to avoid segfault when dereferencing NULL
    ! pointers. Otherwise %DATA and %VIEW will always be unassociated.
    TYPE(FIELD_LOG5D) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: SHAPE(4)

    SELF%PTR => NULL()
    IF (PRESENT(SHAPE)) THEN
      ALLOCATE(SELF%VIEW(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4)))
    END IF
    SELF%ACTIVE = .FALSE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = 0
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG5D_EMPTY

# 259 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_5D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_5D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT3D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT4D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT5D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_LOG3D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG3D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_LOG4D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG4D_WRAP

# 261 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 262 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_WRAP(DATA) RESULT(SELF)
    ! Create FIELD object by wrapping existing data
    TYPE(FIELD_LOG5D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:,:,:)

    SELF%PTR => DATA
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG5D_WRAP

# 277 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_5D), TARGET :: SELF
    REAL(KIND=JPRB), TARGET, INTENT(IN) :: DATA(:,:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_5D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_INT3D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_INT4D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), TARGET, INTENT(IN) :: DATA(:,:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_INT5D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_LOG3D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 3)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_LOG3D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_LOG4D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 4)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_LOG4D_WRAP_PACKED

# 279 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 280 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_WRAP_PACKED(DATA, IDX) RESULT(SELF)
    ! Create FIELD object packed in a multi-field buffer by storing a
    ! contiguous pointer to existing data and an index.
    TYPE(FIELD_LOG5D), TARGET :: SELF
    LOGICAL, TARGET, INTENT(IN) :: DATA(:,:,:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: IDX

    SELF%PTR => DATA(:,:,:,:,IDX,:)
    SELF%ACTIVE = .TRUE.
    SELF%THREAD_BUFFER = .FALSE.
    SELF%OWNED = .FALSE.
    SELF%NBLOCKS = SIZE(SELF%PTR, 5)
    SELF%BASE_PTR => DATA
    SELF%FIDX = IDX
  END FUNCTION FIELD_LOG5D_WRAP_PACKED

# 297 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

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
    ALLOCATE(SELF%DATA(SHAPE(1),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_2D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_3D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(2)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

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
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_3D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_4D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_4D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(3)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

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
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_4D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(4)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_5D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_5D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

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
    ALLOCATE(SELF%DATA(SHAPE(1),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT2D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(2)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_INT3D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT3D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(3)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_INT4D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT4D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(4)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_INT5D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_INT5D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG2D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_LOG2D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(1)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

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
    ALLOCATE(SELF%DATA(SHAPE(1),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 2)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG2D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_LOG3D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(2)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_LOG3D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 3)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG3D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_LOG4D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(3)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_LOG4D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 4)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG4D_ALLOCATE

# 299 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 300 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_ALLOCATE(SHAPE, NBLOCKS, PERSISTENT) RESULT(SELF)
    ! Create FIELD object by explicitly allocating new data
    !
    ! Please note that SHAPE is the conceptual shape without the block dimension
    TYPE(FIELD_LOG5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: SHAPE(4)
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NBLOCKS
    LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
    INTEGER(KIND=JPIM) :: NBLK

    ! By default we allocate thread-local temporaries
    SELF%THREAD_BUFFER = .TRUE.
    NBLK = OML_MAX_THREADS()

    IF (PRESENT(PERSISTENT)) THEN
      IF (PERSISTENT) THEN
        ! Adjust outer dim for full-sized persistent blocked arrays
        IF (.NOT. PRESENT(NBLOCKS)) CALL &
         & ABOR1('FIELD_LOG5D_ALLOCATE : NBLOCKS not given for persistent allocation!')
        SELF%THREAD_BUFFER = .FALSE.
        NBLK = NBLOCKS
      END IF
    END IF

    ! Allocate storage array and store metadata
    ALLOCATE(SELF%DATA(SHAPE(1),SHAPE(2),SHAPE(3),SHAPE(4),NBLK))
    SELF%PTR => SELF%DATA
    SELF%ACTIVE = .TRUE.
    SELF%OWNED = .TRUE.
    SELF%NBLOCKS = SIZE(SELF%DATA, 5)
    SELF%BASE_PTR => NULL()
    SELF%FIDX = -1
  END FUNCTION FIELD_LOG5D_ALLOCATE

# 335 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_5D) :: SELF
    TYPE(FIELD_5D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_5D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_INT3D) :: SELF
    TYPE(FIELD_INT3D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_INT3D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_INT4D) :: SELF
    TYPE(FIELD_INT4D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_INT4D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_INT5D) :: SELF
    TYPE(FIELD_INT5D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_INT5D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_LOG3D) :: SELF
    TYPE(FIELD_LOG3D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_LOG3D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_LOG4D) :: SELF
    TYPE(FIELD_LOG4D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_LOG4D_CLONE

# 337 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 338 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this FIELD object, keeping the DATA pointer
    ! intact, but replicating view pointers.
    CLASS(FIELD_LOG5D) :: SELF
    TYPE(FIELD_LOG5D), POINTER :: NEWOBJ

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
  END FUNCTION FIELD_LOG5D_CLONE

# 359 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,:,:,:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:,:,:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_5D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_INT3D), TARGET :: SELF
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
      IF (ZERO) SELF%VIEW(:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT3D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_INT4D), TARGET :: SELF
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
      IF (ZERO) SELF%VIEW(:,:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT4D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,:,:,:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:,:,:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT5D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_LOG3D), TARGET :: SELF
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
      IF (ZERO) SELF%VIEW(:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG3D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_LOG4D), TARGET :: SELF
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
      IF (ZERO) SELF%VIEW(:,:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG4D_UPDATE_VIEW

# 361 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 362 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_UPDATE_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Sets the view pointer FIELD%MP to the block of the given index
    CLASS(FIELD_LOG5D), TARGET :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      SELF%VIEW => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      SELF%VIEW => SELF%PTR(:,:,:,:,IDX)
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(SELF%VIEW, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) SELF%VIEW(:,:,:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG5D_UPDATE_VIEW

# 389 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_5D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:,:,:) = 0.0_JPRB
    END IF
  END SUBROUTINE FIELD_5D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:)
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
      IF (ZERO) VIEW_PTR(:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT3D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:)
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
      IF (ZERO) VIEW_PTR(:,:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT4D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:,:,:) = 0.0_JPIM
    END IF
  END SUBROUTINE FIELD_INT5D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG3D), TARGET :: SELF
    LOGICAL, POINTER, INTENT(INOUT) :: VIEW_PTR(:,:)
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
      IF (ZERO) VIEW_PTR(:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG3D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG4D), TARGET :: SELF
    LOGICAL, POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:)
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
      IF (ZERO) VIEW_PTR(:,:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG4D_EXTRACT_VIEW

# 391 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 392 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_EXTRACT_VIEW(SELF, VIEW_PTR, BLOCK_INDEX, BLOCK_SIZE, ZERO)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG5D), TARGET :: SELF
    LOGICAL, POINTER, INTENT(INOUT) :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
    ELSE
      VIEW_PTR => SELF%VIEW  ! Set to NaN'd field buffer
    END IF

    IF (PRESENT(BLOCK_SIZE) .AND. BLOCK_INDEX == SELF%NBLOCKS) THEN
      ! Fill the the buffer by replicating the last entry
      CALL FILL_BUFFER(VIEW_PTR, INDEX=BLOCK_SIZE)
    END IF

    IF (PRESENT(ZERO)) THEN
      IF (ZERO) VIEW_PTR(:,:,:,:) = .FALSE.
    END IF
  END SUBROUTINE FIELD_LOG5D_EXTRACT_VIEW

# 422 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_5D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
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
      IF (ZERO) VIEW_PTR(:,:,:,:) = 0.0_JPRB
    END IF
  END FUNCTION FIELD_5D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER :: VIEW_PTR(:,:)
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
      IF (ZERO) VIEW_PTR(:,:) = 0.0_JPIM
    END IF
  END FUNCTION FIELD_INT3D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER :: VIEW_PTR(:,:,:)
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
      IF (ZERO) VIEW_PTR(:,:,:) = 0.0_JPIM
    END IF
  END FUNCTION FIELD_INT4D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
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
      IF (ZERO) VIEW_PTR(:,:,:,:) = 0.0_JPIM
    END IF
  END FUNCTION FIELD_INT5D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG3D), TARGET :: SELF
    LOGICAL, POINTER :: VIEW_PTR(:,:)
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
      IF (ZERO) VIEW_PTR(:,:) = .FALSE.
    END IF
  END FUNCTION FIELD_LOG3D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG4D), TARGET :: SELF
    LOGICAL, POINTER :: VIEW_PTR(:,:,:)
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
      IF (ZERO) VIEW_PTR(:,:,:) = .FALSE.
    END IF
  END FUNCTION FIELD_LOG4D_GET_VIEW

# 424 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 425 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_GET_VIEW(SELF, BLOCK_INDEX, BLOCK_SIZE, ZERO) RESULT(VIEW_PTR)
    ! Updates internal view and exports it to an external pointer
    CLASS(FIELD_LOG5D), TARGET :: SELF
    LOGICAL, POINTER :: VIEW_PTR(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: BLOCK_SIZE
    LOGICAL, OPTIONAL, INTENT(IN) :: ZERO
    INTEGER(KIND=JPIM) :: IDX

    IDX = BLOCK_INDEX
    IF (SELF%THREAD_BUFFER) IDX = OML_MY_THREAD()
    IF (SELF%ACTIVE .AND. SELF%OWNED) THEN
      VIEW_PTR => SELF%DATA(:,:,:,:,IDX)
    ELSEIF (SELF%ACTIVE .AND. .NOT. SELF%OWNED) THEN
      VIEW_PTR => SELF%PTR(:,:,:,:,IDX)
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
      IF (ZERO) VIEW_PTR(:,:,:,:) = .FALSE.
    END IF
  END FUNCTION FIELD_LOG5D_GET_VIEW

# 457 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_2D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_3D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_4D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_4D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_4D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_5D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_5D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT2D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT3D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT3D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT4D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT4D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT5D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT5D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG2D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG2D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG2D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG3D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG3D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG4D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG4D_CREATE_DEVICE

# 459 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 460 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_CREATE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG5D), TARGET :: SELF

    SELF%DEVPTR => SELF%DATA
    !$acc enter data create(SELF%DATA)
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG5D_CREATE_DEVICE

# 470 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_3D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_3D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_5D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_5D), TARGET :: SELF
    REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_5D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT3D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_INT3D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT4D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_INT4D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_INT5D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM), POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_INT5D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG3D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG3D), TARGET :: SELF
    LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_LOG3D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG4D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG4D), TARGET :: SELF
    LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_LOG4D_GET_DEVICE_DATA

# 472 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 473 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  FUNCTION FIELD_LOG5D_GET_DEVICE_DATA(SELF) RESULT(DEVPTR)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG5D), TARGET :: SELF
    LOGICAL, POINTER, CONTIGUOUS :: DEVPTR(:,:,:,:,:)

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF

    IF (SELF%OWNED) THEN
      DEVPTR => SELF%DATA
    ELSE
      DEVPTR => SELF%DEVPTR
    END IF
  END FUNCTION FIELD_LOG5D_GET_DEVICE_DATA

# 490 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,IBL))
      END DO
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

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_3D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_3D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_4D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_4D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,IBL))
      END DO
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

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_5D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_5D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_INT2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,IBL))
      END DO
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

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_INT3D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT3D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_INT4D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,IBL))
      END DO
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
  END SUBROUTINE FIELD_INT4D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_INT5D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_INT5D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG2D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_LOG2D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,IBL))
      END DO
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

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_LOG3D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG3D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_LOG4D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,IBL))
      END DO
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
  END SUBROUTINE FIELD_LOG4D_UPDATE_DEVICE

# 492 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 493 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_UPDATE_DEVICE(SELF)
    ! Create a copy of this field on device and copy data over
    CLASS(FIELD_LOG5D), TARGET :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      !$acc enter data create(SELF%DATA)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      SELF%DEVPTR => SELF%DATA
    ELSE
      ALLOCATE(SELF%DEVPTR, SOURCE=SELF%PTR)
      !$acc enter data create(SELF%DEVPTR)
      DO IBL=1, SELF%NBLOCKS
        !$acc update device(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
    END IF
    SELF%ON_DEVICE = .TRUE.
  END SUBROUTINE FIELD_LOG5D_UPDATE_DEVICE

# 517 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_2D) :: SELF
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
  END SUBROUTINE FIELD_2D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_3D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:) = SELF%DEVPTR(:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_3D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_5D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:,:,:) = SELF%DEVPTR(:,:,:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_5D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_INT3D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:) = SELF%DEVPTR(:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT3D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_INT4D) :: SELF
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
  END SUBROUTINE FIELD_INT4D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_INT5D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:,:,:) = SELF%DEVPTR(:,:,:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT5D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_LOG3D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:) = SELF%DEVPTR(:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG3D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_LOG4D) :: SELF
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
  END SUBROUTINE FIELD_LOG4D_UPDATE_HOST

# 519 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 520 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_UPDATE_HOST(SELF)
    ! Synchronize device data back to host
    CLASS(FIELD_LOG5D) :: SELF
    INTEGER(KIND=JPIM) :: IBL

    IF (SELF%OWNED) THEN
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DATA(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DATA)
    ELSE
      DO IBL=1, SELF%NBLOCKS
        !$acc update host(SELF%DEVPTR(:,:,:,:,IBL))
      END DO
      !$acc wait
      !$acc exit data delete(SELF%DEVPTR)
      SELF%PTR(:,:,:,:,:) = SELF%DEVPTR(:,:,:,:,:)
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG5D_UPDATE_HOST

# 544 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
 
# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_3D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_3D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_5D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_5D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT3D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT3D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT4D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT4D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_INT5D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_INT5D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
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

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG3D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG3D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG4D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG4D_DELETE_DEVICE

# 546 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 547 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_DELETE_DEVICE(SELF)
    ! Initialize a copy of this field on GPU device
    CLASS(FIELD_LOG5D), TARGET :: SELF

    !$acc exit data delete(SELF%DEVPTR)
    IF (SELF%OWNED) THEN
      NULLIFY(SELF%DEVPTR)
    ELSE
      DEALLOCATE(SELF%DEVPTR)
    END IF
    SELF%ON_DEVICE = .FALSE.
  END SUBROUTINE FIELD_LOG5D_DELETE_DEVICE

# 561 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_2D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_3D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_3D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_4D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_4D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_4D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_5D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_5D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_INT2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_INT2D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_INT3D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_INT3D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_INT4D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_INT4D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_INT5D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_INT5D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG2D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_LOG2D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_LOG2D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_LOG3D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_LOG3D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_LOG4D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_LOG4D_ENSURE_HOST

# 563 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 564 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_ENSURE_HOST(SELF)
    ! Ensure that field has been moved back to host
    CLASS(FIELD_LOG5D), TARGET :: SELF

    IF (SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_HOST()
    END IF
  END SUBROUTINE FIELD_LOG5D_ENSURE_HOST

# 574 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_2D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_3D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_3D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_4D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_4D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_4D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_5D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_5D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_INT2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_INT2D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_INT3D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_INT3D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_INT4D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_INT4D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_INT5D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_INT5D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG2D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_LOG2D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_LOG2D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_LOG3D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_LOG3D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_LOG4D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_LOG4D_ENSURE_DEVICE

# 576 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 577 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_ENSURE_DEVICE(SELF)
    ! Ensure that field has been moved over to device
    CLASS(FIELD_LOG5D), TARGET :: SELF

    IF (.NOT. SELF%ON_DEVICE) THEN
      CALL SELF%UPDATE_DEVICE()
    END IF
  END SUBROUTINE FIELD_LOG5D_ENSURE_DEVICE

# 587 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_2D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_3D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_3D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_3D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_4D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_4D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_4D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_5D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_5D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_5D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_INT2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_INT2D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT3D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_INT3D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_INT3D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT4D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_INT4D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_INT4D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_INT5D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_INT5D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_INT5D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG2D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_LOG2D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_LOG2D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG3D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_LOG3D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_LOG3D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG4D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_LOG4D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_LOG4D_FINAL

# 589 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
# 590 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"
  SUBROUTINE FIELD_LOG5D_FINAL(SELF)
    ! Finalizes field and dealloactes owned data
    CLASS(FIELD_LOG5D) :: SELF
    IF (SELF%OWNED) THEN
      DEALLOCATE(SELF%DATA)
    END IF
    NULLIFY(SELF%PTR)
    NULLIFY(SELF%VIEW)
  END SUBROUTINE FIELD_LOG5D_FINAL

# 601 "/local/hdd/naml/ifs-bundle-gpu/source/ifs_dp/ifs/module/field_module.fypp"

END MODULE FIELD_MODULE
