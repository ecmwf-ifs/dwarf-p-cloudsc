module load_state_mod
  USE PARKIND1 , ONLY : JPIM, JPRB
  USE YOMPHYDER, ONLY : STATE_TYPE

  USE m_serialize, ONLY: &
   fs_create_savepoint, &
   fs_add_serializer_metainfo, &
   fs_get_serializer_metainfo, &
   fs_read_field, &
   fs_write_field
  USE utils_ppser, ONLY:  &
   ppser_initialize, &
   ppser_finalize, &
   ppser_serializer, &
   ppser_serializer_ref, &
   ppser_set_mode, &
   ppser_savepoint

  implicit none

  interface expand
     procedure expand_l1, expand_i1, expand_r1, expand_r2, expand_r3
  end interface expand

  interface load_and_expand
     procedure load_and_expand_l1, load_and_expand_i1
     procedure load_and_expand_r1, load_and_expand_r2, load_and_expand_r3
  end interface load_and_expand

contains

  subroutine query_dimensions(KLON, KLEV, KFLDX)
    ! Initial query routine to determine data dimensions
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLON             ! Number of grid points
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLEV             ! Number of levels
    ! INTEGER(KIND=JPIM),INTENT(OUT) :: NCLV
    INTEGER(KIND=JPIM),INTENT(OUT) :: KFLDX

    ! Get dimensions information from stored metadata
    call ppser_initialize(directory='data', prefix='dummy', prefix_ref='input')
    call fs_create_savepoint('input', ppser_savepoint)
    call ppser_set_mode(1)

    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLON', KLON)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLEV', KLEV)
    ! call fs_get_serializer_metainfo(ppser_serializer_ref, 'NCLV', NCLV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KFLDX', KFLDX)

  end subroutine query_dimensions

  subroutine load_and_expand_i1(name, field, nlon, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    integer(kind=jpim), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim) :: buffer(nlon)

    allocate(field(nproma, nblocks))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
  end subroutine load_and_expand_i1

  subroutine load_and_expand_l1(name, field, nlon, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    logical, pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    logical :: buffer(nlon)

    allocate(field(nproma, nblocks))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
  end subroutine load_and_expand_l1

  subroutine load_and_expand_r1(name, field, nlon, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon)

    allocate(field(nproma, nblocks))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
  end subroutine load_and_expand_r1

  subroutine load_and_expand_r2(name, field, nlon, nlev, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev)

    allocate(field(nproma, nlev, nblocks))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
  end subroutine load_and_expand_r2

  subroutine load_and_expand_r3(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev, ndim)

    allocate(field(nproma, nlev, ndim, nblocks))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
  end subroutine load_and_expand_r3

  subroutine load_and_expand_state(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev, 6+ndim)

    allocate(field(nproma, nlev, ndim, nblocks))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_U', buffer(:,:,1))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_V', buffer(:,:,2))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_T', buffer(:,:,3))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_O3', buffer(:,:,4))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_A', buffer(:,:,5))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_Q', buffer(:,:,6))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_CLD', buffer(:,:,7:))

    ! This is a beyond hacky! We basically rely on the underlyingstructure
    ! of the memory block under the the apparent STATE_TYPE array.
    ! call expand(buffer(:,:,1), field(:,:,1,:), nlon, nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,2), field(:,:,2,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field(:,:,3,:), nlon, nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,4), field(:,:,4,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,5), field(:,:,5,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,6), field(:,:,6,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,7:), field(:,:,7:,:), nlon, nproma, nlev, ndim, ngptot, nblocks)
  end subroutine load_and_expand_state

  subroutine expand_l1(buffer, field, nlon, nproma, ngptot, nblocks)
    logical, intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, bsize, bidx

    !omp parallel do default(shared) private(b, bsize, bidx)
    do b=1, nblocks
       bsize = min(nproma, ngptot - (b-1)*nproma)  ! Field block size
       bidx = mod((b-1)*nproma, nlon)  ! Rolling index into input buffer
       field(1:bsize,b) = buffer(1:bsize)

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = .FALSE.
    end do
  end subroutine expand_l1

  subroutine expand_i1(buffer, field, nlon, nproma, ngptot, nblocks)
    integer(kind=jpim), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, bsize, bidx

    !omp parallel do default(shared) private(b, bsize, bidx)
    do b=1, nblocks
       bsize = min(nproma, ngptot - (b-1)*nproma)  ! Field block size
       bidx = mod((b-1)*nproma, nlon)  ! Rolling index into input buffer
       field(1:bsize,b) = buffer(1:bsize)

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0_JPIM
    end do
  end subroutine expand_i1

  subroutine expand_r1(buffer, field, nlon, nproma, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, bsize, bidx

    !omp parallel do default(shared) private(b, bsize, bidx)
    do b=1, nblocks
       bsize = min(nproma, ngptot - (b-1)*nproma)  ! Field block size
       bidx = mod((b-1)*nproma, nlon)  ! Rolling index into input buffer
       field(1:bsize,b) = buffer(1:bsize)

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0.0_JPRB
    end do
  end subroutine expand_r1

  subroutine expand_r2(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon, nlev), field(nproma, nlev, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer :: b, bsize, bidx

    !omp parallel do default(shared) private(b, bsize, bidx)
    do b=1, nblocks
       bsize = min(nproma, ngptot - (b-1)*nproma)  ! Field block size
       bidx = mod((b-1)*nproma, nlon)  ! Rolling index into input buffer
       field(1:bsize,:,b) = buffer(1:bsize,:)

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r2

  subroutine expand_r3(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon, nlev, ndim), field(nproma, nlev, ndim, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer :: b, bsize, bidx

    !omp parallel do default(shared) private(b, bsize, bidx)
    do b=1, nblocks
       bsize = min(nproma, ngptot - (b-1)*nproma)  ! Field block size
       bidx = mod((b-1)*nproma, nlon)  ! Rolling index into input buffer
       field(1:bsize,:,:,b) = buffer(1:bsize,:,:)

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r3

end module load_state_mod
