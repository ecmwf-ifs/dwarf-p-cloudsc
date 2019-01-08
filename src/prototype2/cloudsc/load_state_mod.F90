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
     procedure expand_r1, expand_r2, expand_r3
  end interface expand

  interface load_and_expand
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

  subroutine load_and_expand_r1(name, field, nlon, nproma, ngptot, nblocks)
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon)

    allocate(field(nproma, nblocks))

    ! Load into the local memory buffer and expand to global field
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
  end subroutine load_and_expand_r1

  subroutine load_and_expand_r2(name, field, nlon, nlev, nproma, ngptot, nblocks)
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev)

    allocate(field(nproma, nlev, nblocks))

    ! Load into the local memory buffer and expand to global field
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
  end subroutine load_and_expand_r2

  subroutine load_and_expand_r3(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks)
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev, ndim)

    allocate(field(nproma, nlev, ndim, nblocks))

    ! Load into the local memory buffer and expand to global field
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
  end subroutine load_and_expand_r3

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
