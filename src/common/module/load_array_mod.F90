module load_array_mod
  USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
  USE YOMPHYDER, ONLY : STATE_TYPE

  USE YOECLDP  , ONLY : YRECLDP, TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV
  USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
   & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
   & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
  USE YOEPHLI  , ONLY : YREPHLI, TEPHLI
  use cloudsc_mpi_mod, only : irank, numproc

#ifdef HAVE_SERIALBOX
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
#endif

#ifdef HAVE_HDF5
  USE hdf5_file_mod
#endif

  implicit none

#ifdef HAVE_HDF5
  TYPE(hdf5_file) :: input_file
#endif

  interface expand
     procedure expand_l1, expand_i1, expand_r1, expand_r2, expand_r3
  end interface expand

  interface load_and_expand
     procedure load_and_expand_l1, load_and_expand_i1
     procedure load_and_expand_r1, load_and_expand_r2, load_and_expand_r3
  end interface load_and_expand

  interface load_metainfo
     procedure load_metainfo_scalar_real, load_metainfo_scalar_int, load_metainfo_scalar_log
     procedure load_metainfo_field
  end interface load_metainfo

contains

  subroutine query_dimensions(KLON, KLEV, KFLDX, NAME)
    ! Initial query routine to determine data dimensions
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLON             ! Number of grid points
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLEV             ! Number of levels
    ! INTEGER(KIND=JPIM),INTENT(OUT) :: NCLV
    INTEGER(KIND=JPIM),INTENT(OUT) :: KFLDX
    CHARACTER(*), INTENT(IN) :: NAME

#ifdef HAVE_SERIALBOX
    ! Get dimensions information from stored metadata
    call ppser_initialize(directory='data', prefix='dummy', prefix_ref=NAME)
    call fs_create_savepoint(NAME, ppser_savepoint)
    call ppser_set_mode(1)

    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLON', KLON)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLEV', KLEV)
    ! call fs_get_serializer_metainfo(ppser_serializer_ref, 'NCLV', NCLV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KFLDX', KFLDX)
#elif defined(HAVE_HDF5)
    call input_file%open_file(NAME//'.h5')
    call input_file%load('KLON', KLON)
    call input_file%load('KLEV', KLEV)
    call input_file%load('KFLDX', KFLDX)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif

  end subroutine query_dimensions

  subroutine get_offsets(start, count, nlon, ndim, nlev, ngptot, ngptotg)
    integer(kind=jpim), intent(inout) :: start(3), count(3)
    integer(kind=jpim), intent(in) :: nlon, ndim, nlev, ngptot
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim) :: rankstride
    logical :: use_offset = .false.

    if (present(ngptotg)) use_offset = nlon >= ngptotg
    if (use_offset) then
      rankstride = (ngptotg - 1) / numproc + 1
      start(1) = irank * rankstride
    else
      start(1) = 0
    end if
    start(2) = 0
    start(3) = 0
    count(1) = min(nlon, ngptot)
    count(2) = nlev
    count(3) = ndim
  end subroutine get_offsets

  subroutine load_and_expand_i1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    integer(kind=jpim), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim), allocatable :: buffer(:)
    integer(kind=jpim) :: start(3), count(3)

    allocate(field(nproma, nblocks))
#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, 1, 1, ngptot, ngptotg)
    allocate(buffer(count(1)))
    call input_file%load(name, buffer, start(1), count(1))
    call expand(buffer, field, count(1), nproma, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_and_expand_i1

  subroutine load_and_expand_l1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    logical, pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    logical, allocatable :: buffer(:)
    integer(kind=jpim) :: start(3), count(3)

    allocate(field(nproma, nblocks))
#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, 1, 1, ngptot, ngptotg)
    allocate(buffer(count(1)))
    call input_file%load(name, buffer, start(1), count(1))
    call expand(buffer, field, count(1), nproma, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_and_expand_l1

  subroutine load_and_expand_r1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprd), allocatable :: buffer(:)
    integer(kind=jpim) :: start(3), count(3)

    allocate(field(nproma, nblocks))
#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, ngptot, nblocks)
    deallocate(buffer)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, 1, 1, ngptot, ngptotg)
    allocate(buffer(count(1)))
    call input_file%load(name, buffer, start(1), count(1))
    call expand(buffer, field, count(1), nproma, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_and_expand_r1

  subroutine load_and_expand_r2(name, field, nlon, nlev, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprd), allocatable :: buffer(:,:)
    integer(kind=jpim) :: start(3), count(3)

    allocate(field(nproma, nlev, nblocks))
#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon, nlev))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
    deallocate(buffer)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, 1, nlev, ngptot, ngptotg)
    allocate(buffer(count(1), count(2)))
    call input_file%load(name, buffer, start(1:2), count(1:2))
    call expand(buffer, field, count(1), nproma, nlev, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_and_expand_r2

  subroutine load_and_expand_r3(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprd), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start(3), count(3)

    allocate(field(nproma, nlev, ndim, nblocks))
#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon, nlev, ndim))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    call expand(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, ndim, nlev, ngptot, ngptotg)
    allocate(buffer(count(1), count(2), count(3)))
    call input_file%load(name, buffer, start, count)
    call expand(buffer, field, count(1), nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_and_expand_r3

  subroutine load_and_expand_state(name, state, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    type(state_type), pointer, intent(inout) :: state(:)
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprd), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start(3), count(3)

    integer :: b

    allocate(state(nblocks))
    allocate(field(nproma, nlev, 6+ndim, nblocks))

#ifdef HAVE_SERIALBOX
    allocate(buffer(nlon, nlev, 6+ndim))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_U', buffer(:,:,1))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_V', buffer(:,:,2))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_T', buffer(:,:,3))
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_O3', buffer(:,:,4))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_A', buffer(:,:,5))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_Q', buffer(:,:,6))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name//'_CLD', buffer(:,:,7:))

    ! call expand(buffer(:,:,1), field(:,:,1,:), nlon, nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,2), field(:,:,2,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field(:,:,3,:), nlon, nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,4), field(:,:,4,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,5), field(:,:,5,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,6), field(:,:,6,:), nlon, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,7:), field(:,:,7:,:), nlon, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
#elif defined(HAVE_HDF5)
    call get_offsets(start, count, nlon, ndim, nlev, ngptot, ngptotg)
    allocate(buffer(count(1), count(2), 6+ndim))
    ! call input_file%load(name//'_U', buffer(:,:,1), start(1:2), count(1:2))
    ! call input_file%load(name//'_V', buffer(:,:,2), start(1:2), count(1:2))
    call input_file%load(name//'_T', buffer(:,:,3), start(1:2), count(1:2))
    ! call input_file%load(name//'_O3', buffer(:,:,4), start(1:2), count(1:2))
    call input_file%load(name//'_A', buffer(:,:,5), start(1:2), count(1:2))
    call input_file%load(name//'_Q', buffer(:,:,6), start(1:2), count(1:2))
    call input_file%load(name//'_CLD', buffer(:,:,7:), start, count)

    ! call expand(buffer(:,:,1), field(:,:,1,:), count(1), nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,2), field(:,:,2,:), count(1), nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field(:,:,3,:), count(1), nproma, nlev, ngptot, nblocks)
    ! call expand(buffer(:,:,4), field(:,:,4,:), count(1), nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,5), field(:,:,5,:), count(1), nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,6), field(:,:,6,:), count(1), nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,7:), field(:,:,7:,:), count(1), nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif

    !OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B)
    do b=1, nblocks
       ! state(b)%u => field(:,:,1,b)
       ! state(b)%v => field(:,:,2,b)
       state(b)%t => field(:,:,3,b)
       ! state(b)%o3 => field(:,:,4,b)
       state(b)%a => field(:,:,5,b)
       state(b)%q => field(:,:,6,b)
       state(b)%cld => field(:,:,7:6+ndim,b)
    end do

  end subroutine load_and_expand_state

  subroutine expand_l1(buffer, field, nlon, nproma, ngptot, nblocks)
    logical, intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

    !omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx,nlon)
       bend = min(nlon, mod(gidx, nlon)+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = .FALSE.
    end do
  end subroutine expand_l1

  subroutine expand_i1(buffer, field, nlon, nproma, ngptot, nblocks)
    integer(kind=jpim), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

    !omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx,nlon)
       bend = min(nlon, mod(gidx, nlon)+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0_JPIM
    end do
  end subroutine expand_i1

  subroutine expand_r1(buffer, field, nlon, nproma, ngptot, nblocks)
    real(kind=JPRD), intent(inout) :: buffer(nlon)
    real(kind=JPRB), intent(inout) :: field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

    !omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx,nlon)
       bend = min(nlon, mod(gidx, nlon)+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0.0_JPRB
    end do
  end subroutine expand_r1

  subroutine expand_r2(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
    real(kind=JPRD), intent(inout) :: buffer(nlon, nlev)
    real(kind=JPRB), intent(inout) :: field(nproma, nlev, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

    !omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx,nlon)
       bend = min(nlon, mod(gidx, nlon)+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,:,b) = buffer(bidx:bend,:)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,:,b) = buffer(bidx:bend,:)
       end do

       field(bsize+1:nproma,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r2

  subroutine expand_r3(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    real(kind=JPRD), intent(inout) :: buffer(nlon, nlev, ndim)
    real(kind=JPRB), intent(inout) :: field(nproma, nlev, ndim, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

    !omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx,nlon)
       bend = min(nlon, mod(gidx, nlon)+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,:,:,b) = buffer(bidx:bend,:,:)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,:,:,b) = buffer(bidx:bend,:,:)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r3

  subroutine load_metainfo_scalar_real(name, variable)
    character(len=*), intent(in) :: name
    real(kind=JPRB), intent(inout) :: variable
    real(kind=JPRD) :: buffer

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, buffer)
    variable = buffer
#elif defined(HAVE_HDF5)
    call input_file%load(name, buffer)
    variable = buffer
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_metainfo_scalar_real

  subroutine load_metainfo_scalar_int(name, variable)
    character(len=*), intent(in) :: name
    integer(kind=JPIM), intent(inout) :: variable

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, variable)
#elif defined(HAVE_HDF5)
    call input_file%load(name, variable)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_metainfo_scalar_int

  subroutine load_metainfo_scalar_log(name, variable)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: variable

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, variable)
#elif defined(HAVE_HDF5)
    call input_file%load(name, variable)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_metainfo_scalar_log

  subroutine load_metainfo_field(name, variable)
    character(len=*), intent(in) :: name
    real(kind=JPRB), intent(inout) :: variable(:)
    real(kind=JPRD), allocatable :: buffer(:)

#ifdef HAVE_SERIALBOX
    allocate(buffer(size(variable)))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, buffer)
    variable(:) = buffer(:)
    deallocate(buffer)
#elif defined(HAVE_HDF5)
    allocate(buffer(size(variable)))
    call input_file%load(name, buffer)
    variable(:) = buffer(:)
    deallocate(buffer)
#else
    call abor1('ERROR: Serialbox not found.')
#endif
  end subroutine load_metainfo_field

  subroutine initialise_parameters(PTSPHY, LDSLPHY, LDMAINCALL)
    ! Retrieve parametersand timestep size from the serializer
    ! Note that we use `ppser_serializer_ref` to get the previously stored data
    real(kind=JPRB), intent(inout) :: PTSPHY
    logical, intent(inout) :: LDSLPHY, LDMAINCALL

#if defined(HAVE_SERIALBOX) || defined(HAVE_HDF5)
    call load_metainfo('PTSPHY', PTSPHY)
#ifdef HAVE_SERIALBOX
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'LDSLPHY', LDSLPHY)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'LDMAINCALL', LDMAINCALL)
#else
    call input_file%load('LDSLPHY', LDSLPHY)
    call input_file%load('LDMAINCALL', LDMAINCALL)
#endif

    call load_metainfo('RG', RG)
    call load_metainfo('RD', RD)
    call load_metainfo('RCPD', RCPD)
    call load_metainfo('RETV', RETV)
    call load_metainfo('RLVTT', RLVTT)
    call load_metainfo('RLSTT', RLSTT)
    call load_metainfo('RLMLT', RLMLT)
    call load_metainfo('RTT', RTT)
    call load_metainfo('RV', RV)
    call load_metainfo('R2ES', R2ES)
    call load_metainfo('R3LES', R3LES)
    call load_metainfo('R3IES', R3IES)
    call load_metainfo('R4LES', R4LES)
    call load_metainfo('R4IES', R4IES)
    call load_metainfo('R5LES', R5LES)
    call load_metainfo('R5IES', R5IES)
    call load_metainfo('R5ALVCP', R5ALVCP)
    call load_metainfo('R5ALSCP', R5ALSCP)
    call load_metainfo('RALVDCP', RALVDCP)
    call load_metainfo('RALSDCP', RALSDCP)
    call load_metainfo('RALFDCP', RALFDCP)
    call load_metainfo('RTWAT', RTWAT)
    call load_metainfo('RTICE', RTICE)
    call load_metainfo('RTICECU', RTICECU)
    call load_metainfo('RTWAT_RTICE_R', RTWAT_RTICE_R)
    call load_metainfo('RTWAT_RTICECU_R', RTWAT_RTICECU_R)
    call load_metainfo('RKOOP1', RKOOP1)
    call load_metainfo('RKOOP2', RKOOP2)

    ! Retrieve parameters contained in TECLDP type
    allocate(YRECLDP)
    call load_metainfo('YRECLDP_RAMID', YRECLDP%RAMID)
    call load_metainfo('YRECLDP_RCLDIFF', YRECLDP%RCLDIFF)
    call load_metainfo('YRECLDP_RCLDIFF_CONVI', YRECLDP%RCLDIFF_CONVI)
    call load_metainfo('YRECLDP_RCLCRIT', YRECLDP%RCLCRIT)
    call load_metainfo('YRECLDP_RCLCRIT_SEA', YRECLDP%RCLCRIT_SEA)
    call load_metainfo('YRECLDP_RCLCRIT_LAND', YRECLDP%RCLCRIT_LAND)
    call load_metainfo('YRECLDP_RKCONV', YRECLDP%RKCONV)
    call load_metainfo('YRECLDP_RPRC1', YRECLDP%RPRC1)
    call load_metainfo('YRECLDP_RPRC2', YRECLDP%RPRC2)
    call load_metainfo('YRECLDP_RCLDMAX', YRECLDP%RCLDMAX)
    call load_metainfo('YRECLDP_RPECONS', YRECLDP%RPECONS)
    call load_metainfo('YRECLDP_RVRFACTOR', YRECLDP%RVRFACTOR)
    call load_metainfo('YRECLDP_RPRECRHMAX', YRECLDP%RPRECRHMAX)
    call load_metainfo('YRECLDP_RTAUMEL', YRECLDP%RTAUMEL)
    call load_metainfo('YRECLDP_RAMIN', YRECLDP%RAMIN)
    call load_metainfo('YRECLDP_RLMIN', YRECLDP%RLMIN)
    call load_metainfo('YRECLDP_RKOOPTAU', YRECLDP%RKOOPTAU)

    call load_metainfo('YRECLDP_RCLDTOPP', YRECLDP%RCLDTOPP)
    call load_metainfo('YRECLDP_RLCRITSNOW', YRECLDP%RLCRITSNOW)
    call load_metainfo('YRECLDP_RSNOWLIN1', YRECLDP%RSNOWLIN1)
    call load_metainfo('YRECLDP_RSNOWLIN2', YRECLDP%RSNOWLIN2)
    call load_metainfo('YRECLDP_RICEHI1', YRECLDP%RICEHI1)
    call load_metainfo('YRECLDP_RICEHI2', YRECLDP%RICEHI2)
    call load_metainfo('YRECLDP_RICEINIT', YRECLDP%RICEINIT)
    call load_metainfo('YRECLDP_RVICE', YRECLDP%RVICE)
    call load_metainfo('YRECLDP_RVRAIN', YRECLDP%RVRAIN)
    call load_metainfo('YRECLDP_RVSNOW', YRECLDP%RVSNOW)
    call load_metainfo('YRECLDP_RTHOMO', YRECLDP%RTHOMO)
    call load_metainfo('YRECLDP_RCOVPMIN', YRECLDP%RCOVPMIN)
    call load_metainfo('YRECLDP_RCCN', YRECLDP%RCCN)
    call load_metainfo('YRECLDP_RNICE', YRECLDP%RNICE)
    call load_metainfo('YRECLDP_RCCNOM', YRECLDP%RCCNOM)
    call load_metainfo('YRECLDP_RCCNSS', YRECLDP%RCCNSS)
    call load_metainfo('YRECLDP_RCCNSU', YRECLDP%RCCNSU)
    call load_metainfo('YRECLDP_RCLDTOPCF', YRECLDP%RCLDTOPCF)
    call load_metainfo('YRECLDP_RDEPLIQREFRATE', YRECLDP%RDEPLIQREFRATE)
    call load_metainfo('YRECLDP_RDEPLIQREFDEPTH', YRECLDP%RDEPLIQREFDEPTH)
    Call load_metainfo('YRECLDP_RCL_KKAac', YRECLDP%RCL_KKAac)
    Call load_metainfo('YRECLDP_RCL_KKBac', YRECLDP%RCL_KKBac)
    Call load_metainfo('YRECLDP_RCL_KKAau', YRECLDP%RCL_KKAau)
    Call load_metainfo('YRECLDP_RCL_KKBauq', YRECLDP%RCL_KKBauq)
    Call load_metainfo('YRECLDP_RCL_KKBaun', YRECLDP%RCL_KKBaun)
    Call load_metainfo('YRECLDP_RCL_KK_cloud_num_sea', YRECLDP%RCL_KK_cloud_num_sea)
    Call load_metainfo('YRECLDP_RCL_KK_cloud_num_land', YRECLDP%RCL_KK_cloud_num_land)
    call load_metainfo('YRECLDP_RCL_AI', YRECLDP%RCL_AI)
    call load_metainfo('YRECLDP_RCL_BI', YRECLDP%RCL_BI)
    call load_metainfo('YRECLDP_RCL_CI', YRECLDP%RCL_CI)
    call load_metainfo('YRECLDP_RCL_DI', YRECLDP%RCL_DI)
    call load_metainfo('YRECLDP_RCL_X1I', YRECLDP%RCL_X1I)
    call load_metainfo('YRECLDP_RCL_X2I', YRECLDP%RCL_X2I)
    call load_metainfo('YRECLDP_RCL_X3I', YRECLDP%RCL_X3I)
    call load_metainfo('YRECLDP_RCL_X4I', YRECLDP%RCL_X4I)
    call load_metainfo('YRECLDP_RCL_CONST1I', YRECLDP%RCL_CONST1I)
    call load_metainfo('YRECLDP_RCL_CONST2I', YRECLDP%RCL_CONST2I)
    call load_metainfo('YRECLDP_RCL_CONST3I', YRECLDP%RCL_CONST3I)
    call load_metainfo('YRECLDP_RCL_CONST4I', YRECLDP%RCL_CONST4I)
    call load_metainfo('YRECLDP_RCL_CONST5I', YRECLDP%RCL_CONST5I)
    call load_metainfo('YRECLDP_RCL_CONST6I', YRECLDP%RCL_CONST6I)
    call load_metainfo('YRECLDP_RCL_APB1', YRECLDP%RCL_APB1)
    call load_metainfo('YRECLDP_RCL_APB2', YRECLDP%RCL_APB2)
    call load_metainfo('YRECLDP_RCL_APB3', YRECLDP%RCL_APB3)
    call load_metainfo('YRECLDP_RCL_AS', YRECLDP%RCL_AS)
    call load_metainfo('YRECLDP_RCL_BS', YRECLDP%RCL_BS)
    call load_metainfo('YRECLDP_RCL_CS', YRECLDP%RCL_CS)
    call load_metainfo('YRECLDP_RCL_DS', YRECLDP%RCL_DS)
    call load_metainfo('YRECLDP_RCL_X1S', YRECLDP%RCL_X1S)
    call load_metainfo('YRECLDP_RCL_X2S', YRECLDP%RCL_X2S)
    call load_metainfo('YRECLDP_RCL_X3S', YRECLDP%RCL_X3S)
    call load_metainfo('YRECLDP_RCL_X4S', YRECLDP%RCL_X4S)
    call load_metainfo('YRECLDP_RCL_CONST1S', YRECLDP%RCL_CONST1S)
    call load_metainfo('YRECLDP_RCL_CONST2S', YRECLDP%RCL_CONST2S)
    call load_metainfo('YRECLDP_RCL_CONST3S', YRECLDP%RCL_CONST3S)
    call load_metainfo('YRECLDP_RCL_CONST4S', YRECLDP%RCL_CONST4S)
    call load_metainfo('YRECLDP_RCL_CONST5S', YRECLDP%RCL_CONST5S)
    call load_metainfo('YRECLDP_RCL_CONST6S', YRECLDP%RCL_CONST6S)
    call load_metainfo('YRECLDP_RCL_CONST7S', YRECLDP%RCL_CONST7S)
    call load_metainfo('YRECLDP_RCL_CONST8S', YRECLDP%RCL_CONST8S)
    call load_metainfo('YRECLDP_RDENSWAT', YRECLDP%RDENSWAT)
    call load_metainfo('YRECLDP_RDENSREF', YRECLDP%RDENSREF)
    call load_metainfo('YRECLDP_RCL_AR', YRECLDP%RCL_AR)
    call load_metainfo('YRECLDP_RCL_BR', YRECLDP%RCL_BR)
    call load_metainfo('YRECLDP_RCL_CR', YRECLDP%RCL_CR)
    call load_metainfo('YRECLDP_RCL_DR', YRECLDP%RCL_DR)
    call load_metainfo('YRECLDP_RCL_X1R', YRECLDP%RCL_X1R)
    call load_metainfo('YRECLDP_RCL_X2R', YRECLDP%RCL_X2R)
    call load_metainfo('YRECLDP_RCL_X4R', YRECLDP%RCL_X4R)
    call load_metainfo('YRECLDP_RCL_KA273', YRECLDP%RCL_KA273)
    call load_metainfo('YRECLDP_RCL_CDENOM1', YRECLDP%RCL_CDENOM1)
    call load_metainfo('YRECLDP_RCL_CDENOM2', YRECLDP%RCL_CDENOM2)
    call load_metainfo('YRECLDP_RCL_CDENOM3', YRECLDP%RCL_CDENOM3)
    call load_metainfo('YRECLDP_RCL_SCHMIDT', YRECLDP%RCL_SCHMIDT)
    call load_metainfo('YRECLDP_RCL_DYNVISC', YRECLDP%RCL_DYNVISC)
    call load_metainfo('YRECLDP_RCL_CONST1R', YRECLDP%RCL_CONST1R)
    call load_metainfo('YRECLDP_RCL_CONST2R', YRECLDP%RCL_CONST2R)
    call load_metainfo('YRECLDP_RCL_CONST3R', YRECLDP%RCL_CONST3R)
    call load_metainfo('YRECLDP_RCL_CONST4R', YRECLDP%RCL_CONST4R)
    call load_metainfo('YRECLDP_RCL_FAC1', YRECLDP%RCL_FAC1)
    call load_metainfo('YRECLDP_RCL_FAC2', YRECLDP%RCL_FAC2)
    call load_metainfo('YRECLDP_RCL_CONST5R', YRECLDP%RCL_CONST5R)
    call load_metainfo('YRECLDP_RCL_CONST6R', YRECLDP%RCL_CONST6R)
    call load_metainfo('YRECLDP_RCL_FZRAB', YRECLDP%RCL_FZRAB)
    call load_metainfo('YRECLDP_RCL_FZRBB', YRECLDP%RCL_FZRBB)
    call load_metainfo('YRECLDP_LCLDEXTRA', YRECLDP%LCLDEXTRA)
    call load_metainfo('YRECLDP_LCLDBUDGET', YRECLDP%LCLDBUDGET)
    call load_metainfo('YRECLDP_NSSOPT', YRECLDP%NSSOPT)
    call load_metainfo('YRECLDP_NCLDTOP', YRECLDP%NCLDTOP)
    call load_metainfo('YRECLDP_NAECLBC', YRECLDP%NAECLBC)
    call load_metainfo('YRECLDP_NAECLDU', YRECLDP%NAECLDU)
    call load_metainfo('YRECLDP_NAECLOM', YRECLDP%NAECLOM)
    call load_metainfo('YRECLDP_NAECLSS', YRECLDP%NAECLSS)
    call load_metainfo('YRECLDP_NAECLSU', YRECLDP%NAECLSU)
    call load_metainfo('YRECLDP_NCLDDIAG', YRECLDP%NCLDDIAG)
    call load_metainfo('YRECLDP_NAERCLD', YRECLDP%NAERCLD)
    call load_metainfo('YRECLDP_LAERLIQAUTOLSP', YRECLDP%LAERLIQAUTOLSP)
    call load_metainfo('YRECLDP_LAERLIQAUTOCP', YRECLDP%LAERLIQAUTOCP)
    call load_metainfo('YRECLDP_LAERLIQAUTOCPB', YRECLDP%LAERLIQAUTOCPB)
    call load_metainfo('YRECLDP_LAERLIQCOLL', YRECLDP%LAERLIQCOLL)
    call load_metainfo('YRECLDP_LAERICESED', YRECLDP%LAERICESED)
    call load_metainfo('YRECLDP_LAERICEAUTO', YRECLDP%LAERICEAUTO)
    call load_metainfo('YRECLDP_NSHAPEP', YRECLDP%NSHAPEP)
    call load_metainfo('YRECLDP_NSHAPEQ', YRECLDP%NSHAPEQ)
    call load_metainfo('YRECLDP_NBETA', YRECLDP%NBETA)
    ! The last two are actually arrays, so treat them as fields
    call load_metainfo('YRECLDP_RBETA', YRECLDP%RBETA(0:100))
    call load_metainfo('YRECLDP_RBETAP1', YRECLDP%RBETAP1(0:100))

    ! Retrieve parameters contained in TECLDP type
    allocate(YREPHLI)
    call load_metainfo('YREPHLI_LTLEVOL', YREPHLI%LTLEVOL)
    call load_metainfo('YREPHLI_LPHYLIN', YREPHLI%LPHYLIN)
    call load_metainfo('YREPHLI_LENOPERT', YREPHLI%LENOPERT)
    call load_metainfo('YREPHLI_LEPPCFLS', YREPHLI%LEPPCFLS)
    call load_metainfo('YREPHLI_LRAISANEN', YREPHLI%LRAISANEN)
    call load_metainfo('YREPHLI_RLPTRC', YREPHLI%RLPTRC)
    call load_metainfo('YREPHLI_RLPAL1', YREPHLI%RLPAL1)
    call load_metainfo('YREPHLI_RLPAL2', YREPHLI%RLPAL2)
    call load_metainfo('YREPHLI_RLPBB', YREPHLI%RLPBB)
    call load_metainfo('YREPHLI_RLPCC', YREPHLI%RLPCC)
    call load_metainfo('YREPHLI_RLPDD', YREPHLI%RLPDD)
    call load_metainfo('YREPHLI_RLPMIXL', YREPHLI%RLPMIXL)
    call load_metainfo('YREPHLI_RLPBETA', YREPHLI%RLPBETA)
    call load_metainfo('YREPHLI_RLPDRAG', YREPHLI%RLPDRAG)
    call load_metainfo('YREPHLI_RLPEVAP', YREPHLI%RLPEVAP)
    call load_metainfo('YREPHLI_RLPP00', YREPHLI%RLPP00)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine initialise_parameters

  subroutine finalize()
#ifdef HAVE_SERIALBOX
    call ppser_finalize()
#elif defined(HAVE_HDF5)
    call input_file%close_file()
#else
    call abor1('ERROR: Serialbox not found.')
#endif
  end subroutine finalize
end module load_array_mod
