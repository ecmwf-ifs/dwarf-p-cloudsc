module expand_mod
  use parkind1 , only : jpim, jprb
  use yomphyder, only : state_type

  use cloudsc_mpi_mod, only : irank, numproc
  use file_io_mod, only: input_initialize, load_scalar, load_array

  implicit none

  interface expand
     procedure expand_l1, expand_i1, expand_r1, expand_r2, expand_r3
  end interface expand

  interface load_and_expand
     procedure load_and_expand_l1, load_and_expand_i1
     procedure load_and_expand_r1, load_and_expand_r2, load_and_expand_r3
  end interface load_and_expand

contains

  subroutine get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    integer(kind=jpim), intent(inout) :: start, end, size
    integer(kind=jpim), intent(in) :: nlon, ndim, nlev, ngptot
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim) :: rankstride
    logical :: use_offset = .false.

    if (present(ngptotg)) use_offset = nlon >= ngptotg
    if (use_offset) then
      rankstride = (ngptotg - 1) / numproc + 1
      start = irank * rankstride + 1
    else
      start = 1
    end if
    size = min(nlon, ngptot)
    end = start + size - 1
  end subroutine get_offsets

  subroutine load_and_expand_i1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    integer(kind=jpim), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim), allocatable :: buffer(:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_i1

  subroutine load_and_expand_l1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    logical, pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    logical, allocatable :: buffer(:), rbuf(:)
    integer(kind=jpim) :: start, end, size
    integer(kind=4), allocatable :: tmp(:)

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_l1

  subroutine load_and_expand_r1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r1

  subroutine load_and_expand_r2(name, field, nlon, nlev, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, nlev, ngptot, ngptotg)
    allocate(field(nproma, nlev, nblocks))
    allocate(buffer(size, nlev))
    call load_array(name, start, end, size, nlon, nlev, buffer)
    call expand(buffer, field, size, nproma, nlev, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r2

  subroutine load_and_expand_r3(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    allocate(field(nproma, nlev, ndim, nblocks))
    allocate(buffer(size, nlev, ndim))
    call load_array(name, start, end, size, nlon, nlev, ndim, buffer)
    call expand(buffer, field, size, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r3

  subroutine load_and_expand_state(name, state, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    type(state_type), pointer, intent(inout) :: state(:)
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start, end, size

    integer :: b

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    allocate(state(nblocks))
    allocate(field(nproma, nlev, 6+ndim, nblocks))
    allocate(buffer(size, nlev, 6+ndim))

    call load_array(name//'_T', start, end, size, nlon, nlev, buffer(:,:,3))
    call load_array(name//'_A', start, end, size, nlon, nlev, buffer(:,:,5))
    call load_array(name//'_Q', start, end, size, nlon, nlev, buffer(:,:,6))
    call load_array(name//'_CLD', start, end, size, nlon, nlev, ndim, buffer(:,:,7:))

    call expand(buffer(:,:,3), field(:,:,3,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,5), field(:,:,5,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,6), field(:,:,6,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,7:), field(:,:,7:,:), size, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B) schedule(runtime)
    do b=1, nblocks
       ! state(b)%u => field(:,:,1,b)
       ! state(b)%v => field(:,:,2,b)
       state(b)%t => field(:,:,3,b)
       ! state(b)%o3 => field(:,:,4,b)
       state(b)%a => field(:,:,5,b)
       state(b)%q => field(:,:,6,b)
       state(b)%cld => field(:,:,7:6+ndim,b)
    end do
!$OMP end parallel do

  end subroutine load_and_expand_state

  subroutine expand_l1(buffer, field, nlon, nproma, ngptot, nblocks)
    logical, intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
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
!$omp end parallel do
  end subroutine expand_l1

  subroutine expand_i1(buffer, field, nlon, nproma, ngptot, nblocks)
    integer(kind=jpim), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime) 
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
!$omp end parallel do
  end subroutine expand_i1

  subroutine expand_r1(buffer, field, nlon, nproma, ngptot, nblocks)
    real(kind=jprb), intent(inout) :: buffer(nlon)
    real(kind=jprb), intent(inout) :: field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime) 
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
!$omp end parallel do    
  end subroutine expand_r1

  subroutine expand_r2(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
          use omp_lib
    real(kind=jprb), intent(inout) :: buffer(nlon, nlev)
    real(kind=jprb), intent(inout) :: field(nproma, nlev, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
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
!$omp end parallel do

  end subroutine expand_r2

  subroutine expand_r3(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    real(kind=jprb), intent(inout) :: buffer(nlon, nlev, ndim)
    real(kind=jprb), intent(inout) :: field(nproma, nlev, ndim, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime) 
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
!$omp end parallel do
  end subroutine expand_r3

end module expand_mod
