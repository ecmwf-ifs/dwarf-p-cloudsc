module cloudsc_mpi_mod
  ! Helper module to handle MPI communication and environment

  use parkind1, only : jpim, jprd, jprs

#ifdef HAVE_MPI
  use cloudsc_mpif
#endif

  implicit none

  integer(kind=jpim) :: numproc = 1  ! number of MPI processes
  integer(kind=jpim) :: irank   = 0  ! local process number

  interface cloudsc_mpi_reduce
    procedure cloudsc_mpi_reduce_jprd, cloudsc_mpi_reduce_jprs, cloudsc_mpi_reduce_int
  end interface cloudsc_mpi_reduce

  interface cloudsc_mpi_reduce_sum
#ifdef HAVE_MPI
    procedure cloudsc_mpi_reduce_sum_jprd, cloudsc_mpi_reduce_sum_jprs, cloudsc_mpi_reduce_sum_int
#else
    procedure cloudsc_mpi_reduce_dummy_jprd, cloudsc_mpi_reduce_dummy_jprs, cloudsc_mpi_reduce_dummy_int
#endif
  end interface cloudsc_mpi_reduce_sum

  interface cloudsc_mpi_reduce_min
#ifdef HAVE_MPI
    procedure cloudsc_mpi_reduce_min_jprd, cloudsc_mpi_reduce_min_jprs, cloudsc_mpi_reduce_min_int
#else
    procedure cloudsc_mpi_reduce_dummy_jprd, cloudsc_mpi_reduce_dummy_jprs, cloudsc_mpi_reduce_dummy_int
#endif
  end interface cloudsc_mpi_reduce_min

  interface cloudsc_mpi_reduce_max
#ifdef HAVE_MPI
    procedure cloudsc_mpi_reduce_max_jprd, cloudsc_mpi_reduce_max_jprs, cloudsc_mpi_reduce_max_int
#else
    procedure cloudsc_mpi_reduce_dummy_jprd, cloudsc_mpi_reduce_dummy_jprs, cloudsc_mpi_reduce_dummy_int
#endif
  end interface cloudsc_mpi_reduce_max

  interface cloudsc_mpi_gather
    procedure cloudsc_mpi_gather_jprd, cloudsc_mpi_gather_jprs, cloudsc_mpi_gather_int
  end interface cloudsc_mpi_gather

contains

  subroutine cloudsc_mpi_init(numomp)
    ! cloudsc_mpi_init - initialises the MPI environment

    integer(kind=jpim), intent(in), optional :: numomp   ! number of OpenMP threads
#ifdef HAVE_MPI
    integer(kind=jpim) :: ierror, iprovided, irequired  ! MPI status variables

    ! request threading support if multiple OpenMP threads are used
    iprovided = mpi_thread_single
    irequired = mpi_thread_single
    if (present(numomp)) then
      if (numomp > 1) then
        irequired = mpi_thread_multiple
      end if
    end if

    call mpi_init_thread(irequired, iprovided, ierror)

    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_init_thread failed')
    if (iprovided < irequired) then
      call abor1('cloudsc_mpi: insufficient threading support provided')
    end if

    ! determine communicator size and local rank
    call mpi_comm_rank(MPI_COMM_WORLD, irank, ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_comm_rank failed')

    call mpi_comm_size(MPI_COMM_WORLD, numproc, ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_comm_size failed')
#endif
  end subroutine cloudsc_mpi_init

  subroutine cloudsc_mpi_end()
    ! cloudsc_mpi_end - terminates the MPI environment

#ifdef HAVE_MPI
    integer(kind=jpim) :: ierror

    call mpi_finalize(ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_finalize failed')
#endif

  end subroutine cloudsc_mpi_end

  subroutine cloudsc_mpi_reduce_jprd(buf, nelem, op, iroot)
    ! cloudsc_mpi_reduce_jprd - implementation of cloudsc_mpi_reduce for double

    real(kind=jprd)   , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem
    integer           , intent(in)    :: op
    integer(kind=jpim), intent(in)    :: iroot
#ifdef HAVE_MPI
    integer(kind=jpim)                :: ierror
    logical 			      :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    if (irank == iroot) then
      call mpi_reduce(MPI_IN_PLACE, buf, nelem, &
        & MPI_DOUBLE, op, iroot, MPI_COMM_WORLD, ierror)
    else
      call mpi_reduce(buf, 0, nelem, &
        & MPI_DOUBLE, op, iroot, MPI_COMM_WORLD, ierror)
    end if
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_reduce failed')
#endif
  end subroutine cloudsc_mpi_reduce_jprd

  subroutine cloudsc_mpi_reduce_jprs(buf, nelem, op, iroot)
    ! cloudsc_mpi_reduce_jprs - implementation of cloudsc_mpi_reduce for single

    real(kind=jprs)   , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem
    integer           , intent(in)    :: op
    integer(kind=jpim), intent(in)    :: iroot
#ifdef HAVE_MPI
    integer(kind=jpim)                :: ierror
    logical 			      :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    if (irank == iroot) then
      call mpi_reduce(MPI_IN_PLACE, buf, nelem, &
        & MPI_FLOAT, op, iroot, MPI_COMM_WORLD, ierror)
    else
      call mpi_reduce(buf, 0, nelem, &
        & MPI_FLOAT, op, iroot, MPI_COMM_WORLD, ierror)
    end if
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_reduce failed')
#endif
  end subroutine cloudsc_mpi_reduce_jprs

  subroutine cloudsc_mpi_reduce_int(buf, nelem, op, iroot)
    ! cloudsc_mpi_reduce_int - implementation of cloudsc_mpi_reduce for int

    integer           , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem
    integer           , intent(in)    :: op
    integer(kind=jpim), intent(in)    :: iroot
#ifdef HAVE_MPI
    integer(kind=jpim)                :: ierror
    logical 			      :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    if (irank == iroot) then
      call mpi_reduce(MPI_IN_PLACE, buf, nelem, &
        & MPI_INTEGER, op, iroot, MPI_COMM_WORLD, ierror)
    else
      call mpi_reduce(buf, 0, nelem, &
        & MPI_INTEGER, op, iroot, MPI_COMM_WORLD, ierror)
    end if
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_reduce failed')
#endif
  end subroutine cloudsc_mpi_reduce_int

#ifdef HAVE_MPI

  subroutine cloudsc_mpi_reduce_sum_jprd(buf, nelem, iroot)
    real(kind=jprd), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_SUM, iroot)
  end subroutine cloudsc_mpi_reduce_sum_jprd

  subroutine cloudsc_mpi_reduce_sum_jprs(buf, nelem, iroot)
    real(kind=jprs), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_SUM, iroot)
  end subroutine cloudsc_mpi_reduce_sum_jprs

  subroutine cloudsc_mpi_reduce_sum_int(buf, nelem, iroot)
    integer           , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_SUM, iroot)
  end subroutine cloudsc_mpi_reduce_sum_int

  subroutine cloudsc_mpi_reduce_min_jprd(buf, nelem, iroot)
    real(kind=jprd), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MIN, iroot)
  end subroutine cloudsc_mpi_reduce_min_jprd

  subroutine cloudsc_mpi_reduce_min_jprs(buf, nelem, iroot)
    real(kind=jprs), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MIN, iroot)
  end subroutine cloudsc_mpi_reduce_min_jprs

  subroutine cloudsc_mpi_reduce_min_int(buf, nelem, iroot)
    integer           , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MIN, iroot)
  end subroutine cloudsc_mpi_reduce_min_int

  subroutine cloudsc_mpi_reduce_max_jprd(buf, nelem, iroot)
    real(kind=jprd), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MAX, iroot)
  end subroutine cloudsc_mpi_reduce_max_jprd

  subroutine cloudsc_mpi_reduce_max_jprs(buf, nelem, iroot)
    real(kind=jprs), intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MAX, iroot)
  end subroutine cloudsc_mpi_reduce_max_jprs

  subroutine cloudsc_mpi_reduce_max_int(buf, nelem, iroot)
    integer           , intent(inout) :: buf(*)
    integer(kind=jpim), intent(in) :: nelem, iroot

    call cloudsc_mpi_reduce(buf, nelem, MPI_MAX, iroot)
  end subroutine cloudsc_mpi_reduce_max_int

#else

  subroutine cloudsc_mpi_reduce_dummy_jprd(buf, nelem, iroot)
    real(kind=jprd),    intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem, iroot
  end subroutine cloudsc_mpi_reduce_dummy_jprd

  subroutine cloudsc_mpi_reduce_dummy_jprs(buf, nelem, iroot)
    real(kind=jprs),    intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem, iroot
  end subroutine cloudsc_mpi_reduce_dummy_jprs

  subroutine cloudsc_mpi_reduce_dummy_int(buf, nelem, iroot)
    integer,            intent(inout) :: buf(*)
    integer(kind=jpim), intent(in)    :: nelem, iroot
  end subroutine cloudsc_mpi_reduce_dummy_int

#endif

  subroutine cloudsc_mpi_gather_jprd(sendbuf, sendcount, recvbuf, recvcount, iroot)
    real(kind=jprd),    intent(in)  :: sendbuf(:,:)
    real(kind=jprd),    intent(out) :: recvbuf(:,:,:)
    integer(kind=jpim), intent(in)  :: sendcount, recvcount, iroot

#ifdef HAVE_MPI
    integer(kind=jpim) :: ierror
    logical :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    call mpi_gather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, iroot, MPI_COMM_WORLD, ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_gather failed')
#else
    recvbuf(:,:,1) = sendbuf(:,:)
#endif
  end subroutine cloudsc_mpi_gather_jprd

  subroutine cloudsc_mpi_gather_jprs(sendbuf, sendcount, recvbuf, recvcount, iroot)
    real(kind=jprs),    intent(in)  :: sendbuf(:,:)
    real(kind=jprs),    intent(out) :: recvbuf(:,:,:)
    integer(kind=jpim), intent(in)  :: sendcount, recvcount, iroot

#ifdef HAVE_MPI
    integer(kind=jpim) :: ierror
    logical :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    call mpi_gather(sendbuf, sendcount, MPI_FLOAT, recvbuf, recvcount, MPI_FLOAT, iroot, MPI_COMM_WORLD, ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_gather failed')
#else
    recvbuf(:,:,1) = sendbuf(:,:)
#endif
  end subroutine cloudsc_mpi_gather_jprs

  subroutine cloudsc_mpi_gather_int(sendbuf, sendcount, recvbuf, recvcount, iroot)
    integer,            intent(in)  :: sendbuf(:,:)
    integer,            intent(out) :: recvbuf(:,:,:)
    integer(kind=jpim), intent(in)  :: sendcount, recvcount, iroot

#ifdef HAVE_MPI
    integer(kind=jpim) :: ierror
    logical :: is_initialized

    ! bail out early in case mpi was never initialized
    call mpi_initialized(is_initialized, ierror)
    if (.not. is_initialized) return

    call mpi_gather(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, iroot, MPI_COMM_WORLD, ierror)
    if (ierror /= 0) call abor1('cloudsc_mpi: mpi_gather failed')
#else
    recvbuf(:,:,1) = sendbuf(:,:)
#endif
  end subroutine cloudsc_mpi_gather_int

end module cloudsc_mpi_mod
