MODULE hdf5_file_mod

  USE PARKIND1 , ONLY : JPIM, JPRB, JPRD

  USE hdf5

  implicit none

  TYPE hdf5_library
    LOGICAL :: is_init = .false.
    LOGICAL :: is_closed = .false.

  CONTAINS
    procedure :: init => hdf5_library_init
    final :: hdf5_library_close
  END TYPE hdf5_library

  TYPE(hdf5_library) hdf5_library_status

  TYPE hdf5_file

    LOGICAL :: is_rdwr = .false.
    LOGICAL :: is_open = .false.
    INTEGER(HID_T) :: file_id

  CONTAINS

    procedure :: create_file => hdf5_file_create
    procedure :: open_file => hdf5_file_open 
    procedure :: assert_file_open => hdf5_file_assert_open
    procedure :: close_file => hdf5_file_close

    procedure, private :: hdf5_file_load_l0, hdf5_file_load_i0, hdf5_file_load_r0
    procedure, private :: hdf5_file_load_l1, hdf5_file_load_i1, hdf5_file_load_r1
    procedure, private :: hdf5_file_load_r2, hdf5_file_load_r3

    generic :: load => hdf5_file_load_l0, hdf5_file_load_i0, hdf5_file_load_r0, &
                        & hdf5_file_load_l1, hdf5_file_load_i1, hdf5_file_load_r1, &
                        & hdf5_file_load_r2, hdf5_file_load_r3

  END TYPE hdf5_file

CONTAINS

  subroutine hdf5_library_init(self)
    class(hdf5_library) :: self
    integer :: error

    if (.not. self%is_init) then
      call h5open_f(error)
      if (error /= 0) call abor1('ERROR: HDF5 Fortran interface not initialized.')

      self%is_init = .true.
      self%is_closed = .false.
    end if
  end subroutine hdf5_library_init

  subroutine hdf5_library_close(self)
    type(hdf5_library) :: self
    integer :: error

    if (.not. self%is_closed .and. self%is_init) then
      call h5close_f(error)
      if (error /= 0) call abor1('ERROR: HDF5 Fortran interface not closed.')
      self%is_init = .false.
      self%is_closed = .true.
    end if
  end subroutine hdf5_library_close

  subroutine hdf5_file_create(self, filename)
    class(hdf5_file) :: self
    character(len=*), intent(in) :: filename
    integer :: error

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, self%file_id, error)
    if (error /= 0) call abor1('ERROR: Failed to create HDF5 file.')
  end subroutine hdf5_file_create

  subroutine hdf5_file_open(self, filename)
    class(hdf5_file) :: self
    character(len=*), intent(in) :: filename
    integer :: error

    if (.not. self%is_open) then
      call hdf5_library_status%init()

      if (self%is_rdwr) then
        call self%create_file(filename)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, self%file_id, error)
      else
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%file_id, error)
      end if
      if (error /= 0) call abor1('ERROR: Failed to open HDF5 file.')

      self%is_open = .true.
    end if
  end subroutine hdf5_file_open

  subroutine hdf5_file_assert_open(self)
    class(hdf5_file) :: self

    if (.not. self%is_open) then
      call abor1('ERROR: HDF5 file is not open.')
    end if
  end subroutine hdf5_file_assert_open

  subroutine hdf5_file_close(self)
    class(hdf5_file) :: self
    integer :: error
    
    if (self%is_open) then
      call h5fclose_f(self%file_id, error)
      if (error /= 0) call abor1('ERROR: Failed to close HDF5 input file.')

      self%is_open = .false.
    end if
  end subroutine hdf5_file_close

  subroutine hdf5_dataset_load(file_id, type_id, name, buf, start, count)
    use iso_c_binding, only: c_ptr

    integer(HID_T), intent(in) :: file_id, type_id
    character(len=*), intent(in) :: name
    type(c_ptr), intent(inout) :: buf
    integer(HSIZE_T), intent(in), optional :: start(:), count(:)
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: error

    call h5dopen_f(file_id, name, dset_id, error)
    if (error /= 0) call abor1('ERROR: Failed to open dataset.')

    if (present(start) .and. present(count)) then
      call h5dget_space_f(dset_id, file_space_id, error)
      if (error /= 0) call abor1('Error: Failed to get file dataspace.')
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, start, count, error)
      if (error /= 0) call abor1('Error: Failed to select hyperslab in file dataspace.')
      call h5screate_simple_f(size(count), count, mem_space_id, error) 
      if (error /= 0) call abor1('Error: Failed to create memory dataspace.')
      call h5dread_f(dset_id, type_id, buf, error, mem_space_id, file_space_id)
      if (error /= 0) call abor1('ERROR: Failed to read dataset.')
    else
      call h5dread_f(dset_id, type_id, buf, error)
      if (error /= 0) call abor1('ERROR: Failed to read dataset.')
    end if

    call h5dclose_f(dset_id, error)
    if (error /= 0) call abor1('ERROR: Failed to close dataset.')
  end subroutine hdf5_dataset_load

  subroutine hdf5_file_load_l0(self, name, val)
    use iso_c_binding, only: c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    logical, intent(inout) :: val
    integer(kind=jpim), target :: buffer
    type(c_ptr) :: ptr
    ptr = c_loc(buffer)

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, h5kind_to_type(jpim, H5_INTEGER_KIND), name, ptr)
    if (buffer > 0) then
      val = .true.
    else
      val = .false.
    end if
  end subroutine hdf5_file_load_l0

  subroutine hdf5_file_load_i0(self, name, val)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    integer(kind=jpim), target, intent(inout) :: val
    type(c_ptr) :: ptr
    ptr = c_loc(val)

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, h5kind_to_type(jpim, H5_INTEGER_KIND), name, ptr)
  end subroutine hdf5_file_load_i0

  subroutine hdf5_file_load_r0(self, name, val)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    real(kind=jprd), intent(inout) :: val
    real(kind=jprd), target :: buf
    type(c_ptr) :: ptr
    ptr = c_loc(buf)

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, H5T_NATIVE_DOUBLE, name, ptr)
    val = buf
  end subroutine hdf5_file_load_r0

  subroutine hdf5_file_load_l1(self, name, buffer, start, count)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    logical, intent(inout) :: buffer(:)
    integer(kind=jpim), intent(in), optional :: start, count
    integer(kind=hsize_t) :: hstart(1) = (/ 0 /), hcount(1)
    integer(kind=jpim), target :: tmp(size(buffer))
    type(c_ptr) :: ptr

    if (present(start)) then
      hstart(1) = start
    end if
    if (present(count)) then
      hcount(1) = count
    else
      hcount(1) = size(buffer)
    end if
    ptr = c_loc(tmp(1))

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, h5kind_to_type(jpim, H5_INTEGER_KIND), name, ptr, hstart, hcount)
    where (tmp>0)
      buffer = .true.
    elsewhere
      buffer = .false.
    end where
  end subroutine hdf5_file_load_l1

  subroutine hdf5_file_load_i1(self, name, buffer, start, count)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    integer(kind=jpim), target, intent(inout) :: buffer(:)
    integer(kind=jpim), intent(in), optional :: start, count
    integer(kind=hsize_t) :: hstart(1) = (/ 0 /), hcount(1)
    type(c_ptr) :: ptr

    if (present(start)) then
      hstart(1) = start
    end if
    if (present(count)) then
      hcount(1) = count
    else
      hcount(1) = size(buffer)
    end if
    ptr = c_loc(buffer(1))

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, h5kind_to_type(jpim, H5_INTEGER_KIND), name, ptr, hstart, hcount)
  end subroutine hdf5_file_load_i1

  subroutine hdf5_file_load_r1(self, name, buffer, start, count)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    real(kind=jprd), target, intent(inout) :: buffer(:)
    integer(kind=jpim), intent(in), optional :: start, count
    integer(kind=hsize_t) :: hstart(1), hcount(1)
    type(c_ptr) :: ptr

    if (present(start)) then
      hstart(1) = start
    else
      hstart(1) = 0
    end if
    if (present(count)) then
      hcount(1) = count
    else
      hcount(1) = size(buffer)
    end if
    ptr = c_loc(buffer(1))

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, H5T_NATIVE_DOUBLE, name, ptr, hstart, hcount)
  end subroutine hdf5_file_load_r1

  subroutine hdf5_file_load_r2(self, name, buffer, start, count)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    real(kind=jprd), target, intent(inout) :: buffer(:,:)
    integer(kind=jpim), intent(in), optional :: start(2), count(2)
    integer(kind=hsize_t) :: hstart(2), hcount(2)
    type(c_ptr) :: ptr

    if (present(start)) then
      hstart(:) = start(:)
    else
      hstart(:) = 0
    end if
    if (present(count)) then
      hcount(:) = count(:)
    else
      hcount = shape(buffer, kind=hsize_t)
    end if
    ptr = c_loc(buffer(1,1))

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, H5T_NATIVE_DOUBLE, name, ptr, hstart, hcount)
  end subroutine hdf5_file_load_r2

  subroutine hdf5_file_load_r3(self, name, buffer, start, count)
    use iso_c_binding, only: c_ptr, c_loc

    class(hdf5_file) :: self
    character(len=*), intent(in) :: name
    real(kind=jprd), target, intent(inout) :: buffer(:,:,:)
    integer(kind=jpim), intent(in), optional :: start(3), count(3)
    integer(kind=hsize_t) :: hstart(3), hcount(3)
    type(c_ptr) :: ptr

    if (present(start)) then
      hstart(:) = start(:)
    else
      hstart(:) = 0
    end if
    if (present(count)) then
      hcount(:) = count(:)
    else
      hcount = shape(buffer, kind=hsize_t)
    end if
    ptr = c_loc(buffer(1,1,1))

    call self%assert_file_open()
    call hdf5_dataset_load(self%file_id, H5T_NATIVE_DOUBLE, name, ptr, hstart, hcount)
  end subroutine hdf5_file_load_r3

END MODULE hdf5_file_mod
