! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module expand_atlas_mod
  use atlas_module
  use atlas_fieldset_module
  use atlas_functionspace_blockstructuredcolumns_module

  use parkind1 , only : jpim, jprb
  use yomphyder, only : state_type

  use cloudsc_mpi_mod, only : irank, numproc
  use file_io_mod, only: input_initialize, load_scalar, load_array
  use expand_mod, only: get_offsets, expand

  use, intrinsic :: iso_c_binding, only : c_int

  implicit none

contains

  subroutine loadvar_atlas(fset, fspace, name, nlon, ngptotg)
    ! Load into the local memory buffer and expand to global field
    type(atlas_fieldset), intent(inout) :: fset
    type(atlas_functionspace_blockstructuredcolumns), intent(in) :: fspace
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: nlon
    integer(kind=jpim), intent(in), optional :: ngptotg

    integer(kind=jpim) :: start, end, size, nlev, nproma, ngptot, nblocks, ndim, frank
    type(atlas_field) :: field
    real(kind=jprb), allocatable :: buffer_r1(:), buffer_r2(:,:), buffer_r3(:,:,:)
    integer(kind=jpim), allocatable :: buffer_i1(:)
    logical, allocatable :: buffer_l1(:)
    real(kind=jprb), pointer :: field_r1(:,:), field_r2(:,:,:), field_r3(:,:,:,:)
    integer(c_int), pointer :: field_i1(:,:)
    logical, pointer :: field_l1(:,:)
    logical :: lfield, rfield, ifield
    type(atlas_trace) :: trace
    trace = atlas_trace("expand_atlas_mod.F90", __LINE__, "loadvar_atlas", "IO")

    field = fset%field(name)
    frank = field%rank()
    lfield = (name == "LDCUM")
    ifield = (name == "KTYPE")
    rfield = ((.not. lfield) .and. (.not. ifield))

    nlev = field%levels()
    !nproma = fspace%nproma()
    nproma = field%shape(1)
    ngptot = fspace%size()
    nblocks = fspace%nblks()

    if (frank == 2) then
      call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
      if (rfield) then
        allocate(buffer_r1(size))
        call field%data(field_r1)
        call load_array(name, start, end, size, nlon, buffer_r1)
        call expand(buffer_r1, field_r1, size, nproma, ngptot, nblocks)
        deallocate(buffer_r1)
      else if (lfield) then
        allocate(buffer_l1(size))
        call field%data(field_l1)
        call load_array(name, start, end, size, nlon, buffer_l1)
        call expand(buffer_l1, field_l1, size, nproma, ngptot, nblocks)
        deallocate(buffer_l1)
      else
        allocate(buffer_i1(size))
        call field%data(field_i1)
        call load_array(name, start, end, size, nlon, buffer_i1)
        call expand(buffer_i1, field_i1, size, nproma, ngptot, nblocks)
        deallocate(buffer_i1)
      endif
    else if (frank == 3) then
      call get_offsets(start, end, size, nlon, 1, nlev, ngptot, ngptotg)
      if (rfield) then
        call field%data(field_r2)
        allocate(buffer_r2(size, nlev))
        call load_array(name, start, end, size, nlon, nlev, buffer_r2)
        call expand(buffer_r2, field_r2, size, nproma, nlev, ngptot, nblocks)
        deallocate(buffer_r2)
      endif
    else if (frank == 4) then
      ndim = field%shape(3)
      call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
      if (rfield) then
        call field%data(field_r3)
        allocate(buffer_r3(size, nlev, ndim))
        call load_array(name, start, end, size, nlon, nlev, ndim, buffer_r3)
        call expand(buffer_r3, field_r3, size, nproma, nlev, ndim, ngptot, nblocks)
        deallocate(buffer_r3)
      endif
    endif
    call field%final()
    call trace%final()
  end subroutine loadvar_atlas

  subroutine loadstate_atlas(fset, name, nlon, ngptotg)
    ! Load into the local memory buffer and expand to global field
    type(atlas_fieldset), intent(inout) :: fset
    character(len=*) :: name
    integer(kind=jpim), intent(in) :: nlon
    integer(kind=jpim), intent(in), optional :: ngptotg

    integer(kind=jpim) :: start, end, size, nlev, nproma, ngptot, nblocks, ndim
    type(atlas_field) :: field
    type(atlas_functionspace_blockstructuredcolumns) :: fspace
    type(atlas_trace) :: trace

    real(kind=jprb), allocatable :: buffer(:,:,:)
    real(kind=jprb), pointer :: field_r3(:,:,:,:)

    trace = atlas_trace("expand_atlas_mod.F90", __LINE__, "loadstate_atlas", "IO")

    field = fset%field(name)
    fspace = field%functionspace()
    nlev = field%levels()
    ngptot = fspace%size()
    !nproma = fspace%nproma()
    nproma = field%shape(1)
    nblocks = fspace%nblks()
    ndim = field%shape(3) - 3

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    allocate(buffer(size, nlev, 3+ndim))
    call field%data(field_r3)

    call load_array(name//'_T', start, end, size, nlon, nlev, buffer(:,:,1))
    call load_array(name//'_A', start, end, size, nlon, nlev, buffer(:,:,2))
    call load_array(name//'_Q', start, end, size, nlon, nlev, buffer(:,:,3))
    call load_array(name//'_CLD', start, end, size, nlon, nlev, ndim, buffer(:,:,4:))

    call expand(buffer(:,:,1), field_r3(:,:,1,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,2), field_r3(:,:,2,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field_r3(:,:,3,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,4:), field_r3(:,:,4:,:), size, nproma, nlev, ndim, ngptot, nblocks)

    deallocate(buffer)
    call field%final()
    call fspace%final()
    call trace%final()
  end subroutine loadstate_atlas

end module expand_atlas_mod
