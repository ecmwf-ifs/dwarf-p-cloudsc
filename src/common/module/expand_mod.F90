module expand_mod
  USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
  USE YOMPHYDER, ONLY : STATE_TYPE

  USE YOECLDP  , ONLY : YRECLDP, TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV
  USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
   & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
   & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
  USE YOEPHLI  , ONLY : YREPHLI, TEPHLI
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
    integer(kind=jpim), allocatable :: buffer(:), rbuf(:)
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
    real(kind=jprd), allocatable :: buffer(:), rbuf(:)
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
    real(kind=jprd), allocatable :: buffer(:,:), rbuf(:,:)
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
    real(kind=jprd), allocatable :: buffer(:,:,:), rbuf(:,:,:)
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
    real(kind=jprd), allocatable :: buffer(:,:,:), rbuf(:,:,:)
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
    real(kind=JPRD), intent(inout) :: buffer(nlon)
    real(kind=JPRB), intent(inout) :: field(nproma, nblocks)
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
    real(kind=JPRD), intent(inout) :: buffer(nlon, nlev)
    real(kind=JPRB), intent(inout) :: field(nproma, nlev, nblocks)
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
    real(kind=JPRD), intent(inout) :: buffer(nlon, nlev, ndim)
    real(kind=JPRB), intent(inout) :: field(nproma, nlev, ndim, nblocks)
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


  subroutine initialise_parameters(PTSPHY, LDSLPHY, LDMAINCALL)
    ! Retrieve parametersand timestep size from the serializer
    ! Note that we use `ppser_serializer_ref` to get the previously stored data
    real(kind=JPRB), intent(inout) :: PTSPHY
    logical, intent(inout) :: LDSLPHY, LDMAINCALL

    call load_scalar('PTSPHY', PTSPHY)
    call load_scalar('LDSLPHY', LDSLPHY)
    call load_scalar('LDMAINCALL', LDMAINCALL)

    call load_scalar('RG', RG)
    call load_scalar('RD', RD)
    call load_scalar('RCPD', RCPD)
    call load_scalar('RETV', RETV)
    call load_scalar('RLVTT', RLVTT)
    call load_scalar('RLSTT', RLSTT)
    call load_scalar('RLMLT', RLMLT)
    call load_scalar('RTT', RTT)
    call load_scalar('RV', RV)
    call load_scalar('R2ES', R2ES)
    call load_scalar('R3LES', R3LES)
    call load_scalar('R3IES', R3IES)
    call load_scalar('R4LES', R4LES)
    call load_scalar('R4IES', R4IES)
    call load_scalar('R5LES', R5LES)
    call load_scalar('R5IES', R5IES)
    call load_scalar('R5ALVCP', R5ALVCP)
    call load_scalar('R5ALSCP', R5ALSCP)
    call load_scalar('RALVDCP', RALVDCP)
    call load_scalar('RALSDCP', RALSDCP)
    call load_scalar('RALFDCP', RALFDCP)
    call load_scalar('RTWAT', RTWAT)
    call load_scalar('RTICE', RTICE)
    call load_scalar('RTICECU', RTICECU)
    call load_scalar('RTWAT_RTICE_R', RTWAT_RTICE_R)
    call load_scalar('RTWAT_RTICECU_R', RTWAT_RTICECU_R)
    call load_scalar('RKOOP1', RKOOP1)
    call load_scalar('RKOOP2', RKOOP2)

    ! Retrieve parameters contained in TECLDP type
    allocate(YRECLDP)
    call load_scalar('YRECLDP_RAMID', YRECLDP%RAMID)
    call load_scalar('YRECLDP_RCLDIFF', YRECLDP%RCLDIFF)
    call load_scalar('YRECLDP_RCLDIFF_CONVI', YRECLDP%RCLDIFF_CONVI)
    call load_scalar('YRECLDP_RCLCRIT', YRECLDP%RCLCRIT)
    call load_scalar('YRECLDP_RCLCRIT_SEA', YRECLDP%RCLCRIT_SEA)
    call load_scalar('YRECLDP_RCLCRIT_LAND', YRECLDP%RCLCRIT_LAND)
    call load_scalar('YRECLDP_RKCONV', YRECLDP%RKCONV)
    call load_scalar('YRECLDP_RPRC1', YRECLDP%RPRC1)
    call load_scalar('YRECLDP_RPRC2', YRECLDP%RPRC2)
    call load_scalar('YRECLDP_RCLDMAX', YRECLDP%RCLDMAX)
    call load_scalar('YRECLDP_RPECONS', YRECLDP%RPECONS)
    call load_scalar('YRECLDP_RVRFACTOR', YRECLDP%RVRFACTOR)
    call load_scalar('YRECLDP_RPRECRHMAX', YRECLDP%RPRECRHMAX)
    call load_scalar('YRECLDP_RTAUMEL', YRECLDP%RTAUMEL)
    call load_scalar('YRECLDP_RAMIN', YRECLDP%RAMIN)
    call load_scalar('YRECLDP_RLMIN', YRECLDP%RLMIN)
    call load_scalar('YRECLDP_RKOOPTAU', YRECLDP%RKOOPTAU)

    call load_scalar('YRECLDP_RCLDTOPP', YRECLDP%RCLDTOPP)
    call load_scalar('YRECLDP_RLCRITSNOW', YRECLDP%RLCRITSNOW)
    call load_scalar('YRECLDP_RSNOWLIN1', YRECLDP%RSNOWLIN1)
    call load_scalar('YRECLDP_RSNOWLIN2', YRECLDP%RSNOWLIN2)
    call load_scalar('YRECLDP_RICEHI1', YRECLDP%RICEHI1)
    call load_scalar('YRECLDP_RICEHI2', YRECLDP%RICEHI2)
    call load_scalar('YRECLDP_RICEINIT', YRECLDP%RICEINIT)
    call load_scalar('YRECLDP_RVICE', YRECLDP%RVICE)
    call load_scalar('YRECLDP_RVRAIN', YRECLDP%RVRAIN)
    call load_scalar('YRECLDP_RVSNOW', YRECLDP%RVSNOW)
    call load_scalar('YRECLDP_RTHOMO', YRECLDP%RTHOMO)
    call load_scalar('YRECLDP_RCOVPMIN', YRECLDP%RCOVPMIN)
    call load_scalar('YRECLDP_RCCN', YRECLDP%RCCN)
    call load_scalar('YRECLDP_RNICE', YRECLDP%RNICE)
    call load_scalar('YRECLDP_RCCNOM', YRECLDP%RCCNOM)
    call load_scalar('YRECLDP_RCCNSS', YRECLDP%RCCNSS)
    call load_scalar('YRECLDP_RCCNSU', YRECLDP%RCCNSU)
    call load_scalar('YRECLDP_RCLDTOPCF', YRECLDP%RCLDTOPCF)
    call load_scalar('YRECLDP_RDEPLIQREFRATE', YRECLDP%RDEPLIQREFRATE)
    call load_scalar('YRECLDP_RDEPLIQREFDEPTH', YRECLDP%RDEPLIQREFDEPTH)
    Call load_scalar('YRECLDP_RCL_KKAac', YRECLDP%RCL_KKAac)
    Call load_scalar('YRECLDP_RCL_KKBac', YRECLDP%RCL_KKBac)
    Call load_scalar('YRECLDP_RCL_KKAau', YRECLDP%RCL_KKAau)
    Call load_scalar('YRECLDP_RCL_KKBauq', YRECLDP%RCL_KKBauq)
    Call load_scalar('YRECLDP_RCL_KKBaun', YRECLDP%RCL_KKBaun)
    Call load_scalar('YRECLDP_RCL_KK_cloud_num_sea', YRECLDP%RCL_KK_cloud_num_sea)
    Call load_scalar('YRECLDP_RCL_KK_cloud_num_land', YRECLDP%RCL_KK_cloud_num_land)
    call load_scalar('YRECLDP_RCL_AI', YRECLDP%RCL_AI)
    call load_scalar('YRECLDP_RCL_BI', YRECLDP%RCL_BI)
    call load_scalar('YRECLDP_RCL_CI', YRECLDP%RCL_CI)
    call load_scalar('YRECLDP_RCL_DI', YRECLDP%RCL_DI)
    call load_scalar('YRECLDP_RCL_X1I', YRECLDP%RCL_X1I)
    call load_scalar('YRECLDP_RCL_X2I', YRECLDP%RCL_X2I)
    call load_scalar('YRECLDP_RCL_X3I', YRECLDP%RCL_X3I)
    call load_scalar('YRECLDP_RCL_X4I', YRECLDP%RCL_X4I)
    call load_scalar('YRECLDP_RCL_CONST1I', YRECLDP%RCL_CONST1I)
    call load_scalar('YRECLDP_RCL_CONST2I', YRECLDP%RCL_CONST2I)
    call load_scalar('YRECLDP_RCL_CONST3I', YRECLDP%RCL_CONST3I)
    call load_scalar('YRECLDP_RCL_CONST4I', YRECLDP%RCL_CONST4I)
    call load_scalar('YRECLDP_RCL_CONST5I', YRECLDP%RCL_CONST5I)
    call load_scalar('YRECLDP_RCL_CONST6I', YRECLDP%RCL_CONST6I)
    call load_scalar('YRECLDP_RCL_APB1', YRECLDP%RCL_APB1)
    call load_scalar('YRECLDP_RCL_APB2', YRECLDP%RCL_APB2)
    call load_scalar('YRECLDP_RCL_APB3', YRECLDP%RCL_APB3)
    call load_scalar('YRECLDP_RCL_AS', YRECLDP%RCL_AS)
    call load_scalar('YRECLDP_RCL_BS', YRECLDP%RCL_BS)
    call load_scalar('YRECLDP_RCL_CS', YRECLDP%RCL_CS)
    call load_scalar('YRECLDP_RCL_DS', YRECLDP%RCL_DS)
    call load_scalar('YRECLDP_RCL_X1S', YRECLDP%RCL_X1S)
    call load_scalar('YRECLDP_RCL_X2S', YRECLDP%RCL_X2S)
    call load_scalar('YRECLDP_RCL_X3S', YRECLDP%RCL_X3S)
    call load_scalar('YRECLDP_RCL_X4S', YRECLDP%RCL_X4S)
    call load_scalar('YRECLDP_RCL_CONST1S', YRECLDP%RCL_CONST1S)
    call load_scalar('YRECLDP_RCL_CONST2S', YRECLDP%RCL_CONST2S)
    call load_scalar('YRECLDP_RCL_CONST3S', YRECLDP%RCL_CONST3S)
    call load_scalar('YRECLDP_RCL_CONST4S', YRECLDP%RCL_CONST4S)
    call load_scalar('YRECLDP_RCL_CONST5S', YRECLDP%RCL_CONST5S)
    call load_scalar('YRECLDP_RCL_CONST6S', YRECLDP%RCL_CONST6S)
    call load_scalar('YRECLDP_RCL_CONST7S', YRECLDP%RCL_CONST7S)
    call load_scalar('YRECLDP_RCL_CONST8S', YRECLDP%RCL_CONST8S)
    call load_scalar('YRECLDP_RDENSWAT', YRECLDP%RDENSWAT)
    call load_scalar('YRECLDP_RDENSREF', YRECLDP%RDENSREF)
    call load_scalar('YRECLDP_RCL_AR', YRECLDP%RCL_AR)
    call load_scalar('YRECLDP_RCL_BR', YRECLDP%RCL_BR)
    call load_scalar('YRECLDP_RCL_CR', YRECLDP%RCL_CR)
    call load_scalar('YRECLDP_RCL_DR', YRECLDP%RCL_DR)
    call load_scalar('YRECLDP_RCL_X1R', YRECLDP%RCL_X1R)
    call load_scalar('YRECLDP_RCL_X2R', YRECLDP%RCL_X2R)
    call load_scalar('YRECLDP_RCL_X4R', YRECLDP%RCL_X4R)
    call load_scalar('YRECLDP_RCL_KA273', YRECLDP%RCL_KA273)
    call load_scalar('YRECLDP_RCL_CDENOM1', YRECLDP%RCL_CDENOM1)
    call load_scalar('YRECLDP_RCL_CDENOM2', YRECLDP%RCL_CDENOM2)
    call load_scalar('YRECLDP_RCL_CDENOM3', YRECLDP%RCL_CDENOM3)
    call load_scalar('YRECLDP_RCL_SCHMIDT', YRECLDP%RCL_SCHMIDT)
    call load_scalar('YRECLDP_RCL_DYNVISC', YRECLDP%RCL_DYNVISC)
    call load_scalar('YRECLDP_RCL_CONST1R', YRECLDP%RCL_CONST1R)
    call load_scalar('YRECLDP_RCL_CONST2R', YRECLDP%RCL_CONST2R)
    call load_scalar('YRECLDP_RCL_CONST3R', YRECLDP%RCL_CONST3R)
    call load_scalar('YRECLDP_RCL_CONST4R', YRECLDP%RCL_CONST4R)
    call load_scalar('YRECLDP_RCL_FAC1', YRECLDP%RCL_FAC1)
    call load_scalar('YRECLDP_RCL_FAC2', YRECLDP%RCL_FAC2)
    call load_scalar('YRECLDP_RCL_CONST5R', YRECLDP%RCL_CONST5R)
    call load_scalar('YRECLDP_RCL_CONST6R', YRECLDP%RCL_CONST6R)
    call load_scalar('YRECLDP_RCL_FZRAB', YRECLDP%RCL_FZRAB)
    call load_scalar('YRECLDP_RCL_FZRBB', YRECLDP%RCL_FZRBB)
    call load_scalar('YRECLDP_LCLDEXTRA', YRECLDP%LCLDEXTRA)
    call load_scalar('YRECLDP_LCLDBUDGET', YRECLDP%LCLDBUDGET)
    call load_scalar('YRECLDP_NSSOPT', YRECLDP%NSSOPT)
    call load_scalar('YRECLDP_NCLDTOP', YRECLDP%NCLDTOP)
    call load_scalar('YRECLDP_NAECLBC', YRECLDP%NAECLBC)
    call load_scalar('YRECLDP_NAECLDU', YRECLDP%NAECLDU)
    call load_scalar('YRECLDP_NAECLOM', YRECLDP%NAECLOM)
    call load_scalar('YRECLDP_NAECLSS', YRECLDP%NAECLSS)
    call load_scalar('YRECLDP_NAECLSU', YRECLDP%NAECLSU)
    call load_scalar('YRECLDP_NCLDDIAG', YRECLDP%NCLDDIAG)
    call load_scalar('YRECLDP_NAERCLD', YRECLDP%NAERCLD)
    call load_scalar('YRECLDP_LAERLIQAUTOLSP', YRECLDP%LAERLIQAUTOLSP)
    call load_scalar('YRECLDP_LAERLIQAUTOCP', YRECLDP%LAERLIQAUTOCP)
    call load_scalar('YRECLDP_LAERLIQAUTOCPB', YRECLDP%LAERLIQAUTOCPB)
    call load_scalar('YRECLDP_LAERLIQCOLL', YRECLDP%LAERLIQCOLL)
    call load_scalar('YRECLDP_LAERICESED', YRECLDP%LAERICESED)
    call load_scalar('YRECLDP_LAERICEAUTO', YRECLDP%LAERICEAUTO)
    call load_scalar('YRECLDP_NSHAPEP', YRECLDP%NSHAPEP)
    call load_scalar('YRECLDP_NSHAPEQ', YRECLDP%NSHAPEQ)
    call load_scalar('YRECLDP_NBETA', YRECLDP%NBETA)
    ! The last two are actually arrays, so treat them as fields
    call load_array('YRECLDP_RBETA', 1, 101, 101, 101, YRECLDP%RBETA(0:100))
    call load_array('YRECLDP_RBETAP1', 1, 101, 101, 101, YRECLDP%RBETAP1(0:100))

    ! Retrieve parameters contained in TECLDP type
    allocate(YREPHLI)
    call load_scalar('YREPHLI_LTLEVOL', YREPHLI%LTLEVOL)
    call load_scalar('YREPHLI_LPHYLIN', YREPHLI%LPHYLIN)
    call load_scalar('YREPHLI_LENOPERT', YREPHLI%LENOPERT)
    call load_scalar('YREPHLI_LEPPCFLS', YREPHLI%LEPPCFLS)
    call load_scalar('YREPHLI_LRAISANEN', YREPHLI%LRAISANEN)
    call load_scalar('YREPHLI_RLPTRC', YREPHLI%RLPTRC)
    call load_scalar('YREPHLI_RLPAL1', YREPHLI%RLPAL1)
    call load_scalar('YREPHLI_RLPAL2', YREPHLI%RLPAL2)
    call load_scalar('YREPHLI_RLPBB', YREPHLI%RLPBB)
    call load_scalar('YREPHLI_RLPCC', YREPHLI%RLPCC)
    call load_scalar('YREPHLI_RLPDD', YREPHLI%RLPDD)
    call load_scalar('YREPHLI_RLPMIXL', YREPHLI%RLPMIXL)
    call load_scalar('YREPHLI_RLPBETA', YREPHLI%RLPBETA)
    call load_scalar('YREPHLI_RLPDRAG', YREPHLI%RLPDRAG)
    call load_scalar('YREPHLI_RLPEVAP', YREPHLI%RLPEVAP)
    call load_scalar('YREPHLI_RLPP00', YREPHLI%RLPP00)
  end subroutine initialise_parameters

end module expand_mod
