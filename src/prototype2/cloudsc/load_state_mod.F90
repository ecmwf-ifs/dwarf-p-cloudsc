module load_state_mod
  USE PARKIND1 , ONLY : JPIM, JPRB
  USE YOMPHYDER, ONLY : STATE_TYPE

  ! Argh....
  USE YOECLDP  , ONLY : YRECLDP, TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
  USE YOMMP0   , ONLY : LSCMEC
  USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV
  USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
   & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
   & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
  USE YOEPHLI  , ONLY : YREPHLI, TEPHLI


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

  subroutine query_dimensions(KLON, KLEV, KFLDX, NAME)
    ! Initial query routine to determine data dimensions
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLON             ! Number of grid points
    INTEGER(KIND=JPIM),INTENT(OUT) :: KLEV             ! Number of levels
    ! INTEGER(KIND=JPIM),INTENT(OUT) :: NCLV
    INTEGER(KIND=JPIM),INTENT(OUT) :: KFLDX
    CHARACTER(*), INTENT(IN) :: NAME

    ! Get dimensions information from stored metadata
    call ppser_initialize(directory='data', prefix='dummy', prefix_ref=NAME)
    call fs_create_savepoint(NAME, ppser_savepoint)
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

  subroutine load_and_expand_state(name, state, field, nlon, nlev, ndim, nproma, ngptot, nblocks)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    type(state_type), pointer, intent(inout) :: state(:)
    real(kind=JPRB), pointer, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    real(kind=JPRB) :: buffer(nlon, nlev, 6+ndim)

    integer :: b

    allocate(state(nblocks))
    allocate(field(nproma, nlev, 6+ndim, nblocks))

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
    integer :: b, gidx, bidx, bsize, bend, fsize

    !omp parallel do default(shared) private(b, gidx, bidx, bsize, bend, fsize)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Index of the block in global iteration domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Field block size
       bidx = mod(gidx, nlon)  ! Rolling index into input buffer
       bend = min(bidx+bsize-1, nlon)
       if (bend-bidx+1 < bsize) then
          ! The input buffer does not hold enough data to fill field block;
          ! we need to fill the rest of the block with data from front of buffer.
          fsize = nlon - bidx + 1
          field(1:fsize,b) = buffer(bidx:nlon)
          field(fsize+1:bsize,b) = buffer(1:bsize-fsize)
       else
          ! Simply copy a block of data from the rolling buffer index
          field(1:bsize,b) = buffer(bidx:bend)
       end if

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = .FALSE.
    end do
  end subroutine expand_l1

  subroutine expand_i1(buffer, field, nlon, nproma, ngptot, nblocks)
    integer(kind=jpim), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bidx, bsize, bend, fsize

    !omp parallel do default(shared) private(b, gidx, bidx, bsize, bend, fsize)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Index of the block in global iteration domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Field block size
       bidx = mod(gidx, nlon)  ! Rolling index into input buffer
       bend = min(bidx+bsize-1, nlon)
       if (bend-bidx+1 < bsize) then
          ! The input buffer does not hold enough data to fill field block;
          ! we need to fill the rest of the block with data from front of buffer.
          fsize = nlon - bidx + 1
          field(1:fsize,b) = buffer(bidx:nlon)
          field(fsize+1:bsize,b) = buffer(1:bsize-fsize)
       else
          ! Simply copy a block of data from the rolling buffer index
          field(1:bsize,b) = buffer(bidx:bend)
       end if

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0_JPIM
    end do
  end subroutine expand_i1

  subroutine expand_r1(buffer, field, nlon, nproma, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bidx, bsize, bend, fsize

    !omp parallel do default(shared) private(b, gidx, bidx, bsize, bend, fsize)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Index of the block in global iteration domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Field block size
       bidx = mod(gidx, nlon)  ! Rolling index into input buffer
       bend = min(bidx+bsize-1, nlon)
       if (bend-bidx+1 < bsize) then
          ! The input buffer does not hold enough data to fill field block;
          ! we need to fill the rest of the block with data from front of buffer.
          fsize = nlon - bidx + 1
          field(1:fsize,b) = buffer(bidx:nlon)
          field(fsize+1:bsize,b) = buffer(1:bsize-fsize)
       else
          ! Simply copy a block of data from the rolling buffer index
          field(1:bsize,b) = buffer(bidx:bend)
       end if

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0.0_JPRB
    end do
  end subroutine expand_r1

  subroutine expand_r2(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon, nlev), field(nproma, nlev, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer :: b, gidx, bidx, bsize, bend, fsize

    !omp parallel do default(shared) private(b, gidx, bidx, bsize, bend, fsize)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block
       bidx = mod(gidx, nlon)  ! Rolling index into the input buffer
       bend = min(bidx+bsize-1, nlon)  ! Idealised final index in the input buffer
       if (bend-bidx+1 < bsize) then
          ! The input buffer does not hold enough data to fill field block;
          ! we need to fill the rest of the block with data from front of buffer.
          fsize = nlon - bidx + 1
          field(1:fsize,:,b) = buffer(bidx:nlon,:)
          field(fsize+1:bsize,:,b) = buffer(1:bsize-fsize,:)
       else
          ! Simply copy a block of data from the rolling buffer index
          field(1:bsize,:,b) = buffer(bidx:bend,:)
       end if

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r2

  subroutine expand_r3(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    real(kind=JPRB), intent(inout) :: buffer(nlon, nlev, ndim), field(nproma, nlev, ndim, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer :: b, gidx, bidx, bsize, bend, fsize

    !omp parallel do default(shared) private(b, gidx, bidx, bsize, bend, fsize)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Index of the block in global iteration domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Field block size
       bidx = mod(gidx, nlon)  ! Rolling index into input buffer
       bend = min(bidx+bsize-1, nlon)
       if (bend-bidx+1 < bsize) then
          ! The input buffer does not hold enough data to fill field block;
          ! we need to fill the rest of the block with data from front of buffer.
          fsize = nlon - bidx + 1
          field(1:fsize,:,:,b) = buffer(bidx:nlon,:,:)
          field(fsize+1:bsize,:,:,b) = buffer(1:bsize-fsize,:,:)
       else
          ! Simply copy a block of data from the rolling buffer index
          field(1:bsize,:,:,b) = buffer(bidx:bend,:,:)
       end if

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,:,b) = 0.0_JPRB
    end do
  end subroutine expand_r3


  subroutine initialise_parameters(PTSPHY)
    ! Retrieve parametersand timestep size from the serializer
    ! Note that we use `ppser_serializer_ref` to get the previously stored data
    real(kind=JPRB), intent(inout) :: PTSPHY
    
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'PTSPHY', PTSPHY)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'LSCMEC', LSCMEC)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RG', RG)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RD', RD)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RCPD', RCPD)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RETV', RETV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RLVTT', RLVTT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RLSTT', RLSTT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RLMLT', RLMLT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTT', RTT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RV', RV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R2ES', R2ES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R3LES', R3LES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R3IES', R3IES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R4LES', R4LES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R4IES', R4IES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R5LES', R5LES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R5IES', R5IES)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R5ALVCP', R5ALVCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'R5ALSCP', R5ALSCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RALVDCP', RALVDCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RALSDCP', RALSDCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RALFDCP', RALFDCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTWAT', RTWAT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTICE', RTICE)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTICECU', RTICECU)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTWAT_RTICE_R', RTWAT_RTICE_R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RTWAT_RTICECU_R', RTWAT_RTICECU_R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RKOOP1', RKOOP1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'RKOOP2', RKOOP2)

    ! Retrieve parameters contained in TECLDP type
    allocate(YRECLDP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RAMID', YRECLDP%RAMID)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLDIFF', YRECLDP%RCLDIFF)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLDIFF_CONVI', YRECLDP%RCLDIFF_CONVI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLCRIT', YRECLDP%RCLCRIT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLCRIT_SEA', YRECLDP%RCLCRIT_SEA)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLCRIT_LAND', YRECLDP%RCLCRIT_LAND)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RKCONV', YRECLDP%RKCONV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RPRC1', YRECLDP%RPRC1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RPRC2', YRECLDP%RPRC2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLDMAX', YRECLDP%RCLDMAX)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RPECONS', YRECLDP%RPECONS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RVRFACTOR', YRECLDP%RVRFACTOR)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RPRECRHMAX', YRECLDP%RPRECRHMAX)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RTAUMEL', YRECLDP%RTAUMEL)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RAMIN', YRECLDP%RAMIN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RLMIN', YRECLDP%RLMIN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RKOOPTAU', YRECLDP%RKOOPTAU)

    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLDTOPP', YRECLDP%RCLDTOPP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RLCRITSNOW', YRECLDP%RLCRITSNOW)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RSNOWLIN1', YRECLDP%RSNOWLIN1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RSNOWLIN2', YRECLDP%RSNOWLIN2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RICEHI1', YRECLDP%RICEHI1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RICEHI2', YRECLDP%RICEHI2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RICEINIT', YRECLDP%RICEINIT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RVICE', YRECLDP%RVICE)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RVRAIN', YRECLDP%RVRAIN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RVSNOW', YRECLDP%RVSNOW)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RTHOMO', YRECLDP%RTHOMO)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCOVPMIN', YRECLDP%RCOVPMIN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCCN', YRECLDP%RCCN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RNICE', YRECLDP%RNICE)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCCNOM', YRECLDP%RCCNOM)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCCNSS', YRECLDP%RCCNSS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCCNSU', YRECLDP%RCCNSU)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCLDTOPCF', YRECLDP%RCLDTOPCF)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RDEPLIQREFRATE', YRECLDP%RDEPLIQREFRATE)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RDEPLIQREFDEPTH', YRECLDP%RDEPLIQREFDEPTH)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KKAac', YRECLDP%RCL_KKAac)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KKBac', YRECLDP%RCL_KKBac)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KKAau', YRECLDP%RCL_KKAau)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KKBauq', YRECLDP%RCL_KKBauq)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KKBaun', YRECLDP%RCL_KKBaun)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KK_cloud_num_sea', YRECLDP%RCL_KK_cloud_num_sea)
    Call Fs_get_serializer_metainfo(Ppser_serializer_ref, 'YRECLDP_RCL_KK_cloud_num_land', YRECLDP%RCL_KK_cloud_num_land)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_AI', YRECLDP%RCL_AI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_BI', YRECLDP%RCL_BI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CI', YRECLDP%RCL_CI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_DI', YRECLDP%RCL_DI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X1I', YRECLDP%RCL_X1I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X2I', YRECLDP%RCL_X2I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X3I', YRECLDP%RCL_X3I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X4I', YRECLDP%RCL_X4I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST1I', YRECLDP%RCL_CONST1I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST2I', YRECLDP%RCL_CONST2I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST3I', YRECLDP%RCL_CONST3I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST4I', YRECLDP%RCL_CONST4I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST5I', YRECLDP%RCL_CONST5I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST6I', YRECLDP%RCL_CONST6I)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_APB1', YRECLDP%RCL_APB1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_APB2', YRECLDP%RCL_APB2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_APB3', YRECLDP%RCL_APB3)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_AS', YRECLDP%RCL_AS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_BS', YRECLDP%RCL_BS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CS', YRECLDP%RCL_CS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_DS', YRECLDP%RCL_DS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X1S', YRECLDP%RCL_X1S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X2S', YRECLDP%RCL_X2S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X3S', YRECLDP%RCL_X3S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X4S', YRECLDP%RCL_X4S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST1S', YRECLDP%RCL_CONST1S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST2S', YRECLDP%RCL_CONST2S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST3S', YRECLDP%RCL_CONST3S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST4S', YRECLDP%RCL_CONST4S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST5S', YRECLDP%RCL_CONST5S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST6S', YRECLDP%RCL_CONST6S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST7S', YRECLDP%RCL_CONST7S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST8S', YRECLDP%RCL_CONST8S)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RDENSWAT', YRECLDP%RDENSWAT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RDENSREF', YRECLDP%RDENSREF)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_AR', YRECLDP%RCL_AR)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_BR', YRECLDP%RCL_BR)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CR', YRECLDP%RCL_CR)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_DR', YRECLDP%RCL_DR)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X1R', YRECLDP%RCL_X1R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X2R', YRECLDP%RCL_X2R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_X4R', YRECLDP%RCL_X4R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_KA273', YRECLDP%RCL_KA273)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CDENOM1', YRECLDP%RCL_CDENOM1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CDENOM2', YRECLDP%RCL_CDENOM2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CDENOM3', YRECLDP%RCL_CDENOM3)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_SCHMIDT', YRECLDP%RCL_SCHMIDT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_DYNVISC', YRECLDP%RCL_DYNVISC)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST1R', YRECLDP%RCL_CONST1R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST2R', YRECLDP%RCL_CONST2R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST3R', YRECLDP%RCL_CONST3R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST4R', YRECLDP%RCL_CONST4R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_FAC1', YRECLDP%RCL_FAC1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_FAC2', YRECLDP%RCL_FAC2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST5R', YRECLDP%RCL_CONST5R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_CONST6R', YRECLDP%RCL_CONST6R)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_FZRAB', YRECLDP%RCL_FZRAB)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_RCL_FZRBB', YRECLDP%RCL_FZRBB)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LCLDEXTRA', YRECLDP%LCLDEXTRA)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LCLDBUDGET', YRECLDP%LCLDBUDGET)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NSSOPT', YRECLDP%NSSOPT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NCLDTOP', YRECLDP%NCLDTOP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAECLBC', YRECLDP%NAECLBC)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAECLDU', YRECLDP%NAECLDU)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAECLOM', YRECLDP%NAECLOM)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAECLSS', YRECLDP%NAECLSS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAECLSU', YRECLDP%NAECLSU)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NCLDDIAG', YRECLDP%NCLDDIAG)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NAERCLD', YRECLDP%NAERCLD)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERLIQAUTOLSP', YRECLDP%LAERLIQAUTOLSP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERLIQAUTOCP', YRECLDP%LAERLIQAUTOCP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERLIQAUTOCPB', YRECLDP%LAERLIQAUTOCPB)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERLIQCOLL', YRECLDP%LAERLIQCOLL)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERICESED', YRECLDP%LAERICESED)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_LAERICEAUTO', YRECLDP%LAERICEAUTO)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NSHAPEP', YRECLDP%NSHAPEP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NSHAPEQ', YRECLDP%NSHAPEQ)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YRECLDP_NBETA', YRECLDP%NBETA)
    ! The last two are actually arrays, so treat them as fields
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'YRECLDP_RBETA', YRECLDP%RBETA(0:100))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'YRECLDP_RBETAP1', YRECLDP%RBETAP1(0:100))

    ! Retrieve parameters contained in TECLDP type
    allocate(YREPHLI)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_LTLEVOL', YREPHLI%LTLEVOL)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_LPHYLIN', YREPHLI%LPHYLIN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_LENOPERT', YREPHLI%LENOPERT)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_LEPPCFLS', YREPHLI%LEPPCFLS)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_LRAISANEN', YREPHLI%LRAISANEN)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPTRC', YREPHLI%RLPTRC)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPAL1', YREPHLI%RLPAL1)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPAL2', YREPHLI%RLPAL2)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPBB', YREPHLI%RLPBB)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPCC', YREPHLI%RLPCC)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPDD', YREPHLI%RLPDD)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPMIXL', YREPHLI%RLPMIXL)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPBETA', YREPHLI%RLPBETA)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPDRAG', YREPHLI%RLPDRAG)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPEVAP', YREPHLI%RLPEVAP)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'YREPHLI_RLPP00', YREPHLI%RLPP00)
  end subroutine initialise_parameters

  subroutine finalize()
    call ppser_finalize()
  end subroutine finalize
end module load_state_mod
