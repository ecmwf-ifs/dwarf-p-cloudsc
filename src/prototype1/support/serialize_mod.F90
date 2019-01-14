module serialize_mod
  USE PARKIND1 , ONLY : JPIM, JPRB
  USE YOMPHYDER, ONLY : STATE_TYPE
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
    call fs_create_savepoint('input', ppser_savepoint)
    call ppser_set_mode(1)

    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLON', KLON)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KLEV', KLEV)
    ! call fs_get_serializer_metainfo(ppser_serializer_ref, 'NCLV', NCLV)
    call fs_get_serializer_metainfo(ppser_serializer_ref, 'KFLDX', KFLDX)

  end subroutine query_dimensions

  subroutine serialize( &
   & KLON, KLEV, PTSPHY,&
   & PT, PQ, TENDENCY_CML, TENDENCY_TMP, &
   & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
   & PHRSW, PHRLW, PVERVEL, PAP, PAPH, &
   & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, &
   & LDSLPHY, LDMAINCALL, PA, PCLV, PSUPSAT, &
   & PLCRIT_AER, PICRIT_AER, PRE_ICE, &
   & PCCN, PNICE, PEXTRA, KFLDX &
   & )
   ! Serialization routine used to generate the input data

    INTEGER(KIND=JPIM),INTENT(IN)       :: KLON             ! Number of grid points
    INTEGER(KIND=JPIM),INTENT(IN)       :: KLEV             ! Number of levels
    INTEGER(KIND=JPIM),INTENT(IN)       :: KFLDX

    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PTSPHY           ! Physics timestep
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PT(KLON,KLEV)    ! T at start of callpar
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PQ(KLON,KLEV)    ! Q at start of callpar
    TYPE (STATE_TYPE) ,INTENT(INOUT)    :: TENDENCY_CML     ! cumulative tendency used for final output
    TYPE (STATE_TYPE) ,INTENT(INOUT)    :: TENDENCY_TMP     ! cumulative tendency used as input
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PDYNL(KLON,KLEV) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PDYNI(KLON,KLEV) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVERVEL(KLON,KLEV) !Vertical velocity
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PAP(KLON,KLEV)   ! Pressure on full levels
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PLSM(KLON)       ! Land fraction (0-1) 
    LOGICAL           ,INTENT(INOUT)    :: LDCUM(KLON)      ! Convection active
    INTEGER(KIND=JPIM),INTENT(INOUT)    :: KTYPE(KLON)      ! Convection type 0,1,2
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PLU(KLON,KLEV)   ! Conv. condensate
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PLUDE(KLON,KLEV) ! Conv. detrained water 
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PSNDE(KLON,KLEV) ! Conv. detrained snow
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
    LOGICAL           ,INTENT(INOUT)    :: LDSLPHY 
    LOGICAL           ,INTENT(INOUT)    :: LDMAINCALL       ! T if main call to cloudsc
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCLV(KLON,KLEV,NCLV) 
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PSUPSAT(KLON,KLEV)
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PLCRIT_AER(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PICRIT_AER(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PRE_ICE(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields

    ! Initialize serializer for storing reference input
    call ppser_initialize(directory='data', prefix='input')
    call fs_create_savepoint('input', ppser_savepoint)
    call ppser_set_mode(0)

    ! Store dimensions and timestep size on the serializer
    call fs_add_serializer_metainfo(ppser_serializer, 'KLON', KLON)
    call fs_add_serializer_metainfo(ppser_serializer, 'KLEV', KLEV)
    call fs_add_serializer_metainfo(ppser_serializer, 'KFLDX', KFLDX)
    call fs_add_serializer_metainfo(ppser_serializer, 'PTSPHY', PTSPHY)

    ! Store parameters on the serializer
    call fs_add_serializer_metainfo(ppser_serializer, 'LSCMEC', LSCMEC)
    call fs_add_serializer_metainfo(ppser_serializer, 'RG', RG)
    call fs_add_serializer_metainfo(ppser_serializer, 'RD', RD)
    call fs_add_serializer_metainfo(ppser_serializer, 'RCPD', RCPD)
    call fs_add_serializer_metainfo(ppser_serializer, 'RETV', RETV)
    call fs_add_serializer_metainfo(ppser_serializer, 'RLVTT', RLVTT)
    call fs_add_serializer_metainfo(ppser_serializer, 'RLSTT', RLSTT)
    call fs_add_serializer_metainfo(ppser_serializer, 'RLMLT', RLMLT)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTT', RTT)
    call fs_add_serializer_metainfo(ppser_serializer, 'RV', RV)
    call fs_add_serializer_metainfo(ppser_serializer, 'R2ES', R2ES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R3LES', R3LES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R3IES', R3IES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R4LES', R4LES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R4IES', R4IES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R5LES', R5LES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R5IES', R5IES)
    call fs_add_serializer_metainfo(ppser_serializer, 'R5ALVCP', R5ALVCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'R5ALSCP', R5ALSCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'RALVDCP', RALVDCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'RALSDCP', RALSDCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'RALFDCP', RALFDCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTWAT', RTWAT)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTICE', RTICE)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTICECU', RTICECU)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTWAT_RTICE_R', RTWAT_RTICE_R)
    call fs_add_serializer_metainfo(ppser_serializer, 'RTWAT_RTICECU_R', RTWAT_RTICECU_R)
    call fs_add_serializer_metainfo(ppser_serializer, 'RKOOP1', RKOOP1)
    call fs_add_serializer_metainfo(ppser_serializer, 'RKOOP2', RKOOP2)

    ! Store parameters contained in TECLDP type
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RAMID', YRECLDP%RAMID)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLDIFF', YRECLDP%RCLDIFF)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLDIFF_CONVI', YRECLDP%RCLDIFF_CONVI)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLCRIT', YRECLDP%RCLCRIT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLCRIT_SEA', YRECLDP%RCLCRIT_SEA)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLCRIT_LAND', YRECLDP%RCLCRIT_LAND)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RKCONV', YRECLDP%RKCONV)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RPRC1', YRECLDP%RPRC1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RPRC2', YRECLDP%RPRC2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLDMAX', YRECLDP%RCLDMAX)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RPECONS', YRECLDP%RPECONS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RVRFACTOR', YRECLDP%RVRFACTOR)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RPRECRHMAX', YRECLDP%RPRECRHMAX)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RTAUMEL', YRECLDP%RTAUMEL)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RAMIN', YRECLDP%RAMIN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RLMIN', YRECLDP%RLMIN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RKOOPTAU', YRECLDP%RKOOPTAU)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLDTOPP', YRECLDP%RCLDTOPP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RLCRITSNOW', YRECLDP%RLCRITSNOW)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RSNOWLIN1', YRECLDP%RSNOWLIN1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RSNOWLIN2', YRECLDP%RSNOWLIN2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RICEHI1', YRECLDP%RICEHI1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RICEHI2', YRECLDP%RICEHI2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RICEINIT', YRECLDP%RICEINIT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RVICE', YRECLDP%RVICE)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RVRAIN', YRECLDP%RVRAIN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RVSNOW', YRECLDP%RVSNOW)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RTHOMO', YRECLDP%RTHOMO)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCOVPMIN', YRECLDP%RCOVPMIN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCCN', YRECLDP%RCCN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RNICE', YRECLDP%RNICE)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCCNOM', YRECLDP%RCCNOM)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCCNSS', YRECLDP%RCCNSS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCCNSU', YRECLDP%RCCNSU)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCLDTOPCF', YRECLDP%RCLDTOPCF)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RDEPLIQREFRATE', YRECLDP%RDEPLIQREFRATE)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RDEPLIQREFDEPTH', YRECLDP%RDEPLIQREFDEPTH)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KKAac', YRECLDP%RCL_KKAac)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KKBac', YRECLDP%RCL_KKBac)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KKAau', YRECLDP%RCL_KKAau)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KKBauq', YRECLDP%RCL_KKBauq)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KKBaun', YRECLDP%RCL_KKBaun)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KK_cloud_num_sea', YRECLDP%RCL_KK_cloud_num_sea)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KK_cloud_num_land', YRECLDP%RCL_KK_cloud_num_land)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_AI', YRECLDP%RCL_AI)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_BI', YRECLDP%RCL_BI)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CI', YRECLDP%RCL_CI)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_DI', YRECLDP%RCL_DI)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X1I', YRECLDP%RCL_X1I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X2I', YRECLDP%RCL_X2I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X3I', YRECLDP%RCL_X3I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X4I', YRECLDP%RCL_X4I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST1I', YRECLDP%RCL_CONST1I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST2I', YRECLDP%RCL_CONST2I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST3I', YRECLDP%RCL_CONST3I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST4I', YRECLDP%RCL_CONST4I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST5I', YRECLDP%RCL_CONST5I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST6I', YRECLDP%RCL_CONST6I)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_APB1', YRECLDP%RCL_APB1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_APB2', YRECLDP%RCL_APB2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_APB3', YRECLDP%RCL_APB3)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_AS', YRECLDP%RCL_AS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_BS', YRECLDP%RCL_BS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CS', YRECLDP%RCL_CS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_DS', YRECLDP%RCL_DS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X1S', YRECLDP%RCL_X1S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X2S', YRECLDP%RCL_X2S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X3S', YRECLDP%RCL_X3S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X4S', YRECLDP%RCL_X4S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST1S', YRECLDP%RCL_CONST1S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST2S', YRECLDP%RCL_CONST2S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST3S', YRECLDP%RCL_CONST3S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST4S', YRECLDP%RCL_CONST4S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST5S', YRECLDP%RCL_CONST5S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST6S', YRECLDP%RCL_CONST6S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST7S', YRECLDP%RCL_CONST7S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST8S', YRECLDP%RCL_CONST8S)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RDENSWAT', YRECLDP%RDENSWAT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RDENSREF', YRECLDP%RDENSREF)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_AR', YRECLDP%RCL_AR)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_BR', YRECLDP%RCL_BR)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CR', YRECLDP%RCL_CR)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_DR', YRECLDP%RCL_DR)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X1R', YRECLDP%RCL_X1R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X2R', YRECLDP%RCL_X2R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_X4R', YRECLDP%RCL_X4R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_KA273', YRECLDP%RCL_KA273)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CDENOM1', YRECLDP%RCL_CDENOM1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CDENOM2', YRECLDP%RCL_CDENOM2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CDENOM3', YRECLDP%RCL_CDENOM3)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_SCHMIDT', YRECLDP%RCL_SCHMIDT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_DYNVISC', YRECLDP%RCL_DYNVISC)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST1R', YRECLDP%RCL_CONST1R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST2R', YRECLDP%RCL_CONST2R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST3R', YRECLDP%RCL_CONST3R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST4R', YRECLDP%RCL_CONST4R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_FAC1', YRECLDP%RCL_FAC1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_FAC2', YRECLDP%RCL_FAC2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST5R', YRECLDP%RCL_CONST5R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_CONST6R', YRECLDP%RCL_CONST6R)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_FZRAB', YRECLDP%RCL_FZRAB)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_RCL_FZRBB', YRECLDP%RCL_FZRBB)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LCLDEXTRA', YRECLDP%LCLDEXTRA)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LCLDBUDGET', YRECLDP%LCLDBUDGET)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NSSOPT', YRECLDP%NSSOPT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NCLDTOP', YRECLDP%NCLDTOP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAECLBC', YRECLDP%NAECLBC)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAECLDU', YRECLDP%NAECLDU)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAECLOM', YRECLDP%NAECLOM)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAECLSS', YRECLDP%NAECLSS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAECLSU', YRECLDP%NAECLSU)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NCLDDIAG', YRECLDP%NCLDDIAG)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NAERCLD', YRECLDP%NAERCLD)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERLIQAUTOLSP', YRECLDP%LAERLIQAUTOLSP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERLIQAUTOCP', YRECLDP%LAERLIQAUTOCP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERLIQAUTOCPB', YRECLDP%LAERLIQAUTOCPB)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERLIQCOLL', YRECLDP%LAERLIQCOLL)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERICESED', YRECLDP%LAERICESED)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_LAERICEAUTO', YRECLDP%LAERICEAUTO)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NSHAPEP', YRECLDP%NSHAPEP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NSHAPEQ', YRECLDP%NSHAPEQ)
    call fs_add_serializer_metainfo(ppser_serializer, 'YRECLDP_NBETA', YRECLDP%NBETA)
    ! The last two are actually arrays, so treat them as fields
    call fs_write_field(ppser_serializer, ppser_savepoint, 'YRECLDP_RBETA', YRECLDP%RBETA(0:100))
    call fs_write_field(ppser_serializer, ppser_savepoint, 'YRECLDP_RBETAP1', YRECLDP%RBETAP1(0:100))

    ! Store parameters contained in TECLDP type
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_LTLEVOL', YREPHLI%LTLEVOL)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_LPHYLIN', YREPHLI%LPHYLIN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_LENOPERT', YREPHLI%LENOPERT)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_LEPPCFLS', YREPHLI%LEPPCFLS)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_LRAISANEN', YREPHLI%LRAISANEN)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPTRC', YREPHLI%RLPTRC)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPAL1', YREPHLI%RLPAL1)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPAL2', YREPHLI%RLPAL2)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPBB', YREPHLI%RLPBB)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPCC', YREPHLI%RLPCC)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPDD', YREPHLI%RLPDD)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPMIXL', YREPHLI%RLPMIXL)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPBETA', YREPHLI%RLPBETA)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPDRAG', YREPHLI%RLPDRAG)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPEVAP', YREPHLI%RLPEVAP)
    call fs_add_serializer_metainfo(ppser_serializer, 'YREPHLI_RLPP00', YREPHLI%RLPP00)

    ! Dump field data arrays to savepoint
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PT', PT)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PQ', PQ)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_CML_T', TENDENCY_CML%T)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_CML_Q', TENDENCY_CML%Q)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_CML_A', TENDENCY_CML%A)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_CML_CLD', TENDENCY_CML%CLD)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_TMP_T', TENDENCY_TMP%T)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_TMP_Q', TENDENCY_TMP%Q)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_TMP_A', TENDENCY_TMP%A)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'TENDENCY_TMP_CLD', TENDENCY_TMP%CLD)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PVFA', PVFA)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PVFL', PVFL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PVFI', PVFI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PDYNA', PDYNA)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PDYNL', PDYNL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PDYNI', PDYNI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PHRSW', PHRSW)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PHRLW', PHRLW)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PVERVEL', PVERVEL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PAP', PAP)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PAPH', PAPH)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLSM', PLSM)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'LDCUM', LDCUM)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'KTYPE', KTYPE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLU', PLU)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLUDE', PLUDE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PSNDE', PSNDE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PMFU', PMFU)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PMFD', PMFD)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'LDSLPHY', LDSLPHY)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'LDMAINCALL', LDMAINCALL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PA', PA)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PCLV', PCLV)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PSUPSAT', PSUPSAT)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLCRIT_AER', PLCRIT_AER)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PICRIT_AER', PICRIT_AER)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PRE_ICE', PRE_ICE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PCCN', PCCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PNICE', PNICE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PEXTRA', PEXTRA)

    call ppser_finalize
  end subroutine serialize


  subroutine deserialize( &
   & KLON, KLEV, PTSPHY,&
   & PT, PQ, TENDENCY_CML, TENDENCY_TMP, &
   & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
   & PHRSW, PHRLW, PVERVEL, PAP, PAPH, &
   & PLSM, LDCUM, KTYPE, PLU, PLUDE, PSNDE, PMFU, PMFD, &
   & LDSLPHY, LDMAINCALL, PA, PCLV, PSUPSAT, &
   & PLCRIT_AER, PICRIT_AER, PRE_ICE, &
   & PCCN, PNICE, PEXTRA, KFLDX &
   & )
    ! Deserialization for reading input data from machine-agnostic format

    INTEGER(KIND=JPIM),INTENT(IN)     :: KLON             ! Number of grid points
    INTEGER(KIND=JPIM),INTENT(IN)     :: KLEV             ! Number of levels
    INTEGER(KIND=JPIM),INTENT(IN)     :: KFLDX

    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PTSPHY           ! Physics timestep
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PT(KLON,KLEV)    ! T at start of callpar
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PQ(KLON,KLEV)    ! Q at start of callpar
    TYPE (STATE_TYPE) ,INTENT(OUT)    :: TENDENCY_CML     ! cumulative tendency used for final output
    TYPE (STATE_TYPE) ,INTENT(OUT)    :: TENDENCY_TMP     ! cumulative tendency used as input
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNL(KLON,KLEV) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PDYNI(KLON,KLEV) ! Liq from Dynamics
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PVERVEL(KLON,KLEV) !Vertical velocity
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PAP(KLON,KLEV)   ! Pressure on full levels
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLSM(KLON)       ! Land fraction (0-1) 
    LOGICAL           ,INTENT(OUT)    :: LDCUM(KLON)      ! Convection active
    INTEGER(KIND=JPIM),INTENT(OUT)    :: KTYPE(KLON)      ! Convection type 0,1,2
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLU(KLON,KLEV)   ! Conv. condensate
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLUDE(KLON,KLEV) ! Conv. detrained water 
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSNDE(KLON,KLEV) ! Conv. detrained snow
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
    LOGICAL           ,INTENT(OUT)    :: LDSLPHY 
    LOGICAL           ,INTENT(OUT)    :: LDMAINCALL       ! T if main call to cloudsc
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PCLV(KLON,KLEV,NCLV) 
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSUPSAT(KLON,KLEV)
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLCRIT_AER(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PICRIT_AER(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PRE_ICE(KLON,KLEV) 
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)
    REAL(KIND=JPRB)   ,INTENT(OUT)    :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields

    ! ! Create serializer for storing reference output
    ! ! We need to use prefix_ref to not over-write prvious data.
    ! call ppser_initialize(directory='data', prefix='dummy', prefix_ref='input')
    ! call fs_create_savepoint('input', ppser_savepoint)
    ! call ppser_set_mode(1)

    ! Retrieve parametersand timestep size from the serializer
    ! Note that we use `ppser_serializer_ref` to get the previously stored data
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

    ! Retrieve field data from savepoint
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PT', PT)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PQ', PQ)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_CML_T', TENDENCY_CML%T)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_CML_Q', TENDENCY_CML%Q)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_CML_A', TENDENCY_CML%A)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_CML_CLD', TENDENCY_CML%CLD)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_TMP_T', TENDENCY_TMP%T)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_TMP_Q', TENDENCY_TMP%Q)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_TMP_A', TENDENCY_TMP%A)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'TENDENCY_TMP_CLD', TENDENCY_TMP%CLD)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PVFA', PVFA)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PVFL', PVFL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PVFI', PVFI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PDYNA', PDYNA)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PDYNL', PDYNL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PDYNI', PDYNI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PHRSW', PHRSW)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PHRLW', PHRLW)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PVERVEL', PVERVEL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PAP', PAP)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PAPH', PAPH)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLSM', PLSM)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'LDCUM', LDCUM)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'KTYPE', KTYPE)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLU', PLU)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLUDE', PLUDE)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PSNDE', PSNDE)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PMFU', PMFU)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PMFD', PMFD)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'LDSLPHY', LDSLPHY)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'LDMAINCALL', LDMAINCALL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PA', PA)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PCLV', PCLV)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PSUPSAT', PSUPSAT)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLCRIT_AER', PLCRIT_AER)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PICRIT_AER', PICRIT_AER)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PRE_ICE', PRE_ICE)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PCCN', PCCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PNICE', PNICE)
    ! Note: The 0-sized array (KFLDX=0) seems to create problems when filled with
    ! data from the C-backend, causing memory corruption if enabled.
    ! call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PEXTRA', PEXTRA)

    call ppser_finalize
  end subroutine deserialize

  subroutine serialize_reference( KLON, KLEV, KFLDX, &
       & PLUDE,    PCOVPTOT, PRAINFRAC_TOPRFZ,&
       & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,&
       & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG,&
       & PFSQLTUR, PFSQITUR, PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN)
    ! Serialize reference data for offline validation
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON, KLEV, KFLDX
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: PCOVPTOT(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: PRAINFRAC_TOPRFZ(KLON)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQLF(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQIF(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFCQLNG(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFCQNNG(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQRF(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQSF(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFCQRNG(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFCQSNG(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQLTUR(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFSQITUR(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFPLSL(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFPLSN(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFHPSL(KLON,KLEV+1)
    REAL(KIND=JPRB), INTENT(INOUT) :: PFHPSN(KLON,KLEV+1)

    ! Initialize serializer for storing reference input
    call ppser_initialize(directory='data', prefix='reference')
    call fs_create_savepoint('reference', ppser_savepoint)
    call ppser_set_mode(0)

    ! Store dimensions on the serializer
    call fs_add_serializer_metainfo(ppser_serializer, 'KLON', KLON)
    call fs_add_serializer_metainfo(ppser_serializer, 'KLEV', KLEV)
    call fs_add_serializer_metainfo(ppser_serializer, 'KFLDX', KFLDX)

    ! Store the reference field data
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLUDE', PLUDE)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PCOVPTOT', PCOVPTOT)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PRAINFRAC_TOPRFZ', PRAINFRAC_TOPRFZ)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQLF', PFSQLF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQIF', PFSQIF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFCQLNG', PFCQLNG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFCQNNG', PFCQNNG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQRF', PFSQRF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQSF', PFSQSF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFCQRNG', PFCQRNG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFCQSNG', PFCQSNG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQLTUR', PFSQLTUR)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFSQITUR', PFSQITUR)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFPLSL', PFPLSL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFPLSN', PFPLSN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFHPSL', PFHPSL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PFHPSN', PFHPSN)

    call ppser_finalize
  end subroutine serialize_reference
end module serialize_mod
