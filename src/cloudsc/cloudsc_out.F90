SUBROUTINE CLOUDSC_OUT &
 !---input
 & (KIDIA,    KFDIA,    KLON,    KLEV,&
 & PTSPHY,&
 & PT, PQ, tendency_cml,tendency_tmp,tendency_loc, &
 & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
 & PHRSW,    PHRLW,&
 & PVERVEL,  PAP,      PAPH,&
 & PLSM,     LDCUM,    KTYPE, &
 & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD,&
 & LDSLPHY,  LDMAINCALL, &
 !---prognostic fields
 & PA,&
 & PCLV,  &
 & PSUPSAT,&
!-- arrays for aerosol-cloud interactions
!!! & PQAER,    KAER, &
 & PLCRIT_AER,PICRIT_AER,&
 & PRE_ICE,&
 & PCCN,     PNICE,&
 !---diagnostic output
 & PCOVPTOT, PRAINFRAC_TOPRFZ,&
 !---resulting fluxes
 & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG,&
 & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG,&
 & PFSQLTUR, PFSQITUR , &
 & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN,&
 & PEXTRA,   KFLDX)  

!===============================================================================

USE PARKIND1 , ONLY : JPIM, JPRB
!USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : LSCMEC
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV  
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
USE YOECLDP  , ONLY : YRECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
USE YOMPHYDER ,ONLY : STATE_TYPE
USE YOMJFH   , ONLY : N_VMASS
USE YOEPHLI  , ONLY : YREPHLI

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
! PLCRIT_AER : critical liquid mmr for rain autoconversion process
! PICRIT_AER : critical liquid mmr for snow autoconversion process
! PRE_LIQ : liq Re
! PRE_ICE : ice Re
! PCCN    : liquid cloud condensation nuclei
! PNICE   : ice number concentration (cf. CCN)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PLCRIT_AER(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PICRIT_AER(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_ICE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY            ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)    ! T at start of callpar
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)    ! Q at start of callpar
TYPE (STATE_TYPE) , INTENT (IN)  :: tendency_cml   ! cumulative tendency used for final output
TYPE (STATE_TYPE) , INTENT (IN)  :: tendency_tmp   ! cumulative tendency used as input
TYPE (STATE_TYPE) , INTENT (OUT) :: tendency_loc   ! local tendency from cloud scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNL(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNI(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) !Vertical velocity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)       ! Land fraction (0-1) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON)      ! Convection active
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)      ! Convection type 0,1,2
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV)   ! Conv. condensate
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) ! Conv. detrained water 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNDE(KLON,KLEV) ! Conv. detrained snow
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
LOGICAL           ,INTENT(IN)    :: LDSLPHY 
LOGICAL           ,INTENT(IN)    :: LDMAINCALL       ! T if main call to cloudsc
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLV(KLON,KLEV,NCLV) 

 ! Supersat clipped at previous time level in SLTEND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSUPSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOVPTOT(KLON,KLEV) ! Precip fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAINFRAC_TOPRFZ(KLON) 
! Flux diagnostics for DDH budget
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLF(KLON,KLEV+1)  ! Flux of liquid
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQIF(KLON,KLEV+1)  ! Flux of ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLNG(KLON,KLEV+1) ! -ve corr for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQNNG(KLON,KLEV+1) ! -ve corr for ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQRF(KLON,KLEV+1)  ! Flux diagnostics
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQSF(KLON,KLEV+1)  !    for DDH, generic
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQRNG(KLON,KLEV+1) ! rain
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQSNG(KLON,KLEV+1) ! snow
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLTUR(KLON,KLEV+1) ! liquid flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQITUR(KLON,KLEV+1) ! ice flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSL(KLON,KLEV+1) ! liq+rain sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSN(KLON,KLEV+1) ! ice+snow sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL(KLON,KLEV+1) ! Enthalpy flux for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN(KLON,KLEV+1) ! Enthalp flux for ice

INTEGER(KIND=JPIM), parameter :: iu = 123

#include "abor1.intfb.h"

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('CLOUDSC_OUT',0,ZHOOK_HANDLE)

open(iu,file='cloudsc.bin',status='new',&
     & access='stream', form='unformatted')

write(iu) KLON,KLEV,KFLDX
write(iu) PTSPHY
write(iu) LDSLPHY 
write(iu) LDMAINCALL

write(iu) LSCMEC
write(iu) RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV  
write(iu) R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
write(iu) NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV

!write(iu) YRECLDP ! can't use this or byte-swapping may or may not be performed upon read                                                                           
! TYPE :: TECLDP -- must write each individual member                                                                                                                
WRITE(IU) YRECLDP%RAMID
WRITE(IU) YRECLDP%RCLDIFF
WRITE(IU) YRECLDP%RCLDIFF_CONVI
WRITE(IU) YRECLDP%RCLCRIT
WRITE(IU) YRECLDP%RCLCRIT_SEA
WRITE(IU) YRECLDP%RCLCRIT_LAND
WRITE(IU) YRECLDP%RKCONV
WRITE(IU) YRECLDP%RPRC1
WRITE(IU) YRECLDP%RPRC2
WRITE(IU) YRECLDP%RCLDMAX
WRITE(IU) YRECLDP%RPECONS
WRITE(IU) YRECLDP%RVRFACTOR
WRITE(IU) YRECLDP%RPRECRHMAX
WRITE(IU) YRECLDP%RTAUMEL
WRITE(IU) YRECLDP%RAMIN
WRITE(IU) YRECLDP%RLMIN
WRITE(IU) YRECLDP%RKOOPTAU
WRITE(IU) YRECLDP%RCLDTOPP
WRITE(IU) YRECLDP%RLCRITSNOW
WRITE(IU) YRECLDP%RSNOWLIN1
WRITE(IU) YRECLDP%RSNOWLIN2
WRITE(IU) YRECLDP%RICEHI1
WRITE(IU) YRECLDP%RICEHI2
WRITE(IU) YRECLDP%RICEINIT
WRITE(IU) YRECLDP%RVICE
WRITE(IU) YRECLDP%RVRAIN
WRITE(IU) YRECLDP%RVSNOW
WRITE(IU) YRECLDP%RTHOMO
WRITE(IU) YRECLDP%RCOVPMIN
WRITE(IU) YRECLDP%RCCN
WRITE(IU) YRECLDP%RNICE
WRITE(IU) YRECLDP%RCCNOM
WRITE(IU) YRECLDP%RCCNSS
WRITE(IU) YRECLDP%RCCNSU
WRITE(IU) YRECLDP%RCLDTOPCF
WRITE(IU) YRECLDP%RDEPLIQREFRATE
WRITE(IU) YRECLDP%RDEPLIQREFDEPTH
WRITE(IU) YRECLDP%RCL_KKAac
WRITE(IU) YRECLDP%RCL_KKBac
WRITE(IU) YRECLDP%RCL_KKAau
WRITE(IU) YRECLDP%RCL_KKBauq
WRITE(IU) YRECLDP%RCL_KKBaun
WRITE(IU) YRECLDP%RCL_KK_cloud_num_sea
WRITE(IU) YRECLDP%RCL_KK_cloud_num_land
WRITE(IU) YRECLDP%RCL_AI
WRITE(IU) YRECLDP%RCL_BI
WRITE(IU) YRECLDP%RCL_CI
WRITE(IU) YRECLDP%RCL_DI
WRITE(IU) YRECLDP%RCL_X1I
WRITE(IU) YRECLDP%RCL_X2I
WRITE(IU) YRECLDP%RCL_X3I
WRITE(IU) YRECLDP%RCL_X4I
WRITE(IU) YRECLDP%RCL_CONST1I
WRITE(IU) YRECLDP%RCL_CONST2I
WRITE(IU) YRECLDP%RCL_CONST3I
WRITE(IU) YRECLDP%RCL_CONST4I
WRITE(IU) YRECLDP%RCL_CONST5I
WRITE(IU) YRECLDP%RCL_CONST6I
WRITE(IU) YRECLDP%RCL_APB1
WRITE(IU) YRECLDP%RCL_APB2
WRITE(IU) YRECLDP%RCL_APB3
WRITE(IU) YRECLDP%RCL_AS
WRITE(IU) YRECLDP%RCL_BS
WRITE(IU) YRECLDP%RCL_CS
WRITE(IU) YRECLDP%RCL_DS
WRITE(IU) YRECLDP%RCL_X1S
WRITE(IU) YRECLDP%RCL_X2S
WRITE(IU) YRECLDP%RCL_X3S
WRITE(IU) YRECLDP%RCL_X4S
WRITE(IU) YRECLDP%RCL_CONST1S
WRITE(IU) YRECLDP%RCL_CONST2S
WRITE(IU) YRECLDP%RCL_CONST3S
WRITE(IU) YRECLDP%RCL_CONST4S
WRITE(IU) YRECLDP%RCL_CONST5S
WRITE(IU) YRECLDP%RCL_CONST6S
WRITE(IU) YRECLDP%RCL_CONST7S
WRITE(IU) YRECLDP%RCL_CONST8S
WRITE(IU) YRECLDP%RDENSWAT
WRITE(IU) YRECLDP%RDENSREF
WRITE(IU) YRECLDP%RCL_AR
WRITE(IU) YRECLDP%RCL_BR
WRITE(IU) YRECLDP%RCL_CR
WRITE(IU) YRECLDP%RCL_DR
WRITE(IU) YRECLDP%RCL_X1R
WRITE(IU) YRECLDP%RCL_X2R
WRITE(IU) YRECLDP%RCL_X4R
WRITE(IU) YRECLDP%RCL_KA273
WRITE(IU) YRECLDP%RCL_CDENOM1
WRITE(IU) YRECLDP%RCL_CDENOM2
WRITE(IU) YRECLDP%RCL_CDENOM3
WRITE(IU) YRECLDP%RCL_SCHMIDT
WRITE(IU) YRECLDP%RCL_DYNVISC
WRITE(IU) YRECLDP%RCL_CONST1R
WRITE(IU) YRECLDP%RCL_CONST2R
WRITE(IU) YRECLDP%RCL_CONST3R
WRITE(IU) YRECLDP%RCL_CONST4R
WRITE(IU) YRECLDP%RCL_FAC1
WRITE(IU) YRECLDP%RCL_FAC2
WRITE(IU) YRECLDP%RCL_CONST5R
WRITE(IU) YRECLDP%RCL_CONST6R
WRITE(IU) YRECLDP%RCL_FZRAB
WRITE(IU) YRECLDP%RCL_FZRBB
WRITE(IU) YRECLDP%LCLDEXTRA, YRECLDP%LCLDBUDGET
WRITE(IU) YRECLDP%NSSOPT
WRITE(IU) YRECLDP%NCLDTOP
WRITE(IU) YRECLDP%NAECLBC, YRECLDP%NAECLDU, YRECLDP%NAECLOM, YRECLDP%NAECLSS, YRECLDP%NAECLSU
WRITE(IU) YRECLDP%NCLDDIAG
WRITE(IU) YRECLDP%NAERCLD
WRITE(IU) YRECLDP%LAERLIQAUTOLSP
WRITE(IU) YRECLDP%LAERLIQAUTOCP
WRITE(IU) YRECLDP%LAERLIQAUTOCPB
WRITE(IU) YRECLDP%LAERLIQCOLL
WRITE(IU) YRECLDP%LAERICESED
WRITE(IU) YRECLDP%LAERICEAUTO
WRITE(IU) YRECLDP%NSHAPEP
WRITE(IU) YRECLDP%NSHAPEQ
WRITE(IU) YRECLDP%NBETA
WRITE(IU) YRECLDP%RBETA(0:100)
WRITE(IU) YRECLDP%RBETAP1(0:100)
!END TYPE TECLDP

!write(iu) YREPHLI! can't use this or byte-swapping may or may not be performed upon read                                                                            
! TYPE :: TECLDP                                                                                                                                                     
! TYPE :: TEPHLI -- must write each individual member                                                                                                                
WRITE(IU) YREPHLI%LTLEVOL
WRITE(IU) YREPHLI%LPHYLIN
WRITE(IU) YREPHLI%LENOPERT
WRITE(IU) YREPHLI%LEPPCFLS
WRITE(IU) YREPHLI%LRAISANEN
WRITE(IU) YREPHLI%RLPTRC
WRITE(IU) YREPHLI%RLPAL1
WRITE(IU) YREPHLI%RLPAL2
WRITE(IU) YREPHLI%RLPBB
WRITE(IU) YREPHLI%RLPCC
WRITE(IU) YREPHLI%RLPDD
WRITE(IU) YREPHLI%RLPMIXL
WRITE(IU) YREPHLI%RLPBETA
WRITE(IU) YREPHLI%RLPDRAG
WRITE(IU) YREPHLI%RLPEVAP
WRITE(IU) YREPHLI%RLPP00
!END TYPE TEPHLI

write(iu) N_VMASS

!USE YOMPHYDER ,ONLY : STATE_TYPE
!type state_type
!  REAL(KIND=JPRB), dimension(:,:), pointer :: u,v,T   ! GMV fields
!  REAL(KIND=JPRB), dimension(:,:), pointer :: o3,q,a  ! GFL fields
!  REAL(KIND=JPRB), dimension(:,:,:), pointer :: cld   ! composed cloud array
!  !REAL(KIND=JPRB), dimension(:,:), pointer :: qsat    ! spec. humidity at saturation
!end type state_type

write(iu) PLCRIT_AER 
write(iu) PICRIT_AER 
write(iu) PRE_ICE 
write(iu) PCCN
write(iu) PNICE
write(iu) PT
write(iu) PQ
write(iu) tendency_cml%u,tendency_cml%v,tendency_cml%T, &
     & tendency_cml%o3,tendency_cml%q,tendency_cml%a, &
     & tendency_cml%cld
write(iu) tendency_tmp%u,tendency_tmp%v,tendency_tmp%T, &
     & tendency_tmp%o3,tendency_tmp%q,tendency_tmp%a, &
     & tendency_tmp%cld
write(iu) PVFA
write(iu) PVFL
write(iu) PVFI
write(iu) PDYNA
write(iu) PDYNL
write(iu) PDYNI
write(iu) PHRSW
write(iu) PHRLW
write(iu) PVERVEL
write(iu) PAP
write(iu) PAPH
write(iu) PLSM
write(iu) LDCUM
write(iu) KTYPE
write(iu) PLU
write(iu) PLUDE
write(iu) PSNDE
write(iu) PMFU
write(iu) PMFD
write(iu) PA
write(iu) PEXTRA
write(iu) PCLV
write(iu) PSUPSAT

close(iu)
call abor1('Captured CLOUDSC output')

!IF (LHOOK) CALL DR_HOOK('CLOUDSC_OUT',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUDSC_OUT
