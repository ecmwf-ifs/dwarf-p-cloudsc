! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module yomphyder

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE
SAVE

type dimension_type
 INTEGER(KIND=JPIM), pointer :: KIDIA,KFDIA,KBL,KLON,KSMAX,KLEVS,KTILES,KSTGLO,&
 & KLEVX,KFLDX,KFLDX2,KTDIA,KLEV,KVCLIS,KLEVI,KLEVO,KLEVSN,&
 & KDHVTLS,KDHFTLS,KDHVTSS,KDHFTSS,KDHVTTS,KDHFTTS,KDHVTIS,&
 & KDHFTIS,KDHVSSS,KDHFSSS,KDHVIIS,KDHFIIS,KDHVWLS,KDHFWLS
 INTEGER(KIND=JPIM), dimension(:), pointer  :: KGPLAT  ! DM-global number of the latitude of point jrof=KSTART
end type dimension_type

type state_type
  REAL(KIND=JPRB), dimension(:,:), pointer :: u,v,T   ! GMV fields
  REAL(KIND=JPRB), dimension(:,:), pointer :: o3,q,a  ! GFL fields
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: cld   ! composed cloud array
  !REAL(KIND=JPRB), dimension(:,:), pointer :: qsat    ! spec. humidity at saturation
end type state_type

type model_state_type
  REAL(KIND=JPRB), dimension(:,:), pointer ::  U, V, T ! GMV
  REAL(KIND=JPRB), dimension(:,:,:), pointer ::  GFL   ! GFL
end type model_state_type

type mask_GFL_type
  LOGICAL, pointer :: o3,q,a,ql,qi,qr,qs
end type mask_GFL_type

type aux_type
  REAL(KIND=JPRB), dimension(:,:), pointer :: PAPHI,PAPHIF ! half and full level geopotential
  REAL(KIND=JPRB), dimension(:,:), pointer :: PGEOM1,PGEOMH ! half and full level distance from surface (in gpm)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PAPRS,PAPRSF ! half and full level pressure
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDELP        ! layer thickness  (in pressure units)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PRS1,PRSF1   ! provisional t+dt pressure on half and full levels
  REAL(KIND=JPRB), dimension(:,:), pointer :: PVERVEL      ! diagnostic vertical velocity
  !  geometry
  REAL(KIND=JPRB), dimension(:), pointer   :: PGELAM, PGELAT  ! LONGITUDE, LATITUDE (RADIANS)
  REAL(KIND=JPRB), dimension(:), pointer   :: PGAW         ! GAUSSIAN WEIGHTS - Reduced Grid - ~ grid box area
  REAL(KIND=JPRB), dimension(:), pointer   :: PCLON, PSLON ! cosine, sine of longitude
  REAL(KIND=JPRB), dimension(:), pointer   :: PMU0, PMU0M  ! local cosine of instantaneous (mean) solar zenith angle
  REAL(KIND=JPRB), dimension(:), pointer   :: PGEMU        ! sine of latitude 
  REAL(KIND=JPRB), dimension(:), pointer   :: POROG        ! orography
  REAL(KIND=JPRB), dimension(:), pointer   :: PGNORDL,PGNORDM ! Longitudial/latitudial derivatives of orography
  REAL(KIND=JPRB), dimension(:), pointer   :: PGSQM2
end type aux_type

type surf_and_more_type
  REAL(KIND=JPRB), dimension(:,:), pointer :: PSP_SG, PSP_SL, PSP_RR, PSP_X2, PSD_WS, PSD_VF, &
    &                                         PSD_VN, PSD_V2, PSD_VD, PSD_X2, PSD_WW
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PSP_OM, PSP_SB, PSP_EP, PSD_V3, PSD_XA
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PEXTRD 
  REAL(KIND=JPRB), dimension(:,:), pointer :: PCOVPTOT  !Precip fraction
  REAL(KIND=JPRB), dimension(:), pointer :: PQCFL
  ! T star tiles
  REAL(KIND=JPRB), dimension(:,:), pointer :: PUSTRTI  ! E-W (INSTANTANEOUS) SURFACE STRESS FOR EACH TILE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PVSTRTI  ! N-S (INSTANTANEOUS) SURFACE STRESS FOR EACH TILE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PAHFSTI  ! (INSTANTANEOUS) SURFACE SENSIBLE HEAT FLUX FOR EACH TILE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PEVAPTI  ! (INSTANTANEOUS) EVAPORATION FOR EACH TILE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTSKTI   ! SKIN TEMPERATURE FOR EACH TILE
  ! other 
  REAL(KIND=JPRB), dimension(:), pointer ::  PEMIS      ! MODEL SURFACE LONGWAVE EMISSIVITY.
  ! GPP/REC flux adjustment coefficients
  REAL(KIND=JPRB), dimension(:), pointer :: PCGPP, PCREC ! to store bias correction coefficients
  !  Tendencies from surface scheme
  REAL(KIND=JPRB), dimension(:), pointer :: PSNSE1     ! OF SNOW MASS PER UNIT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PASNE1     ! OF SNOW ALBEDO
  REAL(KIND=JPRB), dimension(:), pointer :: PRSNE1     ! OF SNOW DENSITY
  REAL(KIND=JPRB), dimension(:), pointer :: PTSNE1     ! OF SNOW TEMPERATURE
  REAL(KIND=JPRB), dimension(:), pointer :: PWLE1      ! OF SKIN RESERVOIR WATER CONTENT
  REAL(KIND=JPRB), dimension(:), pointer :: PTLE1      ! OF SKIN TEMPERATURE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTSAE1   ! OF MULTI-LAYER SOIL TEMPERATURE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PWSAE1   ! OF MULTI-LAYER SOIL WETNESS
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTIAE1   ! OF MULTI-LAYER SEA ICE TEMPERATURE
  REAL(KIND=JPRB), dimension(:,:), pointer :: PUOE1    ! VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2
  REAL(KIND=JPRB), dimension(:,:), pointer :: PVOE1    ! VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTOE1    ! TEMPERATURE TENDENCY OF OCEAN MIXED LAYER MODEL  K/S
  REAL(KIND=JPRB), dimension(:,:), pointer :: PSOE1    ! SALINITY TENDENCY OF OCEAN MIXED LAYER MODEL  psu/S
  REAL(KIND=JPRB), dimension(:), pointer :: PHSTD      ! STANDARD DEVIATION OF SUBGRID OROGRAPHY
  INTEGER(KIND=JPIM), dimension(:), pointer :: ITVH, ITVL, ISOTY ! to store vegetation indices and soil type index
  REAL(KIND=JPRB), dimension(:), pointer :: PCVL, PCVH       ! to store vegetation covers
  REAL(KIND=JPRB), dimension(:), pointer :: PLAIL, PLAIH     ! to store variable LAI
  REAL(KIND=JPRB), dimension(:), pointer :: PEVAPMU    ! potential evaporation
  REAL(KIND=JPRB), dimension(:),  pointer  :: PTLICEE1       ! tendency of lake ice temperature
  REAL(KIND=JPRB), dimension(:),  pointer  :: PTLMNWE1       ! tendency of lake totat layer temperature
  REAL(KIND=JPRB), dimension(:),  pointer  :: PTLWMLE1       ! tendency of lake mixed layer temperature
  REAL(KIND=JPRB), dimension(:),  pointer  :: PTLBOTE1       ! tendency of lake bottom layer temperature
  REAL(KIND=JPRB), dimension(:),  pointer  :: PTLSFE1        ! tendency of lake shape factor - 
  REAL(KIND=JPRB), dimension(:),  pointer  :: PHLICEE1       ! tendency of lake ice depth m
  REAL(KIND=JPRB), dimension(:),  pointer  :: PHLMLE1        ! tendency of lake mixed layer depth m/s
end type surf_and_more_type

type perturb_in_type
  REAL(KIND=JPRB), dimension(:), pointer ::  PSTOPHU,PSTOPHV,PSTOPHT,PSTOPHQ ! random number for defining stochastic
                                                                !  perturbation for U, V, T, and Q diabatic tendency.
  REAL(KIND=JPRB), dimension(:), pointer ::  PSTOPHCA  ! CA pattern 
  REAL(KIND=JPRB), dimension(:,:), pointer :: PGP2DSDT
  REAL(KIND=JPRB), dimension(:,:), pointer :: PVORT, PVORTGRADX, PVORTGRADY ! vorticity and its horizontal gradients
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTOTDISS_SMOOTH ! smoothed total dissipation rate
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFORCEU, PFORCEV, PFORCET, PFORCEQ ! nonlinear stochastic forcing terms 
end type perturb_in_type


type perturb_inout_type
  REAL(KIND=JPRB), dimension(:), pointer ::  PCUCONVCA  ! CA state
  REAL(KIND=JPRB), dimension(:), pointer ::  PNLCONVCA  ! indicates if CA was alive at a gridpoint before
  REAL(KIND=JPRB), dimension(:), pointer ::  PCAPECONVCA ! 0 or NLIVES dependeing on threshold (CIN defined at the moment)
  REAL(KIND=JPRB), dimension(:,:), pointer ::  PSTREAM  ! Streamfunction forcing/spectral pattern
  REAL(KIND=JPRB), dimension(:,:), pointer ::  PTEMP    ! Temperature forcing/spectral pattern
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTOTDISS  ! total dissipation rate
end type perturb_inout_type

type aux_rad_type
  REAL(KIND=JPRB), dimension(:,:), pointer :: PNEB     ! FRACTIONAL CLOUDINESS FOR RADIATION.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PQICE    ! (KLON,KLEV) SPECIFIC HUMIDITY OF SOLID WATER FOR RADIATION.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PQLI     ! (KLON,KLEV) SPECIFIC HUMIDITY OF LIQUID WATER FOR RADIATION
  REAL(KIND=JPRB), dimension(:,:), pointer :: PEMTD    ! (KLON,0:KLEV) Net longwave flux
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTRSW    ! (KLON,0:KLEV) SW transmissivity (scale by incoming SW for flux)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PEMTC    ! (KLON,0:KLEV)  Clear-sky net longwave flux
  REAL(KIND=JPRB), dimension(:,:), pointer :: PTRSC    ! (KLON,0:KLEV)  Clear-sky shortwave transmissivity
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PTAUAER ! (KLON,KLEV,6)  OPTICAL THICKNESS FOR 6 AEROSOL TYPES
  REAL(KIND=JPRB), dimension(:), pointer :: PSRLWD     ! (KLON)  Surface downward longwave flux
  REAL(KIND=JPRB), dimension(:), pointer :: PSRLWDC    ! (KLON)  SURFACE DOWNWARD CLEAR-SKY LONGWAVE 
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWD     ! (KLON)  SURFACE SHORTWAVE DOWNWARDS
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWDC    ! (KLON)  SURFACE DOWNWARD CLEAR-SKY SHORTWAVE
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWDCS   ! (KLON)  SURFACE NET SHORTWAVE CLEAR-SKY
  REAL(KIND=JPRB), dimension(:), pointer :: PEDRO      ! (KLON)  STORED SKIN (BRIGHTNESS) TEMPERATURE FROM LAST RAD. STEP
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWPAR   ! (KLON) downward SW PAR radiation at the surface
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWUVB   ! (KLON) downward UV-B radiation at the surface
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWPARC  ! (KLON) downward clear-sky SW PAR radiation at the surface
  REAL(KIND=JPRB), dimension(:), pointer :: PSRSWTINC  ! (KLON) TOA incident solar radiation
  REAL(KIND=JPRB), dimension(:), pointer :: PSRFDIR    ! (KLON)  total sky direct downward SW radiation
  REAL(KIND=JPRB), dimension(:), pointer :: PSRCDIR    ! (KLON)  clear-sky direct downward SW radiation
  REAL(KIND=JPRB), dimension(:), pointer :: PTINCF     ! TOA INCIDENT SOLAR RADIATION
  REAL(KIND=JPRB), dimension(:), pointer :: PFDIR, PCDIR   ! TOTAL (clear-SKY) DIRECT SW RADIATION AT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PISUND     ! SUNSHINE DURATION (FRACTION OF TIME STEP)
  REAL(KIND=JPRB), dimension(:), pointer :: PINCSR,PINCSR0 ! INCIDENT AND ANUAL INCIDENT SOLAR RADIATION
  REAL(KIND=JPRB), dimension(:,:), pointer :: PHRSW,PHRLW,PHRSC,PHRLC ! SW/LW all-sky/clear-sky heating rates
  REAL(KIND=JPRB), dimension(:,:), pointer :: DerivativeLw ! (KLON,0:KLEV) Derivative to update LW radiation each timestep

end type aux_rad_type

! Other type characteristics
type aux_diag_type
  REAL(KIND=JPRB), dimension(:), pointer :: PCAPE      !  CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
  REAL(KIND=JPRB), dimension(:), pointer :: PCBASE     !  CLOUD BASE
  REAL(KIND=JPRB), dimension(:), pointer :: P0DEGL     ! 0 DEGREE CELSIUS LEVEL
  REAL(KIND=JPRB), dimension(:), pointer :: PVISIH     ! HORIZONTAL VISIBILITY
  REAL(KIND=JPRB), dimension(:), pointer :: PCIN       ! CIN
  REAL(KIND=JPRB), dimension(:), pointer :: PVDIS      ! TURBULENT DISSIPATION
  REAL(KIND=JPRB), dimension(:), pointer :: PVDISG     ! GRAVITY WAVE DRAG DISSIPATION
  REAL(KIND=JPRB), dimension(:), pointer :: PUSTRG     ! GRAVITY WAVE DRAG SURFACE U-STRESS
  REAL(KIND=JPRB), dimension(:), pointer :: PVSTRG     ! GRAVITY WAVE DRAG SURFACE V-STRESS
  REAL(KIND=JPRB), dimension(:), pointer :: PI10FG     ! WIND GUST AT 10 m FOR CURRENT TIME LEVEL (m/s)
  REAL(KIND=JPRB), dimension(:), pointer :: PPRECTYPE  ! Precipitation type
  REAL(KIND=JPRB), dimension(:), pointer :: PFZRA      ! Freezing rain rate
  !INTEGER(KIND=JPIM),dimension(:), pointer :: KTYPE    ! Convection type (0,1,2)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PCONVIND ! CONVECTIVE INDICES
  REAL(KIND=JPRB), dimension(:,:), pointer :: PQSAT    ! SPECIFIC HUMIDITY AT SATURATION.
  ! Stored tendencies before convection (used in 1D var)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PZTENT, PZTENQ
  !    3D DIAGNOSTICS FOR ERA40
  REAL(KIND=JPRB), dimension(:,:), pointer :: PMFUDE_RATE   ! UD detrainmnet rate (KG/(M3*S))
  REAL(KIND=JPRB), dimension(:,:), pointer :: PMFDDE_RATE   ! DD detrainmnet rate (KG/(M3*S))
  REAL(KIND=JPRB), dimension(:,:), pointer :: PKH_VDF       ! turbulent diffusion coefficient for heat 
  !    array for precipitation fraction
  REAL(KIND=JPRB), dimension(:,:), pointer :: PCOVPTOT     !  PRECIPITATION FRACTION IN EACH LAYER
  ! Convection and PBL types
  INTEGER(KIND=JPIM), dimension(:), pointer :: ITYPE, ICBOT, ICTOP
  INTEGER(KIND=JPIM), dimension(:), pointer :: IPBLTYPE
  ! Convection to cloud communication
  REAL(KIND=JPRB), dimension(:,:), pointer :: PMFU,PMFD    !  Conv. mass flux up/down
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZLUDE      ! DETRAINED LIQUID WATER
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZLU        ! LIQUID WATER CONTENT IN UPDRAFTS
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZSNDE      ! DETRAINED SNOW
end type aux_diag_type

type flux_type
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFRSOC,PFRTHC ! (KLON,0:1) SHORTWAVE and LONGWAVE CLEAR SKY RADIATIVE FLUX
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDIFCQ     ! CONVECTIVE FLUX OF SPECIFIC HUMIDITY (NOT RAIN/SNOW).
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDIFCS     ! CONVECTIVE FLUX OF ENTHALPY (NOT RAIN/SNOW).
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDIFTQ     ! TURBULENT FLUX (INC. Q NEGATIVE) OF SPECIFIC HUMIDITY.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDIFTI,PDIFTL ! TURBULENT FLUX (INC. Q NEGATIVE) OF SOLID/LIQUID WATER
  REAL(KIND=JPRB), dimension(:,:), pointer :: PDIFTS     ! TURBULENT FLUX OF ENTHALPY (OR DRY STATIC ENERGY).
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFCQNG     ! PSEUDO-FLUX OF WATER TO CORRECT FOR Q<0.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFCQLNG,PFCQNNG  ! PSEUDO-FLUX OF LIQUID WATER/ICE TO CORRECT FOR Q<0.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFCQRNG,PFCQSNG  ! PSEUDO-FLUX OF RAIN/SNOW TO CORRECT FOR Q<0.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFPLCL     ! CONVECTIVE PRECIPITATION AS RAIN.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFPLCN     ! CONVECTIVE PRECIPITATION AS SNOW.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFPLSL     ! STRATIFORM PRECIPITATION AS RAIN.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFPLSN     ! STRATIFORM PRECIPITATION AS SNOW.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFHPCL     ! ENTHALPY FLUX OF CONVECTIVE PRECIPITATION AS RAIN
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFHPCN     ! ENTHALPY FLUX OF CONVECTIVE PRECIPITATION AS SNOW.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFHPSL     ! ENTHALPY FLUX OF STRATIFORM PRECIPITATION AS RAIN.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFHPSN     ! ENTHALPY FLUX OF STRATIFORM PRECIPITATION AS SNOW.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFRSO      ! SHORTWAVE RADIATIVE FLUX.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFRTH      ! LONGWAVE RADIATIVE FLUX.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PSTRCU,PSTRCV ! CONVECTIVE FLUX OF MOMENTUM (U and V)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PSTRDU,PSTRDV ! GRAVITY WAVE DRAG FLUX OF MOMENTUM (U and V)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PSTRTU,PSTRTV ! TURBULENT FLUX OF MOMENTUM (U and V)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFCCQL,PFCCQN ! CONVECTIVE CONDENSATION FLUX FOR LIQUID WATER/ICE.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFCSQL,PFCSQN ! STRATIFORM CONDENSATION FLUX FOR LIQUID WATER/ICE.
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFSQRF,PFSQSF ! RAIN/SNOW FLUX FROM LS CLOUD SCHEME
  REAL(KIND=JPRB), dimension(:,:), pointer :: PFSQLTUR,PFSQITUR ! LIQUID/ICE FLUX DUE TO VDF
  REAL(KIND=JPRB), dimension(:), pointer :: PFRSOD,PFRSODC  ! DOWNWARD (CLEAR-SKY) SHORTWAVE FLUX AT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PFRTHD,PFRTHDC  ! DOWNWARD (CLEAR-SKY) LONGWAVE FLUX AT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PUVDF,PPARF     ! DOWNWARD UV,PAR FLUX AT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PPARCF          ! CLEAR-SKY DOWNWARD PAR FLUX AT SURFACE
  REAL(KIND=JPRB), dimension(:), pointer :: PTOFDU,PTOFDV   ! TOFD COMP. OF TURBULENT FLUX OF U/V-MOMEMTUM
  REAL(KIND=JPRB), dimension(:), pointer :: PFTG12          ! HEAT FLUX FROM LAYER 1 TO 2  ! NOT USED
  REAL(KIND=JPRB), dimension(:), pointer :: PFWG12          ! WATER FLUX FROM LAYER 1 TO 2  ! NOT USED
  REAL(KIND=JPRB), dimension(:), pointer :: PFWEV           ! SURFACE WATER VAPOUR FLUX (EVAPORATION)
  REAL(KIND=JPRB), dimension(:), pointer :: PFWSB           ! SURFACE WATER VAPOUR FLUX (SUBLIMATION)
  REAL(KIND=JPRB), dimension(:), pointer :: PFWMLT          ! WATER FLUX CORRESPONDING TO SURFACE SNOW MELT.
  REAL(KIND=JPRB), dimension(:), pointer :: PFWROD          ! RUN-OFF FLUX AT DEEPER LAYERS (2-4)
  REAL(KIND=JPRB), dimension(:), pointer :: PFWRO1          ! RUN-OFF FLUX AT LAYER 1
  REAL(KIND=JPRB), dimension(:), pointer :: PFTLHEV         ! SURFACE LATENT HEAT FLUX (SNOW FREE FRACTION)
  REAL(KIND=JPRB), dimension(:), pointer :: PFTLHSB         ! SURFACE LATENT HEAT FLUX (SNOW COVERED FRACTION)
  REAL(KIND=JPRB), dimension(:,:), pointer :: PENES         !  ???
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PAERODDF    ! Diagnostic array with aerosol fluxes (sources and sinks)
end type flux_type

type ddh_surf_type
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHTLS      ! Diagnostic array for tiles (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHTSS      ! Diagnostic array for snow T (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHTTS      ! Diagnostic array for soil T (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHTIS      ! Diagnostic array for ice T (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHSSS      ! Diagnostic array for snow mass (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:),   pointer :: PDHIIS      ! Diagnostic array for interception layer (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: PDHWLS      ! Diagnostic array for soil water (see module yomcdh)
  REAL(KIND=JPRB), dimension(:,:),   pointer :: PDHFSS
end type ddh_surf_type

type gems_local_type
  INTEGER(KIND=JPIM), pointer                :: ITRAC
  INTEGER(KIND=JPIM), dimension(:),  pointer :: ITRACER    ! Index array for flexible tracers
  INTEGER(KIND=JPIM), dimension(:),  pointer :: IAERO      ! Index array for aerosol tracers
  INTEGER(KIND=JPIM), dimension(:),  pointer :: ICHEM      ! Index array for chemical species
  REAL(KIND=JPRB), dimension(:,:) ,  pointer :: ZLRCH4     ! temporal storage for AD code
!-- basic GEMS fields
  REAL(KIND=JPRB), dimension(:,:,:), POINTER :: ZCEN,ZTENC,ZTENC_SKF
  REAL(KIND=JPRB), dimension(:,:)  , POINTER :: ZCFLX,ZCHEMDV
  REAL(KIND=JPRB), dimension(:)    , POINTER :: ZSCAV
! - aerosols transport arrays
  REAL(KIND=JPRB), dimension(:,:,:), POINTER :: ZTAUAER, ZAERFLX, ZAERLISI, ZCAERO, ZAEROP
  REAL(KIND=JPRB), dimension(:,:)  , POINTER :: ZAERSRC, ZAERDDP, ZAERSDM, ZAERLIF
  REAL(KIND=JPRB), dimension(:)    , POINTER :: ZAZ0M, ZAZ0H
! - DMS-related aerosol local arrays
  REAL(KIND=JPRB), dimension(:)    , POINTER :: ZDMSO, ZLDAY, ZLISS, ZSO2, ZTDMS, ZDMSI, ZODMS
!-- other prognostic aerosol-related quantities
  REAL(KIND=JPRB), dimension(:)    , pointer ::  ZAERWS, ZAERGUST, ZAERUST, ZDIST
  REAL(KIND=JPRB), dimension(:,:)  , pointer ::  ZAERMAP ! (KDIM%KLON,5)
!-- combined aerosol optical properties at 2 wavelenghts
  REAL(KIND=JPRB), dimension(:,:,:), pointer ::  ZAERTAULT, ZAERTAULB,ZAERASYL ! (KDIM%KLON,KDIM%KLEV,2)

! Visibility
  REAL(KIND=JPRB), dimension(:)    , pointer :: ZCLAERS,  ZPRAERS ! zpraers not used for the moment
  REAL(KIND=JPRB), dimension(:)    , pointer :: ZVISICL
! Preparation for chemistry
  REAL(KIND=JPRB), dimension(:,:,:), pointer :: ZKOZO ! PHOTOCHEMICAL COEFFICIENTS COMPUTED FROM A 2D PHOTOCHEMICAL MODEL
end type gems_local_type

type perturb_local_type
 REAL(KIND=JPRB), dimension(:)  , pointer :: ZWMEAN   ! vertically averaged conv. updraught speed (M/S)
 REAL(KIND=JPRB), dimension(:,:), pointer :: ZDISSCU  ! Dissipation (J/M2) (convection)
 REAL(KIND=JPRB), dimension(:)  , pointer :: ZCUCONVCA ! scaled CA to pass to convection
 REAL(KIND=JPRB), dimension(:,:), pointer :: ZDISSGW  ! kinetic energy/unit mass dissipated in one t-step (gravity wave drag)
end type perturb_local_type

type surf_and_more_local_type
   REAL(KIND=JPRB), dimension(:),  pointer  :: ZLAILI, ZLAIHI   ! Arrays for interactive LAI (CITESSEL option)
   REAL(KIND=JPRB), dimension(:),  pointer  :: ZSNM1M, ZRSNM1M, ZWLM1M
   REAL(KIND=JPRB), dimension(:),  pointer  :: ZHSDFOR
   REAL(KIND=JPRB), dimension(:),  pointer  :: ZPCLAKE           ! lake fraction
   REAL(KIND=JPRB), dimension(:),  pointer  :: ZTHKICE, ZSNTICE  ! sea ice
   REAL(KIND=JPRB), dimension(:),  pointer  :: PCGPP, PCREC ! to store flux adjustment coefficients
   REAL(KIND=JPRB), dimension(:),  pointer  ::  ZWLMX
   REAL(KIND=JPRB), dimension(:),  pointer  ::  ZEMIR, ZEMIW
   REAL(KIND=JPRB), dimension(:),  pointer  ::  ZEVAPSNW
   REAL(KIND=JPRB), dimension(:,:),  pointer  :: ZFRTI,ZALBTI ! KLON,KTILES
   REAL(KIND=JPRB), dimension(:,:),  pointer  :: ZFRSOTI
   REAL(KIND=JPRB), dimension(:,:),  pointer  :: ZAHFTRTI
   REAL(KIND=JPRB), dimension(:,:),  pointer  :: ZALBD,ZALBP ! KLON,NTSW
   !  CTESSEL: Carbon model 
   REAL(KIND=JPRB), dimension(:,:),  pointer  :: ZANDAYVT, ZANFMVT
   REAL(KIND=JPRB), dimension(:,:,:),pointer  :: ZDHVEGS
end type surf_and_more_local_type

type aux_diag_local_type
  INTEGER(KIND=JPIM), pointer :: IEXT3D   ! position in extra field for diagnostics 
  REAL(KIND=JPRB), dimension(:), pointer :: ZWND     ! horizontal wind in the lowest model level
  REAL(KIND=JPRB), dimension(:), pointer :: ZCCNL, ZCCNO
  INTEGER(KIND=JPIM), dimension(:), pointer :: ITOPC, IBASC, IBOTSC
  ! radiation scheme for convenience
  REAL(KIND=JPRB), dimension(:), pointer ::  ZEMIT ! surface longwave char
  REAL(KIND=JPRB), dimension(:), pointer ::  ZTRSOD, ZSUDU, ZDSRP !  radiation working arrays
  ! Skin temperature from radiation
  REAL(KIND=JPRB), dimension(:), pointer :: ZTSKRAD   !  Skin temperature from radiation
  REAL(KIND=JPRB), dimension(:,:), pointer ::  ZTHT ! half level temperature
  REAL(KIND=JPRB), dimension(:,:), pointer ::  ZCTRSO,ZCEMTR ! radiation working arrays
  !     LOCAL 3D ARRAYS TO PASS OUTPUT FROM GWD TO IMPLICIT SOLVER OF VDF
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZSOTEU    ! Explicit part of U-tendency from subgrid orography scheme
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZSOTEV    ! Explicit part of V-tendency from subgrid orography scheme
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZSOBETA   ! Implicit part of subgrid orography
  ! aerosols in microphysics
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZLCRIT_AER ! critical liquid mmr for rain autoconversion process 
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZICRIT_AER ! critical liquid mmr for snow autoconversion process 
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZRE_LIQ    ! effective radius liquid
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZRE_ICE    ! effective radius ice
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZCCN       ! CCN (prognostic, diagnostic)
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZNICE      ! ice no concentration
  !    Local arrays for new cloud scheme of the linearized model
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZPATMP     ! CLOUD COVER OF INDIVIDUAL LAYER
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZPLSMP     ! ICE WATER CONTENT
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZPLTMP     ! LIQUID WATER CONTENT
  ! liquid water loading in environment neglected in convection
  REAL(KIND=JPRB), dimension(:,:), pointer :: ZLISUM     ! Local array for sum of liquid water and ice to be available
  ! local arrays variables for lightning
  ! REAL(KIND=JPRB)    :: ZLIGH_TOT(KDIM%KLON), ZLIGH_CTG(KDIM%KLON), ZLIGH_EMI(KDIM%KLON), ZWMFU(KDIM%KLON)
  ! REAL(KIND=JPRB)    :: ZCTOPH(KDIM%KLON), ZPRECMX(KDIM%KLON), ZICETOT(KDIM%KLON), ZCDEPTH(KDIM%KLON)
end type aux_diag_local_type

type keys_local_type
  LOGICAL, dimension(:), pointer  :: LLLAND, LLSICE, LLNH, LLCUM, LLSC, LLSHCV
  LOGICAL, dimension(:), pointer  :: LLLAKE, LLOCN_KPP
  LOGICAL, pointer ::  LLCLDCOVER     !  key deciding to compute/not compute the cloud cover
  LOGICAL, pointer ::  LLMAINCALL     !  key controlling the complexity of cloudsc scheme
end type keys_local_type

end module yomphyder

