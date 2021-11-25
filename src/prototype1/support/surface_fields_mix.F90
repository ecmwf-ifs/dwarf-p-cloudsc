! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SURFACE_FIELDS_MIX

!     Purpose.
!     --------

!      SURFACE_FIELDS_MIX contains data structures and manipulation routines
!      for the surface (physics) fields in the IFS            

!      This module is a mix of declarations, type definitions and 
!       subroutines linked with surface fields. There are four parts:
!       1/ Declaration of dimensions (including some parameter variables).
!       2/ Definition of types.
!       3/ Declarations:
!          Declaration of variables SP_[group], YSP_[group]D, YSP_[group]
!           (prognostic surface fields).
!          Declaration of variables SD_[group], YSD_[group]D, YSD_[group]
!           (diagnostic surface fields).
!       4/ Some routines linked to the surface data flow:
!          * INI_SFLP3: Initialize 3-D surface field group
!          * SETUP_SFLP3: Setup 3-D surface field
!          * FINALISE_SFLP3: Finalise 3-D surface field
!          * INI_SFLP2: Initialize 2-D surface field group
!          * SETUP_SFLP2: Setup 2-D surface field
!          * FINALISE_SFLP2: Finalise 2-D surface field
!          * GPPOPER: Operations on prognostic surface fields
!          * GPOPER: Operations on ALL surface groups
!          * GPOPER_2: Operations on 2-D surface groups
!          * GPOPER_3: Operations on 3-D surface groups
!          * SURF_STORE: Store all surface fields
!          * SURF_RESTORE: Restore all surface fields
!          * ALLO_SURF: Allocate surface field arrays

!     Author.
!     -------
!     Mats Hamrud(ECMWF)

!     Modifications.
!     --------------
!        Original : 2006-07-01
!        Modifications:
!        K. Yessad (25 Oct 2006): rephase ALARO0 contribution.
!        K. Yessad (26 Oct 2006): add missing comments.
!        G. Balsamo (14 Mar 2007): add soil type pointer (SOTY).
!        Jean Bidlot (June 2007):  named pointer for fields to wave model.
!        S. Serrar (17 Jul 2007) methane surface fields pointers
!        Y. Takaya (?? Nov 2008): add ocean mixed layer model fields
!        A. Alias  (07 Aug 2009) move field Nudging mask from group VCLIA to VARSF
!        JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface 
!        H. Hersbach (04-Dec-2009): 10m-neutral wind and friction velocity,
!                                   introduce YDUPD, SETPERTOVAL in GPOPEN
!        S. Boussetta/G.Balsamo (May 2009): add high/low vegetation LAI pointer (LAIH/LAIL)
!        Y. Takaya/P. de Rosnay May 2010: SSS for SMOS
!        R. Forbes (01-Mar-2010): Added TCRW,TCSW diagnostics
!        JJMorcrette 20100212 PP of CBASE, 0DEGL and VISIH
!        T. Wilhelmsson (25 Mar 2010) add 6 hourly min/max fields
!        A. Alias  (07 Aug 2009) move field Nudging mask from group VCLIA to VARSF
!        Y. Bouteloup (20 Oct 2010) Surface forcing for 1D model : SFORC
!        Y. Bouteloup (04 Jan 2011) Surface flux for EDKF and 1D model surface forcing : SFLUX
!        G.Balsamo/S.Boussetta (Apr 2011): add land carbon dioxide fluxes
!        H. Hersbach (01 April 2011): auxiliary diagnostic radiation fields
!        P. Marguinaud (07 November 2012): Fix uninitialized component
!        P. Bechtold (9 Aug 2011): add CIN, Convective Indices
!        P.Marguinaud (11 Sep 2012): Initialize TYPE_SURF_MTL_2D%CNAME (avoid valgrind warning)
!        M. Ahlgrimm 31 Oct 2011: clear-sky downward radiation at surface
!        L. Jones (25-Oct-2011): Created FINALISE_SFLP2/3, removing need for 
!                                pre-counting number of SETUP_SFLP2/3 calls
!        M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!        JJMorcrette 20130213 PP optical depths GEMS/MACC aerosols
!        G. Balsamo  14-Jun-2013 Introduce lake buffer group SL=LAKEB
!        R. Forbes 01-March-2014 Add precip rates/type,TCSLW,I10FG,PEV
!        M. Ahlgrimm Apr 2014: Add precip fraction for DDH output
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        JJMorcrette 20130730 15-variable aerosol model + oceanic DMS
!        R. Forbes 10-Jan-2015 Add freezing rain FZRA
!        S. Rémy, 23/1/2015, injection height for biomass burning aerosol emissions
!        A. Agusti-Panareda (31 Oct 2013): add GPP/REC flux adjustment coefficient pointer (CGPP,CREC)
!-------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
!USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDIM   , ONLY : YRDIM
USE YOMLUN   , ONLY : NULOUT, NULERR
USE YOMCT0   , ONLY : LTWOTL, NUNDEFLD
USE YOMDYN   , ONLY : YRDYN
IMPLICIT NONE
SAVE

#include "abor1.intfb.h"
!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPMAXSFLDS=610 ! Max number of fields in individual group  !KPP
INTEGER(KIND=JPIM), PARAMETER :: JPMAXSTRAJ=100 ! Dimension of NSTRAJGRIB
INTEGER(KIND=JPIM) :: NOFFTRAJ                  ! Offset in surf trajectory
INTEGER(KIND=JPIM) :: NOFFTRAJ_CST              ! Offset in "constant" surf trajectory
INTEGER(KIND=JPIM) :: NSTRAJGRIB(JPMAXSTRAJ)    ! Used in trajectory setup
INTEGER(KIND=JPIM), PRIVATE :: NPTRSURF         ! Used by routine GPOPER

! General type defintions

! 2D surface field structure
TYPE TYPE_SURF_MTL_2D
  INTEGER(KIND=JPIM) :: MP           ! Basic field pointer
  INTEGER(KIND=JPIM) :: MP0          ! Field pointer timelevel  0 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP9          ! Field pointer timelevel -1 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP1          ! Field pointer timelevel +1 (prognostic fields)
  INTEGER(KIND=JPIM) :: IGRBCODE     ! GRIB parameter code (default: -999)
  CHARACTER(LEN=16)  :: CNAME = '                '
                                     ! ARPEGE field name   (default: all spaces)

  REAL(KIND=JPRB)    :: REFVALI      ! Default value       (default: 0.0)
  INTEGER(KIND=JPIM) :: NREQIN       ! -1 - initial value from default (default)
                                     ! +1 - initial value from reading file
                                     !  0 - no initial value
  INTEGER(KIND=JPIM) :: ITRAJ        !  0 not in trajectory (default)
                                     !  1 in trajectory
                                     !  2 in "constant" trajectory
  INTEGER(KIND=JPIM) :: IBUGFINDER   ! Debug use only. Can be removed if confident.
  LOGICAL            :: LSET         ! True if structure has been set up
END TYPE TYPE_SURF_MTL_2D

! 3D surface field structure
TYPE TYPE_SURF_MTL_3D
  INTEGER(KIND=JPIM) :: MP   ! Basic field pointer
  INTEGER(KIND=JPIM) :: MP0  ! Field pointer timelevel  0 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP9  ! Field pointer timelevel -1 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP1  ! Field pointer timelevel +1 (prognostic fields)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IGRBCODE(:)  ! GRIB parameter code (default: -999)
  CHARACTER(LEN=16) , ALLOCATABLE :: CNAME(:)     ! ARPEGE field name   (default: all spaces)
  REAL(KIND=JPRB)   , ALLOCATABLE :: REFVALI(:)   ! Default value       (default: 0.0)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NREQIN(:)    ! -1 - initial value from default (default)
                                                  ! +1 - initial value from reading file
                                                  !  0 - no initial value
  INTEGER(KIND=JPIM) :: ITRAJ                     !  0 not in trajectory (default)
                                                  !  1 in trajectory
                                                  !  2 in "constant" trajectory
  INTEGER(KIND=JPIM) :: IBUGFINDER                ! Debug use only. Can be removed if confident.
  LOGICAL            :: LSET                      ! True if structure has been set up
END TYPE TYPE_SURF_MTL_3D

! Descriptor pertaining to group
TYPE TYPE_SURF_GEN
  INTEGER(KIND=JPIM) :: NUMFLDS         ! Number of field in group
  INTEGER(KIND=JPIM) :: NDIM            ! Field dimension
  INTEGER(KIND=JPIM) :: NLEVS           ! Number of levels (for multi level groups)
  INTEGER(KIND=JPIM) :: IPTR            ! Internal use
  INTEGER(KIND=JPIM) :: IPTR5           ! Internal use
  INTEGER(KIND=JPIM) :: NDIM5           ! Dimension of trajectory array
  INTEGER(KIND=JPIM) :: NOFFTRAJ        ! Internal use
  INTEGER(KIND=JPIM) :: NOFFTRAJ_CST    ! Internal use
  CHARACTER(LEN=16)  :: CGRPNAME        ! Name of group (for prints)
  LOGICAL            :: L3D             ! TRUE if multi-level field (3-D)
  LOGICAL            :: LMTL            ! TRUE if prognostic field (multi time level)
  LOGICAL            :: FINALISED       ! TRUE if group finalised & no more fields allowed
END TYPE TYPE_SURF_GEN

! Type descriptor for derived type for communicating with GPOPER (see below)
TYPE TYPE_SFL_COMM
  INTEGER(KIND=JPIM) :: IGRBCODE
  LOGICAL            :: L_OK 
  CHARACTER(LEN=16)  :: CNAME
  INTEGER(KIND=JPIM) :: IFLDNUM
  REAL(KIND=JPRB)    :: VALUE
  INTEGER(KIND=JPIM) :: IPTRSURF
  INTEGER(KIND=JPIM) :: ICODES(JPMAXSFLDS)
  INTEGER(KIND=JPIM) :: ICOUNT
END TYPE TYPE_SFL_COMM

! Group specific type definitions: pronostic groups (SB, SG, RR, CL, OM, EP, X2, CI)

! * Group SB=SOILB: soil prognostic quantities for the different reservoirs
!    (four reservoirs at ECMWF, deep reservoir at METEO-FRANCE):
TYPE TYPE_SFL_SOILB
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YT    ! temperature
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YQ    ! liquid water content
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YTL   ! ice water content (for MF)
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YSB(:) => NULL()
END TYPE TYPE_SFL_SOILB

! * Group SG=SNOWG: surface snow prognostic quantities:
TYPE TYPE_SFL_SNOWG
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YF    ! content of surface snow
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YA    ! snow albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YR    ! snow density
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YT    ! total albedo (diagnostic for MF for LVGSN)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSG(:) => NULL()
END TYPE TYPE_SFL_SNOWG

! * Group SL=LAKEB: Lake (FLAKE Model) prognostic quantities:
TYPE TYPE_SFL_LAKEB
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLICT  ! lake ice temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLMLT  ! lake mixed-layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLTLT  ! lake total layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLBLT  ! lake bottom layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSHF  ! lake shape factor
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLICD  ! lake ice depth
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLMLD  ! lake mixed-layer depth
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YSL(:) => NULL()
END TYPE TYPE_SFL_LAKEB

! * Group RR=RESVR: surface prognostic quantities (ECMWF) or
!   surface + superficial reservoir prognostic quantities (MF):
!   Remark:
!    at ECMWF there are 4 soil reservoirs and there is a
!    clear distinction between the soil reservoirs (group SOILB)
!    and the surface (group RESVR);
!    at METEO-FRANCE there is a deep reservoir (group SOILB) and a
!    superficial reservoir (group RESVR):
!    - there is a skin surface temperature (Ts) which is the temperature at the
!      interface surface/superficial reservoir (and not two separate quantities
!      for superficial reservoir and surface)
!    - there is a skin surface water content (denoted by Wl) and a superficial
!      reservoir water content (denoted by Ws).
!    - there is a superficial reservoir ice content but no surface ice content.
!    (remark k.y.: it would have been more logical to use group name
!    RESVR for internal reservoirs and group name SOILB for surface!).
TYPE TYPE_SFL_RESVR
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YT    ! skin temperature (Ts)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YW    ! skin water content (Wskin) at ECMWF
                                          ! superficial reservoir water content (Ws) at MF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFC   ! skin water content (Wl) at MF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIC   ! superficial reservoir ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFP1  ! interpolated Ts for 2nd part of 927-FULLPOS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRR(:) => NULL()
END TYPE TYPE_SFL_RESVR

! * Group CL=CLS: surface boundary layer prognostic quantities:
TYPE TYPE_SFL_CLS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCLS    ! 2m temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHUCLS   ! 2m humidity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCL(:) => NULL()
END TYPE TYPE_SFL_CLS

! * Group OM=OML: prognostic quantities for ocean mixed layer model (KPP/TKE):
TYPE TYPE_SFL_OML
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTO     ! temperature
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YSO     ! salinity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YUO     ! U velocity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YVO     ! V velocity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YOM(:) => NULL()
END TYPE TYPE_SFL_OML 

! * Group EP=EXTRP: extra 3-d prognostic fields:
TYPE TYPE_SFL_EXTRP
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YEP(:) => NULL()
END TYPE TYPE_SFL_EXTRP

! * Group X2=XTRP2: extra 2-d prognostic fields:
!   (is used for precipitation fields in CANARI)
TYPE TYPE_SFL_XTRP2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YX2(:) => NULL()
END TYPE TYPE_SFL_XTRP2

! * Group CI=CANRI: 2-d prognostic fields for CANARI:
TYPE TYPE_SFL_CANRI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCI(:) => NULL()
END TYPE TYPE_SFL_CANRI

! Group specific type definitions: diagnostic groups 

! * Group VF=VARSF: climatological/geographical diagnostic fields:
TYPE TYPE_SFL_VARSF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0F    ! gravity * surface roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALBF   ! surface shortwave albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEMISF  ! surface longwave emissivity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGETRL  ! standard deviation of orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSM    ! land-sea mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVEG    ! vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVRLAN  ! anisotropy of the sub-grid scale orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVRLDI  ! angle of the direction of orography with the x axis
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSIG    ! characteristic orographic slope
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALBSF  ! soil shortwave albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAN    ! fraction of land
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSST    ! (open) sea surface temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSS    ! sea surface salinity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLZ0H   ! logarithm of roughness length for heat
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCVL    ! low vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCVH    ! high vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTVL    ! low vegetation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTVH    ! high vegetation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAIL   ! low vegetation LAI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAIH   ! high vegetation LAI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOTY   ! soil type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCLK    ! lake cover 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDL     ! lake depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCI     ! sea ice fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUCUR   ! U-component of the ocean current
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVCUR   ! V-component of the ocean current
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0RLF  ! gravity * vegetation roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCO2O   ! oceanic CO2 flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCO2B   ! biosphere CO2 flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCO2A   ! anthropogenic CO2 flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCO2F   ! CO2 fire emissions
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCGPP   ! GPP bias correction factor
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCREC   ! REC bias correction factor
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCH4AG  ! CH4 surface fluxes - aggregated field
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCH4F   ! CH4 fire emissions
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSDFOR  ! SD filtered orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVP  ! MODIS-derived parallel albedo for shortwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVD  ! MODIS-derived diffuse albedo for shortwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNIP  ! MODIS-derived parallel albedo for longwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNID  ! MODIS-derived diffuse albedo for longwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFP1    ! surface orography in the 2nd part of FULLPOS-927
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBCBF   ! black carbon biogenic
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBCFF   ! black carbon fossil fuel
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBCGF   ! black carbon GFED
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YOMBF   ! organic matter biogenic
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YOMFF   ! organic matter fossil fuel
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YOMGF   ! organic matter GFED
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YINJF   ! height of maximum injection for biomass burning emissions
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSO2L   ! sulphate low-level 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSO2H   ! sulphate higher-level 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOGF   ! sulphate GFED
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOA    ! secondary organic
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVOLC   ! volcanic continuous
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVOLE   ! volcanic explosive
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDMSO   ! oceanic DMS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERDEP ! dust emission potential 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERLTS ! dust lifting threshold speed 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERSCC ! dust soil clay content
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTRACFLX(:) ! TRAC emissions 1:NTRAC  
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMFLX(:) ! chemistry emissions
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMDV(:) ! chemistry deposition velocity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNUDM   ! nudging mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVF(:) => NULL()
END TYPE TYPE_SFL_VARSF

! * Group VP=VCLIP: deep soil diagnostic fields:
TYPE TYPE_SFL_VCLIP
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPC    ! climatological deep layer temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWPC    ! climatological deep layer moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVP(:) => NULL()
END TYPE TYPE_SFL_VCLIP

! * Group VV=VCLIV: vegetation diagnostic fields:
TYPE TYPE_SFL_VCLIV
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YARG    ! silt percentage within soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSAB    ! percentage of sand within the soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YD2     ! soil depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIVEG   ! type of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRSMIN  ! stomatal minimum resistance
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAI    ! leaf area index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHV     ! resistance to evapotranspiration
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0H    ! gravity * roughness length for heat
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALS    ! albedo of bare ground
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALV    ! albedo of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVV(:) => NULL()
END TYPE TYPE_SFL_VCLIV

! * Group VN=VCLIN: cloudiness diagnostic predictors:
TYPE TYPE_SFL_VCLIN
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTOP    ! index of convective cloud top
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBAS    ! index of convective cloud base
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YACPR   ! averaged convective precipitaion rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YACCPR  ! accumulated total precipitaion for assimilation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVN(:) => NULL()
END TYPE TYPE_SFL_VCLIN

! * Group VH=VCLIH: convective cloud diagnostic fields:
TYPE TYPE_SFL_VCLIH
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCCH  => NULL() ! total convective cloudiness
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSCCH  => NULL() ! convective cloud summit
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBCCH  => NULL() ! convective cloud base
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPBLH  => NULL() ! PBL height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPSH  => NULL() ! variable for prognostic convection scheme (ALARO)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YQSH   => NULL() ! surface moisture historic variable (used by TOUCANS)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVH(:) => NULL()
END TYPE TYPE_SFL_VCLIH

! * Group VA=VCLIA: aerosol diagnostic fields:
TYPE TYPE_SFL_VCLIA
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSEA    ! aerosol: sea
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAN    ! aerosol: land
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOO    ! aerosol: soot
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDES    ! aerosol: desert
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUL    ! aerosol: sulfate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVOL    ! aerosol: volcano
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVA(:) => NULL()
END TYPE TYPE_SFL_VCLIA

! * Group VG=VCLIG: ice-coupler diagnostic fields:
!   ky: currently not used, missing setup in su_surf_flds.F90
TYPE TYPE_SFL_VCLIG
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YICFR   ! sea-ice fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOUP   ! upward solar flux over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIRUP   ! upward IR flux over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHSS   ! sensible heat over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEVAP   ! evaporation over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUX   ! U-component of stress over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUY   ! V-component of stress over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVG(:) => NULL()
END TYPE TYPE_SFL_VCLIG

! * Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields:
TYPE TYPE_SFL_VO3ABC
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YA      ! A climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YB      ! B climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YC      ! C climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVC(:) => NULL()
END TYPE TYPE_SFL_VO3ABC

! * Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
TYPE TYPE_SFL_VDIAGO2
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YOCDEP  ! bottom layer depth
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YUSTRC  ! taux clim.
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YVSTRC  ! tauy clim.
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YV2(:) => NULL()
END TYPE TYPE_SFL_VDIAGO2 

! * Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
TYPE TYPE_SFL_VDIAGO3 
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFM   ! viscosity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFT   ! diff. coef. of temp
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFS   ! diff. coef. of salinity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YADVT   ! correction term for temp.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YADVS   ! correction term for sal.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTRI0   ! coef. for solving matrix.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTRI1   ! coef. for solving matrix.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YSWDK   ! radiation term
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YZO     ! depth of layer
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YHO     ! depth of interface layer
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDO     ! layer thickness
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YHO_INV ! 1 / YHO
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YUOC    ! U velocity clim.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YVOC    ! V velocity clim.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YV3(:) => NULL()
END TYPE TYPE_SFL_VDIAGO3 

! * Group VD=VDIAG: (ECMWF) diagnostic fields:
TYPE TYPE_SFL_VDIAG
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSP       ! Large scale precipitation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCP        ! Convective precipitation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSF        ! Snowfall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFZRA      ! Freezing rain
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBLD       ! Boundary layer dissipation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSHF      ! Surface sensible heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSLHF      ! Surface latent heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNEE       ! Surface net ecosystem exchange of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGPP       ! Surface gross primary production of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YREC       ! Surface ecosystem respiration of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMSL       ! Mean sea level pressure
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSP        ! Surface pressure
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCC       ! Total cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10U       ! U-wind at 10 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10V       ! V-wind at 10 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y2T        ! Temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y2D        ! Dewpoint temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSR       ! Surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTR       ! Surface thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSR       ! Top solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTR       ! Top thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEWSS      ! Instantaneous surface U-wind stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNSSS      ! Instantaneous surface V-wind stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YE         ! Water evaporation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPEV       ! Potential evaporation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCCC       ! Convective cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLCC       ! Low cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMCC       ! Medium cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHCC       ! High cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLGWS      ! Zonal gravity wave stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMGWS      ! Meridian gravity wave stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGWD       ! Gravity wave dissipation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMX2T      ! Maximum temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMN2T      ! Minimum temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMX2T6(:)  ! Bins for maximum temperature at 2 m since last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMN2T6(:)  ! Bins for minimum temperature at 2 m since last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRO        ! Runoff (total)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSRO       ! Runoff surface
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRO      ! Runoff sub-surface
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALB       ! (surface shortwave) albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIEWSS     ! Instantaneous surface zonal component of stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YINSSS     ! Instantaneous surface meridian component of stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YISSHF     ! Instantaneous surface heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIE        ! Instantaneous surface moisture flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YINEE      ! Instantaneous net ecosystem exchange of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIGPP      ! Instantaneous gross primary production of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIREC      ! Instantaneous ecosystem respiration of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCSF       ! Convective snow fall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSSF      ! Large scale snowfall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXTPR     ! Max precip rate since last post-processing
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMNTPR     ! Min precip rate since last post-processing
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXTPR6(:) ! Max precip rate in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMNTPR6(:) ! Min precip rate in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSRR      ! Large scale rain rate 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCRR       ! Convective rain rate 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSSFR     ! Large scale snowfall rate 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCSFR      ! Convective snowfall rate 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPTYPE     ! Precipitation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YILSPF     ! Large-scale precipitation fraction (inst.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0F       ! Gravity * surface roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLZ0H      ! Logarithm of z0 times heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCW       ! Total water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCWV      ! Total water vapor content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCLW      ! Total liquid water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCIW      ! Total ice water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCRW      ! Total rain water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCSW      ! Total snow water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCSLW     ! Total supercooled liquid water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRD      ! Downward surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRD      ! Downward surface thermic radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRDC     ! Clear-sky downward surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRDC     ! Claer-sky downward surface thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBLH       ! Height of boundary layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUND      ! Sunshine duration
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPAR      ! Surface downward PARadiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUVB      ! Surface downward UV-B radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFDIR     ! Surface total sky direct downward SW radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSCDIR     ! Surface clear-sky direct downward SW radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCAPE      ! Conv.avail.potential energy (CAPE)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSRC      ! Top solar radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTRC      ! Top thermal radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRC      ! Surface solar radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRC      ! Surface thermal radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YES        ! Evaporation of snow
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSMLT      ! Snow melt
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10FG      ! Wind gust at 10 m (max since previous pp)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10FG6(:)  ! Bins for wind gust at 10 m (max since last 6 hours)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YI10FG     ! Wind gust at 10 m ("instantaneous")
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSPF      ! Large scale precipitation fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCO3      ! Total ozone content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVIMD      ! Vertically integrated mass divergence
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPARC     ! Surface clear-sky parallel radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTINC     ! Top of atmosphere incident solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCBASE     ! Cloud base level
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y0DEGL     ! Zero deg. level
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVISIH     ! Horizontal visibility
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCIN       ! CIN
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YKINDEX    ! Convective K-Index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTINDEX   ! Convective TT-Index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCGHG(:)  ! Total column greenhouse gases
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCTRAC(:) ! Total column tracers
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCCHEM(:) ! Total column chemistry
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y100U      ! 100m zonal wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y100V      ! 100m meridional wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZUST      ! Friction velocity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10NU      ! 10m zonal neutral wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10NV      ! 10m meridional neutral wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDNDZN     ! Minimum vertical refractivity gradient
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDNDZA     ! Mean vertical refractivity gradient
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDCTB      ! Duct base height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPLB      ! Trapping layer base height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPLT      ! Trapping layer top height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODTO      ! optical depth total aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODSS      ! optical depth sea salt aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODDU      ! optical depth dust aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODOM      ! optical depth organic m. aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODBC      ! optical depth black C aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODSU      ! optical depth sulphate aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODVFA     ! optical depth volcanic flying ash    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODVSU     ! optical depth volcanic sulphate aerosols    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM1     ! particulate matter le 1 um    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM25    ! particulate matter le 2.5um    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM10    ! particulate matter le 10 um    
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVD(:) => NULL()
END TYPE TYPE_SFL_VDIAG

! * Group WS=WAVES: surface prognostic quantities over sea (used by IFS):
TYPE TYPE_SFL_WAVES
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWS(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHAR      ! Charnock parameter as modified by the wave model.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUSTOKES   ! U-component of the surface Stokes drift.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVSTOKES   ! V-component of the surface Stokes drift.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPHIOC     ! Energy flux to ocean.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPHIAW     ! Energy flux to ocean waves.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUOC     ! Momentum flux to ocean.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEMEAN     ! Wave variance.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFMEAN     ! Wave mean frequency.
END TYPE TYPE_SFL_WAVES

! * Group WW=WAM: surface prognostic quantities over sea (used by WAM):
TYPE TYPE_SFL_WAM
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWW(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YU10N      ! 10m neutral wind U-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YV10N      ! 10m neutral wind V-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRHO       ! surface density passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZIL       ! ZI/L passed to the wave model (used for gustiness in WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCIF       ! Sea ice fraction passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUCURW     ! Ocean current    U-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVCURW     ! Ocean current    V-component passed to the wave model (WAM).
END TYPE TYPE_SFL_WAM

! * Group VX=VCLIX: auxilary climatological diagnostic fields:
TYPE TYPE_SFL_VCLIX
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YORO    ! climatological surface geopotential
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSC    ! climatological surface temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPWS    ! climatological surface max. prop. moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPWP    ! climatological deep soil max. prop. moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSNO    ! climatological snow cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPC    ! climatological deep soil temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSAB    ! climatologic percentage of sand within the soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YXD2    ! climatologic soil depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSM    ! climatologic land sea mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIVEG   ! climatologic type of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVX(:) => NULL()
END TYPE TYPE_SFL_VCLIX

! * Group XA=VEXTRA: extra 3-d diagnostic fields:
TYPE TYPE_SFL_VEXTRA
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXA(:) => NULL()
END TYPE TYPE_SFL_VEXTRA

! * Group X2=VEXTR2: extra 2-d diagnostic fields:
TYPE TYPE_SFL_VEXTR2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YX2(:) => NULL()
END TYPE TYPE_SFL_VEXTR2

! * Group SFL:SFLUX Surface flux for EDKF 
TYPE TYPE_SFL_SFLUX
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFL(:) => NULL()
END TYPE TYPE_SFL_SFLUX

! * Group SFO:SFORC Surface forcing for 1D model (MUSC)
TYPE TYPE_SFL_SFORC
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFO(:) => NULL()
END TYPE TYPE_SFL_SFORC

! End of type definitions

TYPE :: TSURF

INTEGER(KIND=JPIM) :: NSURF=0                   ! Number of surf var.
INTEGER(KIND=JPIM) :: NSURFL=0                  ! Number of surf flds (fields*levels)
INTEGER(KIND=JPIM) :: NDIMSURF=0                ! Total of surf var (includes timelevels etc)
INTEGER(KIND=JPIM) :: NDIMSURFL=0               ! Total dimension of all surface variables
INTEGER(KIND=JPIM) :: NPROGSURF=0               ! Number of prognostic surf var.
INTEGER(KIND=JPIM) :: NPROGSURFL=0              ! Number of prognostic surf flds (fields*levels)

REAL(KIND=JPRB), ALLOCATABLE :: STORE_ARRAY(:,:,:) ! Backup array for surf (see routine SURF_STORE )

! Information on status of tl perturbations
TYPE(TYPE_SFL_COMM) :: YDUPD

! Data structures

! Prognostic (multi time level) fields (SB, SG, SL, RR, CL, OM, EP, X2, CI)

! Soilb
REAL(KIND=JPRB), ALLOCATABLE :: SP_SB (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_SBD
TYPE(TYPE_SFL_SOILB)         :: YSP_SB

! Snowg
REAL(KIND=JPRB), ALLOCATABLE :: SP_SG (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_SGD
TYPE(TYPE_SFL_SNOWG)         :: YSP_SG

! Lakeb
REAL(KIND=JPRB),ALLOCATABLE :: SP_SL (:,:,:) 
TYPE(TYPE_SURF_GEN)     :: YSP_SLD
TYPE(TYPE_SFL_LAKEB)    :: YSP_SL

! Resvr
REAL(KIND=JPRB), ALLOCATABLE :: SP_RR (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_RRD
TYPE(TYPE_SFL_RESVR)         :: YSP_RR

! Cls
REAL(KIND=JPRB), ALLOCATABLE :: SP_CL (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_CLD
TYPE(TYPE_SFL_CLS)           :: YSP_CL

! Oml (used by ocean mixed layer model (KPP))
REAL(KIND=JPRB), ALLOCATABLE :: SP_OM (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_OMD 
TYPE(TYPE_SFL_OML)           :: YSP_OM

! Extrp
REAL(KIND=JPRB), ALLOCATABLE :: SP_EP (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_EPD
TYPE(TYPE_SFL_EXTRP)         :: YSP_EP 

! Xtrp2
REAL(KIND=JPRB), ALLOCATABLE :: SP_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_X2D
TYPE(TYPE_SFL_XTRP2)         :: YSP_X2

! Canri
REAL(KIND=JPRB), ALLOCATABLE :: SP_CI (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_CID
TYPE(TYPE_SFL_CANRI)         :: YSP_CI

! Diagnostic (one time level) fields

! Varsf
REAL(KIND=JPRB), ALLOCATABLE :: SD_VF (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VFD
TYPE(TYPE_SFL_VARSF)         :: YSD_VF

! Vclip
REAL(KIND=JPRB), ALLOCATABLE :: SD_VP (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VPD
TYPE(TYPE_SFL_VCLIP)         :: YSD_VP

! Vcliv
REAL(KIND=JPRB), ALLOCATABLE :: SD_VV (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VVD
TYPE(TYPE_SFL_VCLIV)         :: YSD_VV

! Vclin
REAL(KIND=JPRB), ALLOCATABLE :: SD_VN (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VND
TYPE(TYPE_SFL_VCLIN)         :: YSD_VN

! Vclih
REAL(KIND=JPRB), ALLOCATABLE :: SD_VH (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VHD
TYPE(TYPE_SFL_VCLIH)         :: YSD_VH

! Vclia
REAL(KIND=JPRB), ALLOCATABLE :: SD_VA (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VAD
TYPE(TYPE_SFL_VCLIA)         :: YSD_VA

! Vclig
! currently nothing declared

! Vo3abc
REAL(KIND=JPRB), ALLOCATABLE :: SD_VC (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VCD
TYPE(TYPE_SFL_VO3ABC)        :: YSD_VC

! Vdiago2 (used by ocean mixed layer model (KPP)) 
REAL(KIND=JPRB), ALLOCATABLE :: SD_V2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_V2D 
TYPE(TYPE_SFL_VDIAGO2)       :: YSD_V2 

! Vdiago3 (used by ocean mixed layer model (KPP))
REAL(KIND=JPRB), ALLOCATABLE :: SD_V3 (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_V3D 
TYPE(TYPE_SFL_VDIAGO3)       :: YSD_V3 

! Vdiag
REAL(KIND=JPRB), ALLOCATABLE :: SD_VD (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VDD
TYPE(TYPE_SFL_VDIAG)         :: YSD_VD

! Waves (used by IFS)
REAL(KIND=JPRB), ALLOCATABLE :: SD_WS (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_WSD
TYPE(TYPE_SFL_WAVES)         :: YSD_WS

! Waves (used by WAM)
REAL(KIND=JPRB), ALLOCATABLE :: SD_WW (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_WWD
TYPE(TYPE_SFL_WAM)           :: YSD_WW

! Vclix
REAL(KIND=JPRB), ALLOCATABLE :: SD_VX (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VXD
TYPE(TYPE_SFL_VCLIX)         :: YSD_VX

! Vextra
REAL(KIND=JPRB), ALLOCATABLE :: SD_XA (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XAD
TYPE(TYPE_SFL_VEXTRA)        :: YSD_XA

! Vextra-radiation
REAL(KIND=JPRB), ALLOCATABLE :: SD_XR (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XRD
TYPE(TYPE_SFL_VEXTRA)        :: YSD_XR

! Vextr2
REAL(KIND=JPRB), ALLOCATABLE :: SD_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_X2D
TYPE(TYPE_SFL_VEXTR2)        :: YSD_X2

! SFLUX
REAL(KIND=JPRB), ALLOCATABLE :: SD_SFL (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_SFLD
TYPE(TYPE_SFL_SFLUX)         :: YSD_SFL

! SFORC
REAL(KIND=JPRB), ALLOCATABLE :: SD_SFO (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_SFOD
TYPE(TYPE_SFL_SFORC)         :: YSD_SFO

! Precip fraction
REAL(KIND=JPRB),ALLOCATABLE :: SD_PF (:,:,:,:)
TYPE(TYPE_SURF_GEN)     :: YSD_PFD
TYPE(TYPE_SFL_VEXTRA)   :: YSD_PF

END TYPE TSURF

!-------------------------------------------------------------------------

TYPE(TSURF), POINTER :: YRSURF => NULL()

CONTAINS

!=========================================================================

SUBROUTINE INI_SFLP3(YDSC,YD,KLEVS,LDMTL,CDGRPNAME)
! Initialize 3-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT) :: YD(:)
INTEGER(KIND=JPIM),INTENT(IN)        :: KLEVS
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP3',0,ZHOOK_HANDLE)

YDSC%FINALISED = .FALSE.
YDSC%NLEVS = KLEVS
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST

DO JFLD=1,SIZE(YD)
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
  YD(JFLD)%LSET = .FALSE.
  YD(JFLD)%IBUGFINDER = JFLD  ! For debug use only
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 3-D SURFACE FIELD GROUP ', YDSC%CGRPNAME

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP3

!=========================================================================

SUBROUTINE SETUP_SFLP3(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 3-D surface field
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB(:)
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN(:)

INTEGER(KIND=JPIM) :: IPTR,JLEV,J
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP3',0,ZHOOK_HANDLE)

IPTR = YDSC%IPTR
IF(IPTR > JPMAXSFLDS-1) THEN
  ! We are about to set the last element which is generally used as a
  ! "missing field" placeholder.
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSIONED - GROUP ',&
   & YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('IPTR > JPMAXSFLDS-1')
ENDIF

IF (YDSC%FINALISED) THEN
  WRITE(NULERR,*) 'ERROR - SETUP_SFLP3 CALLED FOR A GROUP WHICH HAS BEEN FINALISED',&
 &  YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('SETUP_SFLP3 CALLED FOR A GROUP WHICH HAS BEEN FINALISED')
ENDIF

! If the YD(:) elements are not being set in order it's a possible sign of a 
! bug - the caller may have passed in a pointer to the wrong element. Or it
! might be entirely intentional and an abort is unnecessary.
! This if-block and all other references to IBUGFINDER can be commented if
! they cause a problem.
IF (YD%IBUGFINDER /= IPTR) THEN
  WRITE(NULERR,*) 'SETUP_SFLP3 NOT BEING CALLED FOR YD ELEMENTS IN ORDER. '//&
    & '-POSSIBLY- A BUG.',YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('SETUP_SFLP3 NOT CALLED FOR YD ELEMENTS IN ORDER. POSSIBLE BUG.')
ENDIF

ALLOCATE(YD%IGRBCODE(YDSC%NLEVS))
ALLOCATE(YD%CNAME(YDSC%NLEVS))
ALLOCATE(YD%REFVALI(YDSC%NLEVS))
ALLOCATE(YD%NREQIN(YDSC%NLEVS))

IF(PRESENT(KGRIB)) THEN
  DO J=1,SIZE(YD%IGRBCODE)
    YD%IGRBCODE(J) = KGRIB(J)
  ENDDO
ELSE
  YD%IGRBCODE(:) = -999
ENDIF

IF(PRESENT(KREQIN)) THEN
  DO J=1,SIZE(YD%NREQIN)
    YD%NREQIN(J) = KREQIN(J)
  ENDDO
ELSE
  YD%NREQIN(:) = -1
ENDIF 

IF(PRESENT(CDNAME)) THEN
  DO J=1,SIZE(YD%CNAME)
    YD%CNAME(J) = CDNAME(J)
  ENDDO
ELSE
  YD%CNAME(:) = ''
ENDIF

IF(PRESENT(PDEFAULT)) THEN
  DO J=1,SIZE(YD%REFVALI)
    YD%REFVALI(J) = PDEFAULT(J)
  ENDDO
ELSE
  YD%REFVALI(:) = 0.0_JPRB
ENDIF

IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    IF (NOFFTRAJ+YDSC%NLEVS > JPMAXSTRAJ) CALL ABOR1('JPMAXSTRAJ MUST BE INCREASED!')
    DO JLEV=1,YDSC%NLEVS
      NSTRAJGRIB(NOFFTRAJ+JLEV) = YD%IGRBCODE(JLEV)
    ENDDO
    NOFFTRAJ = NOFFTRAJ+YDSC%NLEVS
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+YDSC%NLEVS
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP3 - UNKNOWN KTRAJ')
  ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
ELSE
  YD%ITRAJ = 0
ENDIF

YD%LSET=.TRUE.
YDSC%IPTR = YDSC%IPTR+1
  
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP3

!=========================================================================

SUBROUTINE FINALISE_SFLP3(YDSURF,YDSC,YD)
! Finalise 3-D surface field group
TYPE(TSURF), INTENT(INOUT)           :: YDSURF
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT) :: YD(:)

INTEGER(KIND=JPIM) :: JFLD, JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP3',0,ZHOOK_HANDLE)

YDSC%FINALISED=.TRUE.
YDSC%NUMFLDS = YDSC%IPTR-1

WRITE(NULOUT,*) '3-D SURFACE FIELD GROUP INITIALIZED ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' NLEVS=',YDSC%NLEVS,' LMTL=',YDSC%LMTL 

YDSURF%NSURF = YDSURF%NSURF+YDSC%NUMFLDS
YDSURF%NSURFL = YDSURF%NSURFL+YDSC%NUMFLDS*YDSC%NLEVS
IF(YDSC%LMTL) THEN
  YDSURF%NPROGSURF = YDSURF%NPROGSURF+YDSC%NUMFLDS
  YDSURF%NPROGSURFL = YDSURF%NPROGSURFL+YDSC%NUMFLDS*YDSC%NLEVS
ENDIF

IF(YDSC%LMTL) THEN
  IF (LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF
YDSURF%NDIMSURF = YDSURF%NDIMSURF + YDSC%NDIM    
YDSURF%NDIMSURFL = YDSURF%NDIMSURFL + YDSC%NDIM*YDSC%NLEVS

DO JFLD=1,YDSC%NUMFLDS

  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF

  DO JLEV=1,YDSC%NLEVS
    IF(YDSC%LMTL) THEN
      WRITE(NULOUT,'(1X,A,2I4,1X,A,6I6)') &
       & YDSC%CGRPNAME(1:6),JFLD,JLEV,YD(JFLD)%CNAME(JLEV),YD(JFLD)%IGRBCODE(JLEV),&
       & YD(JFLD)%MP0,YD(JFLD)%MP9,YD(JFLD)%MP1,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN(JLEV)
    ELSE
      WRITE(NULOUT,'(1X,A,2I4,1X,A,4I6)') &
       & YDSC%CGRPNAME(1:6),JFLD,JLEV,YD(JFLD)%CNAME(JLEV),YD(JFLD)%IGRBCODE(JLEV),&
       & YD(JFLD)%MP,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN(JLEV)
    ENDIF
  ENDDO

ENDDO

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE FINALISE_SFLP3

!=========================================================================

SUBROUTINE INI_SFLP2(YDSC,YD,LDMTL,CDGRPNAME)
! Initialize 2-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT) :: YD(:)
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP2',0,ZHOOK_HANDLE)

YDSC%FINALISED = .FALSE.
YDSC%NLEVS = -1
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST

DO JFLD=1,SIZE(YD)
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
  YD(JFLD)%LSET = .FALSE.
  YD(JFLD)%IBUGFINDER = JFLD  ! For debug use only.
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 2-D SURFACE FIELD GROUP ', YDSC%CGRPNAME

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP2

!=========================================================================

SUBROUTINE SETUP_SFLP2(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 2-D surface field
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN

INTEGER(KIND=JPIM) :: IPTR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP2',0,ZHOOK_HANDLE)

IPTR = YDSC%IPTR
IF(IPTR > JPMAXSFLDS-1) THEN
  ! We are about to set the last element which is generally used as a
  ! "missing field" placeholder.
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSINED - GROUP ',YDSC%CGRPNAME,&
   &  KGRIB,CDNAME
  CALL ABOR1('IPTR > JPMAXSFLDS-1')
ENDIF

IF (YDSC%FINALISED) THEN
  WRITE(NULERR,*) 'ERROR - SETUP_SFLP2 CALLED FOR A GROUP WHICH HAS BEEN FINALISED',&
 &  YDSC%CGRPNAME,KGRIB,CDNAME
  CALL ABOR1('SETUP_SFLP2 CALLED FOR A GROUP WHICH HAS BEEN FINALISED')
ENDIF

! If the YD(:) elements are not being set in order it's a possible sign of a 
! bug - the caller may have passed in a pointer to the wrong element. Or it
! might be entirely intentional and an abort is unnecessary.
! This if-block and all other references to IBUGFINDER can be commented if
! they cause a problem.
IF (YD%IBUGFINDER /= IPTR) THEN
  WRITE(NULERR,*) 'SETUP_SFLP2 NOT BEING CALLED FOR YD ELEMENTS IN ORDER. '//&
    & '-POSSIBLY- A BUG.',YDSC%CGRPNAME,KGRIB,CDNAME
  CALL ABOR1('SETUP_SFLP2 NOT CALLED FOR YD ELEMENTS IN ORDER. POSSIBLE BUG.')
ENDIF

IF(PRESENT(KGRIB)) THEN
  YD%IGRBCODE = KGRIB
ELSE
  YD%IGRBCODE = -999
ENDIF

IF(PRESENT(KREQIN)) THEN
  YD%NREQIN = KREQIN
ELSE
  YD%NREQIN = -1
ENDIF

IF(PRESENT(CDNAME)) THEN
  YD%CNAME = CDNAME
ELSE
  YD%CNAME = ''
ENDIF

IF(PRESENT(PDEFAULT)) THEN
  YD%REFVALI = PDEFAULT
ELSE
  YD%REFVALI = 0.0_JPRB
ENDIF

IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    IF (NOFFTRAJ+1 > JPMAXSTRAJ) CALL ABOR1('JPMAXSTRAJ MUST BE INCREASED!')
    NSTRAJGRIB(NOFFTRAJ+1) = YD%IGRBCODE
    NOFFTRAJ = NOFFTRAJ+1
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+1
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP2 - UNKNOWN KTRAJ')
  ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
ELSE
  YD%ITRAJ = 0
ENDIF

YD%LSET=.TRUE.
YDSC%IPTR = YDSC%IPTR+1

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP2

!=========================================================================

SUBROUTINE FINALISE_SFLP2(YDSURF,YDSC,YD)
! Finalise 2-D surface field group
TYPE(TSURF), INTENT(INOUT)           :: YDSURF
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT) :: YD(:)

INTEGER(KIND=JPIM) :: JFLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP2',0,ZHOOK_HANDLE)

YDSC%FINALISED=.TRUE.
YDSC%NUMFLDS = YDSC%IPTR-1

WRITE(NULOUT,*) '2-D SURFACE FIELD GROUP INITIALIZED ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' LMTL=',YDSC%LMTL 

YDSURF%NSURF = YDSURF%NSURF+YDSC%NUMFLDS
YDSURF%NSURFL = YDSURF%NSURFL+YDSC%NUMFLDS
IF(YDSC%LMTL) THEN
  YDSURF%NPROGSURF = YDSURF%NPROGSURF+YDSC%NUMFLDS
  YDSURF%NPROGSURFL = YDSURF%NPROGSURFL+YDSC%NUMFLDS
ENDIF

IF(YDSC%LMTL) THEN
  IF (LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF
YDSURF%NDIMSURF = YDSURF%NDIMSURF + YDSC%NDIM    
YDSURF%NDIMSURFL = YDSURF%NDIMSURFL + YDSC%NDIM

DO JFLD=1,YDSC%NUMFLDS

  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF

  IF(YDSC%LMTL) THEN
    WRITE(NULOUT,'(1X,A,I4,1X,A,1x,I6,x,5I6)') &
     & YDSC%CGRPNAME(1:6),JFLD,YD(JFLD)%CNAME,YD(JFLD)%IGRBCODE,&
     & YD(JFLD)%MP0,YD(JFLD)%MP9,YD(JFLD)%MP1,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN
  ELSE
    WRITE(NULOUT,'(1X,A,I4,1X,A,1x,I6,x,3I6)') &
     & YDSC%CGRPNAME(1:6),JFLD,YD(JFLD)%CNAME,YD(JFLD)%IGRBCODE,&
     & YD(JFLD)%MP,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN
  ENDIF

ENDDO

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE FINALISE_SFLP2

!=========================================================================

SUBROUTINE GPPOPER(CDACT,YDSURF,KBL,PSP_SB,PSP_SG,PSP_SL,PSP_RR,PSP_CL,PSP_OM,PSP_EP,PSP_X2,YDCOM) !KPP
! Operations on prognostic surface fields
CHARACTER(LEN=*),            INTENT(IN)    :: CDACT
TYPE(TSURF),                 INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)    :: KBL
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_CL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_OM(:,:,:)                          !KPP
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_EP(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_X2(:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPPOPER',0,ZHOOK_HANDLE)
IF(PRESENT(KBL)) THEN
  CALL GPOPER_3(CDACT,YDSURF%SP_SB(:,:,:,KBL),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM)
  CALL GPOPER_2(CDACT,YDSURF%SP_SG(:,:,KBL)  ,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM)
  CALL GPOPER_2(CDACT,YDSURF%SP_SL(:,:,KBL)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM)
  CALL GPOPER_2(CDACT,YDSURF%SP_RR(:,:,KBL)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM)
  CALL GPOPER_2(CDACT,YDSURF%SP_CL(:,:,KBL)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM)
  CALL GPOPER_3(CDACT,YDSURF%SP_OM(:,:,:,KBL),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM) !KPP
  CALL GPOPER_3(CDACT,YDSURF%SP_EP(:,:,:,KBL),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM)
  CALL GPOPER_2(CDACT,YDSURF%SP_X2(:,:,KBL)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM)
ELSE
  IF(PRESENT(PSP_SB)) CALL GPOPER_3(CDACT,PSP_SB(:,:,:),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM)
  IF(PRESENT(PSP_SG)) CALL GPOPER_2(CDACT,PSP_SG(:,:)  ,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM)
  IF(PRESENT(PSP_SL)) CALL GPOPER_2(CDACT,PSP_SL(:,:)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM)
  IF(PRESENT(PSP_RR)) CALL GPOPER_2(CDACT,PSP_RR(:,:)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM)
  IF(PRESENT(PSP_CL)) CALL GPOPER_2(CDACT,PSP_CL(:,:)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM)
  IF(PRESENT(PSP_OM)) CALL GPOPER_3(CDACT,PSP_OM(:,:,:),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM) !KPP
  IF(PRESENT(PSP_EP)) CALL GPOPER_3(CDACT,PSP_EP(:,:,:),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM)
  IF(PRESENT(PSP_X2)) CALL GPOPER_2(CDACT,PSP_X2(:,:)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM)
ENDIF
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPPOPER

!=========================================================================

SUBROUTINE GPOPER(CDACT,YDSURF,KBL,PSP_SB,PSP_SG,PSP_SL,PSP_RR,PSP_CL,PSP_OM,&
 & PSD_VF,PSD_VD,PSD_VV,PSD_WS,PSD_V2,PSD_V3,YDCOM,PFIELD,PFIELD2)

!Operations on ALL surface groups
CHARACTER(LEN=*),            INTENT(IN)    :: CDACT
TYPE(TSURF),                 INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)    :: KBL
! pronostic groups:
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_CL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_OM(:,:,:)
! diagnostic groups:
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VF(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VD(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VV(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_WS(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_V2(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_V3(:,:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT)   :: YDCOM
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PFIELD(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PFIELD2(:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YRDIM%NPROMA)
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS'   .OR.&
 & CDACT == 'TRAJSTORE'  .OR. CDACT == 'TRAJSTORECST' .OR. &
 & CDACT == 'SET0TOTRAJ' .OR. CDACT == 'GETTRAJ'          ) THEN
  IF(.NOT.PRESENT(PFIELD)) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - PFIELD MISSING')
  IF(SIZE(PFIELD,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD,1) < NPROMA)')
ENDIF
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS') THEN
  IF(SIZE(PFIELD,2) < YDSURF%NPROGSURFL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD,2) < NPROGSURFL)')
ENDIF
IF(CDACT == 'GETTRAJ') THEN
  IF(.NOT.PRESENT(PFIELD2)) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - PFIELD2 MISSING')
  IF(SIZE(PFIELD2,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD2,1) < NPROMA)')
ENDIF
IF(PRESENT(YDCOM)) THEN
  YDCOM%L_OK = .FALSE.
  YDCOM%IPTRSURF = 0
  YDCOM%ICOUNT = 0
ENDIF

NPTRSURF = 0
IF(PRESENT(KBL)) THEN
  ! pronostic groups:
  IF(YDSURF%YSP_SBD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,YDSURF%SP_SB(:,:,:,KBL),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SGD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SP_SG(:,:,KBL)  ,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SLD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SP_SL(:,:,KBL)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_RRD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SP_RR(:,:,KBL)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_CLD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SP_CL(:,:,KBL)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_OMD%NDIM > 0) THEN !KPP
    CALL GPOPER_3(CDACT,YDSURF%SP_OM(:,:,:,KBL),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM,PFIELD,PFIELD2) 
  ENDIF 
  IF(YDSURF%YSP_EPD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,YDSURF%SP_EP(:,:,:,KBL),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_X2D%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SP_X2(:,:,KBL)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  ! diagnostic groups:
  IF(YDSURF%YSD_VFD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VF(:,:,KBL)  ,YDSURF%YSD_VFD,YDSURF%YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VPD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VP(:,:,KBL)  ,YDSURF%YSD_VPD,YDSURF%YSD_VP%YVP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VVD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VV(:,:,KBL)  ,YDSURF%YSD_VVD,YDSURF%YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VND%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VN(:,:,KBL)  ,YDSURF%YSD_VND,YDSURF%YSD_VN%YVN,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VHD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VH(:,:,KBL)  ,YDSURF%YSD_VHD,YDSURF%YSD_VH%YVH,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VAD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VA(:,:,KBL)  ,YDSURF%YSD_VAD,YDSURF%YSD_VA%YVA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VCD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VC(:,:,KBL)  ,YDSURF%YSD_VCD,YDSURF%YSD_VC%YVC,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_V2D%NDIM > 0) THEN !KPP
    CALL GPOPER_2(CDACT,YDSURF%SD_V2(:,:,KBL)  ,YDSURF%YSD_V2D,YDSURF%YSD_V2%YV2,YDCOM,PFIELD,PFIELD2)
  ENDIF  
  IF(YDSURF%YSD_V3D%NDIM > 0) THEN !KPP
    CALL GPOPER_3(CDACT,YDSURF%SD_V3(:,:,:,KBL),YDSURF%YSD_V3D,YDSURF%YSD_V3%YV3,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VDD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VD(:,:,KBL)  ,YDSURF%YSD_VDD,YDSURF%YSD_VD%YVD,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WSD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_WS(:,:,KBL)  ,YDSURF%YSD_WSD,YDSURF%YSD_WS%YWS,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WWD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_WW(:,:,KBL)  ,YDSURF%YSD_WWD,YDSURF%YSD_WW%YWW,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VXD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_VX(:,:,KBL)  ,YDSURF%YSD_VXD,YDSURF%YSD_VX%YVX,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XAD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,YDSURF%SD_XA(:,:,:,KBL),YDSURF%YSD_XAD,YDSURF%YSD_XA%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XRD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,YDSURF%SD_XR(:,:,:,KBL),YDSURF%YSD_XRD,YDSURF%YSD_XR%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_X2D%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_X2(:,:,KBL)  ,YDSURF%YSD_X2D,YDSURF%YSD_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SFLD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_SFL(:,:,KBL)  ,YDSURF%YSD_SFLD,YDSURF%YSD_SFL%YSFL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SFOD%NDIM > 0) THEN
    CALL GPOPER_2(CDACT,YDSURF%SD_SFO(:,:,KBL)  ,YDSURF%YSD_SFOD,YDSURF%YSD_SFO%YSFO,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_PFD%NDIM > 0) THEN
    CALL GPOPER_3(CDACT,YDSURF%SD_PF(:,:,:,KBL),YDSURF%YSD_PFD,YDSURF%YSD_PF%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
ELSE
  ! pronostic groups:
  IF(YDSURF%YSP_SBD%NDIM > 0) THEN
    IF(PRESENT(PSP_SB)) &
     & CALL GPOPER_3(CDACT,PSP_SB,YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SGD%NDIM > 0) THEN
    IF(PRESENT(PSP_SG)) &
     & CALL GPOPER_2(CDACT,PSP_SG,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SLD%NDIM > 0) THEN             !FLAKE
    IF(PRESENT(PSP_SL)) &
     & CALL GPOPER_2(CDACT,PSP_SL,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM,PFIELD,PFIELD2)
  ENDIF  
  IF(YDSURF%YSP_RRD%NDIM > 0) THEN
    IF(PRESENT(PSP_RR)) &
     & CALL GPOPER_2(CDACT,PSP_RR,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_CLD%NDIM > 0) THEN
    IF(PRESENT(PSP_CL)) &
     & CALL GPOPER_2(CDACT,PSP_CL,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_OMD%NDIM > 0) THEN
    IF(PRESENT(PSP_OM)) &
    & CALL GPOPER_3(CDACT,YDSURF%SP_OM(:,:,:,KBL),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM,PFIELD,PFIELD2)
  ENDIF
  ! diagnostic groups:
  IF(YDSURF%YSD_VFD%NDIM > 0) THEN
    IF(PRESENT(PSD_VF)) &
     & CALL GPOPER_2(CDACT,PSD_VF,YDSURF%YSD_VFD,YDSURF%YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VVD%NDIM > 0) THEN
    IF(PRESENT(PSD_VV)) &
     & CALL GPOPER_2(CDACT,PSD_VV,YDSURF%YSD_VVD,YDSURF%YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VDD%NDIM > 0) THEN
    IF(PRESENT(PSD_VD)) &
     & CALL GPOPER_2(CDACT,PSD_VD,YDSURF%YSD_VDD,YDSURF%YSD_VD%YVD,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WSD%NDIM > 0) THEN
    IF(PRESENT(PSD_WS)) &
     & CALL GPOPER_2(CDACT,PSD_WS,YDSURF%YSD_WSD,YDSURF%YSD_WS%YWS,YDCOM,PFIELD,PFIELD2)
  ENDIF
ENDIF
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER

!=========================================================================

SUBROUTINE GPOPER_2(CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2)
! Operations on 2-D surface groups
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB),INTENT(INOUT)              :: PFLD(:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)

INTEGER(KIND=JPIM) :: J,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_2',0,ZHOOK_HANDLE)
ASSOCIATE(REPSP1=>YRDYN%REPSP1)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP9)+PFLD(:,YD(J)%MP1)
    PFLD(:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = REPSP1*PFLD(:,YD(J)%MP1)+ZZPHY*PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+REPSP1*PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP0)+ZZPHY *PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETPERTOVAL') THEN
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    IF(YD(J)%NREQIN == -1) THEN
      PFLD(:,YD(J)%MP) = YD(J)%REFVALI
    ENDIF
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFIELD(:,IPTR) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFIELD(:,IPTR2) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP) = PFIELD(:,IPTR)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP) = PFIELD(:,IPTR)
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFLD(:,YD(J)%MP) = PFIELD2(:,IPTR2)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFIELD(:,NPTRSURF) = PFLD(:,J)
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFLD(:,J) = PFIELD(:,NPTRSURF)
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%IGRBCODE == YDCOM%IGRBCODE) THEN
      YDCOM%IFLDNUM  = YDCOM%IPTRSURF
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
      PFIELD(:,1) = PFLD(:,J)
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  DO J=1,YDSC%NUMFLDS  
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%NREQIN == 1) THEN
      YDCOM%ICOUNT = YDCOM%ICOUNT+1
      YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE
    ENDIF
  ENDDO
ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPOPER_2 UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPOPER_2 - UNKNOWN ACTION')
ENDIF
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_2',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_2

!=========================================================================

SUBROUTINE GPOPER_3(CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2)
! Operations on 3-D surface groups
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB),INTENT(INOUT)              :: PFLD(:,:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)

INTEGER(KIND=JPIM) :: J,JLEV,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_3',0,ZHOOK_HANDLE)
ASSOCIATE(REPSP1=>YRDYN%REPSP1)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP9)+PFLD(:,:,YD(J)%MP1)
    PFLD(:,:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = REPSP1*PFLD(:,:,YD(J)%MP1)+ZZPHY*PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+REPSP1*PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP0)+ZZPHY *PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETPERTOVAL') THEN
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    DO JLEV=1,YDSC%NLEVS
      IF(YD(J)%NREQIN(JLEV) == -1) THEN
        PFLD(:,JLEV,YD(J)%MP) = YD(J)%REFVALI(JLEV)
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFIELD(:,IPTR) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFIELD(:,IPTR2) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD(:,IPTR)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD(:,IPTR)
        ENDDO
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD2(:,IPTR2)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFIELD(:,NPTRSURF) = PFLD(:,JLEV,J)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFLD(:,JLEV,J) = PFIELD(:,NPTRSURF)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YD(J)%IGRBCODE(JLEV) == YDCOM%IGRBCODE) THEN
        YDCOM%IFLDNUM  = YDCOM%IPTRSURF
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
        PFIELD(:,1) = PFLD(:,JLEV,J)
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  DO J=1,YDSC%NUMFLDS  
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YD(J)%NREQIN(JLEV) == 1) THEN
        YDCOM%ICOUNT = YDCOM%ICOUNT+1
        YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE(JLEV)
      ENDIF
    ENDDO
  ENDDO
ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPOPER_3 UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPOPER_3 - UNKNOWN ACTION')
ENDIF
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_3',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_3

!=========================================================================

SUBROUTINE SURF_STORE(YDSURF)
! Store all surface fields
TYPE(TSURF), INTENT(INOUT) :: YDSURF

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_STORE',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YRDIM%NPROMA, NGPBLKS=>YRDIM%NGPBLKS)
ALLOCATE(YDSURF%STORE_ARRAY(NPROMA,YDSURF%NDIMSURFL,NGPBLKS))
DO JBL=1,NGPBLKS
  CALL GPOPER('GETALLFLDS',YDSURF,KBL=JBL,PFIELD=YDSURF%STORE_ARRAY(:,:,JBL))
ENDDO
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_STORE',1,ZHOOK_HANDLE)
END SUBROUTINE SURF_STORE

!=========================================================================

SUBROUTINE SURF_RESTORE(YDSURF)
! Restore all surface fields
TYPE(TSURF), INTENT(INOUT) :: YDSURF

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_RESTORE',0,ZHOOK_HANDLE)
ASSOCIATE(NGPBLKS=>YRDIM%NGPBLKS)
IF(.NOT. ALLOCATED(YDSURF%STORE_ARRAY)) &
 & CALL ABOR1('SURFACE_FIELDS_MIX:SURF_RESTORE - SURF_STORE NOT ALLOCATED')
DO JBL=1,NGPBLKS
  CALL GPOPER('PUTALLFLDS',YDSURF,KBL=JBL,PFIELD=YDSURF%STORE_ARRAY(:,:,JBL))
ENDDO
DEALLOCATE(YDSURF%STORE_ARRAY)
END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_RESTORE',1,ZHOOK_HANDLE)

END SUBROUTINE SURF_RESTORE

!=========================================================================

SUBROUTINE ALLO_SURF(YDSURF)
! Allocate surface field arrays
TYPE(TSURF), POINTER, INTENT(INOUT) :: YDSURF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ALLO_SURF',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YRDIM%NPROMA, NGPBLKS=>YRDIM%NGPBLKS)

IF (.NOT. ASSOCIATED(YDSURF)) ALLOCATE(YDSURF)

! pronostic groups:
ALLOCATE(YDSURF%SP_SB(NPROMA,YDSURF%YSP_SBD%NLEVS,YDSURF%YSP_SBD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_SG(NPROMA,YDSURF%YSP_SGD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_SL(NPROMA,YDSURF%YSP_SLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_RR(NPROMA,YDSURF%YSP_RRD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_CL(NPROMA,YDSURF%YSP_CLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_OM(NPROMA,YDSURF%YSP_OMD%NLEVS,YDSURF%YSP_OMD%NDIM,NGPBLKS)) !KPP
ALLOCATE(YDSURF%SP_EP(NPROMA,YDSURF%YSP_EPD%NLEVS,YDSURF%YSP_EPD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_X2(NPROMA,YDSURF%YSP_X2D%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_CI(NPROMA,YDSURF%YSP_CID%NDIM,NGPBLKS))
! diagnostic groups:
ALLOCATE(YDSURF%SD_VF(NPROMA,YDSURF%YSD_VFD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VP(NPROMA,YDSURF%YSD_VPD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VV(NPROMA,YDSURF%YSD_VVD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VN(NPROMA,YDSURF%YSD_VND%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VH(NPROMA,YDSURF%YSD_VHD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VA(NPROMA,YDSURF%YSD_VAD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VC(NPROMA,YDSURF%YSD_VCD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_V2(NPROMA,YDSURF%YSD_V2D%NDIM,NGPBLKS)) !KPP
ALLOCATE(YDSURF%SD_V3(NPROMA,YDSURF%YSD_V3D%NLEVS,YDSURF%YSD_V3D%NDIM,NGPBLKS)) !KPP 
ALLOCATE(YDSURF%SD_VD(NPROMA,YDSURF%YSD_VDD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_WS(NPROMA,YDSURF%YSD_WSD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_WW(NPROMA,YDSURF%YSD_WWD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VX(NPROMA,YDSURF%YSD_VXD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XA(NPROMA,YDSURF%YSD_XAD%NLEVS,YDSURF%YSD_XAD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XR(NPROMA,YDSURF%YSD_XRD%NLEVS,YDSURF%YSD_XRD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_X2(NPROMA,YDSURF%YSD_X2D%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_SFL(NPROMA,YDSURF%YSD_SFLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_SFO(NPROMA,YDSURF%YSD_SFOD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_PF(NPROMA,YDSURF%YSD_PFD%NLEVS,YDSURF%YSD_PFD%NDIM,NGPBLKS))

IF(ALLOCATED(YDSURF%SP_RR)) YDSURF%SP_RR = 0.0_JPRB
IF(ALLOCATED(YDSURF%SP_SG)) YDSURF%SP_SG = 0.0_JPRB

END ASSOCIATE
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ALLO_SURF',1,ZHOOK_HANDLE)
END SUBROUTINE ALLO_SURF

!=========================================================================

SUBROUTINE ZERO_SURF(SELF)
 
IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF

REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ZERO_SURF',0,ZHOOK_HANDLE)

! pronostic groups:
IF (ALLOCATED(SELF%SP_SB))  SELF%SP_SB  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_SG))  SELF%SP_SG  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_SL))  SELF%SP_SL  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_RR))  SELF%SP_RR  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_CL))  SELF%SP_CL  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_OM))  SELF%SP_OM  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_EP))  SELF%SP_EP  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_X2))  SELF%SP_X2  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_CI))  SELF%SP_CI  = 0.0_JPRB

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF))  SELF%SD_VF  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VP))  SELF%SD_VP  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VV))  SELF%SD_VV  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VN))  SELF%SD_VN  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VH))  SELF%SD_VH  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VA))  SELF%SD_VA  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VC))  SELF%SD_VC  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_V2))  SELF%SD_V2  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_V3))  SELF%SD_V3  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VD))  SELF%SD_VD  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_WS))  SELF%SD_WS  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_WW))  SELF%SD_WW  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VX))  SELF%SD_VX  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XA))  SELF%SD_XA  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XR))  SELF%SD_XR  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_X2))  SELF%SD_X2  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_SFL)) SELF%SD_SFL = 0.0_JPRB
IF (ALLOCATED(SELF%SD_SFO)) SELF%SD_SFO = 0.0_JPRB
IF (ALLOCATED(SELF%SD_PF)) SELF%SD_PF = 0.0_JPRB

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ZERO_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE ZERO_SURF

!=========================================================================

SUBROUTINE COPY_SURF(SELF,RHS)
 
IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF
TYPE(TSURF), INTENT(IN)    :: RHS

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_SURF',0,ZHOOK_HANDLE)
 
! pronostic groups:
IF (ALLOCATED(SELF%SP_SB) .AND. ALLOCATED(RHS%SP_SB)) THEN
  IF (ALL(SHAPE(SELF%SP_SB) == SHAPE(RHS%SP_SB))) THEN
    SELF%SP_SB = RHS%SP_SB
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SB different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SG) .AND. ALLOCATED(RHS%SP_SG)) THEN
  IF (ALL(SHAPE(SELF%SP_SG) == SHAPE(RHS%SP_SG))) THEN
    SELF%SP_SG = RHS%SP_SG
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SG different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SL) .AND. ALLOCATED(RHS%SP_SL)) THEN
  IF (ALL(SHAPE(SELF%SP_SL) == SHAPE(RHS%SP_SL))) THEN
    SELF%SP_SL = RHS%SP_SL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SL different shapes")
  ENDIF
ENDIF


IF (ALLOCATED(SELF%SP_RR) .AND. ALLOCATED(RHS%SP_RR)) THEN
  IF (ALL(SHAPE(SELF%SP_RR) == SHAPE(RHS%SP_RR))) THEN
    SELF%SP_RR = RHS%SP_RR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_RR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CL) .AND. ALLOCATED(RHS%SP_CL)) THEN
  IF (ALL(SHAPE(SELF%SP_CL) == SHAPE(RHS%SP_CL))) THEN
    SELF%SP_CL = RHS%SP_CL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_CL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_OM) .AND. ALLOCATED(RHS%SP_OM)) THEN
  IF (ALL(SHAPE(SELF%SP_OM) == SHAPE(RHS%SP_OM))) THEN
    SELF%SP_OM = RHS%SP_OM
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_OM different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_EP) .AND. ALLOCATED(RHS%SP_EP)) THEN
  IF (ALL(SHAPE(SELF%SP_EP) == SHAPE(RHS%SP_EP))) THEN
    SELF%SP_EP = RHS%SP_EP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_EP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_X2) .AND. ALLOCATED(RHS%SP_X2)) THEN
  IF (ALL(SHAPE(SELF%SP_X2) == SHAPE(RHS%SP_X2))) THEN
    SELF%SP_X2 = RHS%SP_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CI) .AND. ALLOCATED(RHS%SP_CI)) THEN
  IF (ALL(SHAPE(SELF%SP_CI) == SHAPE(RHS%SP_CI))) THEN
    SELF%SP_CI = RHS%SP_CI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_CI different shapes")
  ENDIF
ENDIF

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF) .AND. ALLOCATED(RHS%SD_VF)) THEN
  IF (ALL(SHAPE(SELF%SD_VF) == SHAPE(RHS%SD_VF))) THEN
    SELF%SD_VF= RHS%SD_VF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VF different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VP) .AND. ALLOCATED(RHS%SD_VP)) THEN
  IF (ALL(SHAPE(SELF%SD_VP) == SHAPE(RHS%SD_VP))) THEN
    SELF%SD_VP= RHS%SD_VP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VV) .AND. ALLOCATED(RHS%SD_VV)) THEN
  IF (ALL(SHAPE(SELF%SD_VV) == SHAPE(RHS%SD_VV))) THEN
    SELF%SD_VV = RHS%SD_VV
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VV different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VN) .AND. ALLOCATED(RHS%SD_VN)) THEN
  IF (ALL(SHAPE(SELF%SD_VN) == SHAPE(RHS%SD_VN))) THEN
    SELF%SD_VN = RHS%SD_VN
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VN different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VH) .AND. ALLOCATED(RHS%SD_VH)) THEN
  IF (ALL(SHAPE(SELF%SD_VH) == SHAPE(RHS%SD_VH))) THEN
    SELF%SD_VH = RHS%SD_VH
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VH different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VA) .AND. ALLOCATED(RHS%SD_VA)) THEN
  IF (ALL(SHAPE(SELF%SD_VA) == SHAPE(RHS%SD_VA))) THEN
    SELF%SD_VA = RHS%SD_VA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VC) .AND. ALLOCATED(RHS%SD_VC)) THEN
  IF (ALL(SHAPE(SELF%SD_VC) == SHAPE(RHS%SD_VC))) THEN
    SELF%SD_VC = RHS%SD_VC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V2) .AND. ALLOCATED(RHS%SD_V2)) THEN
  IF (ALL(SHAPE(SELF%SD_V2) == SHAPE(RHS%SD_V2))) THEN
    SELF%SD_V2 = RHS%SD_V2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_V2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V3) .AND. ALLOCATED(RHS%SD_V3)) THEN
  IF (ALL(SHAPE(SELF%SD_V3) == SHAPE(RHS%SD_V3))) THEN
    SELF%SD_V3 = RHS%SD_V3
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_V3 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VD) .AND. ALLOCATED(RHS%SD_VD)) THEN
  IF (ALL(SHAPE(SELF%SD_VD) == SHAPE(RHS%SD_VD))) THEN
    SELF%SD_VD = RHS%SD_VD
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VD different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WS) .AND. ALLOCATED(RHS%SD_WS)) THEN
  IF (ALL(SHAPE(SELF%SD_WS) == SHAPE(RHS%SD_WS))) THEN
    SELF%SD_WS = RHS%SD_WS
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_WS different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WW) .AND. ALLOCATED(RHS%SD_WW)) THEN
  IF (ALL(SHAPE(SELF%SD_WW) == SHAPE(RHS%SD_WW))) THEN
    SELF%SD_WW = RHS%SD_WW
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_WW different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VX) .AND. ALLOCATED(RHS%SD_VX)) THEN
  IF (ALL(SHAPE(SELF%SD_VX) == SHAPE(RHS%SD_VX))) THEN
    SELF%SD_VX = RHS%SD_VX
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VX different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XA) .AND. ALLOCATED(RHS%SD_XA)) THEN
  IF (ALL(SHAPE(SELF%SD_XA) == SHAPE(RHS%SD_XA))) THEN
    SELF%SD_XA = RHS%SD_XA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XR) .AND. ALLOCATED(RHS%SD_XR)) THEN
  IF (ALL(SHAPE(SELF%SD_XR) == SHAPE(RHS%SD_XR))) THEN
    SELF%SD_XR = RHS%SD_XR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_X2) .AND. ALLOCATED(RHS%SD_X2)) THEN
  IF (ALL(SHAPE(SELF%SD_X2) == SHAPE(RHS%SD_X2))) THEN
    SELF%SD_X2 = RHS%SD_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFL) .AND. ALLOCATED(RHS%SD_SFL)) THEN
  IF (ALL(SHAPE(SELF%SD_SFL) == SHAPE(RHS%SD_SFL))) THEN
    SELF%SD_SFL = RHS%SD_SFL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_SFL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFO) .AND. ALLOCATED(RHS%SD_SFO)) THEN
  IF (ALL(SHAPE(SELF%SD_SFO) == SHAPE(RHS%SD_SFO))) THEN
    SELF%SD_SFO = RHS%SD_SFO
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_SFO different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_PF) .AND. ALLOCATED(RHS%SD_PF)) THEN
  IF (ALL(SHAPE(SELF%SD_PF) == SHAPE(RHS%SD_PF))) THEN
    SELF%SD_PF = RHS%SD_PF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_PF different shapes")
  ENDIF
ENDIF


!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_SURF

!=========================================================================

SUBROUTINE MUL_SURF(SELF,PZ)
 
IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF
REAL(KIND=JPRB),   INTENT(IN)    :: PZ

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:MUL_SURF',0,ZHOOK_HANDLE)
! pronostic groups:
IF (ALLOCATED(SELF%SP_SB))  SELF%SP_SB  = PZ*SELF%SP_SB
IF (ALLOCATED(SELF%SP_SG))  SELF%SP_SG  = PZ*SELF%SP_SG
IF (ALLOCATED(SELF%SP_SL))  SELF%SP_SL  = PZ*SELF%SP_SL
IF (ALLOCATED(SELF%SP_RR))  SELF%SP_RR  = PZ*SELF%SP_RR
IF (ALLOCATED(SELF%SP_CL))  SELF%SP_CL  = PZ*SELF%SP_CL
IF (ALLOCATED(SELF%SP_OM))  SELF%SP_OM  = PZ*SELF%SP_OM
IF (ALLOCATED(SELF%SP_EP))  SELF%SP_EP  = PZ*SELF%SP_EP
IF (ALLOCATED(SELF%SP_X2))  SELF%SP_X2  = PZ*SELF%SP_X2
IF (ALLOCATED(SELF%SP_CI))  SELF%SP_CI  = PZ*SELF%SP_CI

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF))  SELF%SD_VF  = PZ*SELF%SD_VF
IF (ALLOCATED(SELF%SD_VP))  SELF%SD_VP  = PZ*SELF%SD_VP
IF (ALLOCATED(SELF%SD_VV))  SELF%SD_VV  = PZ*SELF%SD_VV
IF (ALLOCATED(SELF%SD_VN))  SELF%SD_VN  = PZ*SELF%SD_VN
IF (ALLOCATED(SELF%SD_VH))  SELF%SD_VH  = PZ*SELF%SD_VH
IF (ALLOCATED(SELF%SD_VA))  SELF%SD_VA  = PZ*SELF%SD_VA
IF (ALLOCATED(SELF%SD_VC))  SELF%SD_VC  = PZ*SELF%SD_VC
IF (ALLOCATED(SELF%SD_V2))  SELF%SD_V2  = PZ*SELF%SD_V2
IF (ALLOCATED(SELF%SD_V3))  SELF%SD_V3  = PZ*SELF%SD_V3
IF (ALLOCATED(SELF%SD_VD))  SELF%SD_VD  = PZ*SELF%SD_VD
IF (ALLOCATED(SELF%SD_WS))  SELF%SD_WS  = PZ*SELF%SD_WS
IF (ALLOCATED(SELF%SD_WW))  SELF%SD_WW  = PZ*SELF%SD_WW
IF (ALLOCATED(SELF%SD_VX))  SELF%SD_VX  = PZ*SELF%SD_VX
IF (ALLOCATED(SELF%SD_XA))  SELF%SD_XA  = PZ*SELF%SD_XA
IF (ALLOCATED(SELF%SD_XR))  SELF%SD_XR  = PZ*SELF%SD_XR
IF (ALLOCATED(SELF%SD_X2))  SELF%SD_X2  = PZ*SELF%SD_X2
IF (ALLOCATED(SELF%SD_SFL)) SELF%SD_SFL = PZ*SELF%SD_SFL
IF (ALLOCATED(SELF%SD_SFO)) SELF%SD_SFO = PZ*SELF%SD_SFO
IF (ALLOCATED(SELF%SD_PF))  SELF%SD_PF  = PZ*SELF%SD_PF

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:MUL_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE MUL_SURF

!=========================================================================

SUBROUTINE AXPBY_SURF(SELF,PA,RHS,PB)
 
IMPLICIT NONE
TYPE(TSURF),     INTENT(INOUT) :: SELF
REAL(KIND=JPRB), INTENT(IN)    :: PA
TYPE(TSURF),     INTENT(IN)    :: RHS
REAL(KIND=JPRB), INTENT(IN)    :: PB

REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
#include "abor1.intfb.h"

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:AXPBY_SURF',0,ZHOOK_HANDLE)
 
! pronostic groups:
IF (ALLOCATED(SELF%SP_SB) .AND. ALLOCATED(RHS%SP_SB)) THEN
  IF (ALL(SHAPE(SELF%SP_SB) == SHAPE(RHS%SP_SB))) THEN
    SELF%SP_SB = PA*SELF%SP_SB + PB*RHS%SP_SB
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_SB different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SG) .AND. ALLOCATED(RHS%SP_SG)) THEN
  IF (ALL(SHAPE(SELF%SP_SG) == SHAPE(RHS%SP_SG))) THEN
    SELF%SP_SG = PA*SELF%SP_SG + PB*RHS%SP_SG
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_SG different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SL) .AND. ALLOCATED(RHS%SP_SL)) THEN
  IF (ALL(SHAPE(SELF%SP_SL) == SHAPE(RHS%SP_SL))) THEN
    SELF%SP_SL = PA*SELF%SP_SL + PB*RHS%SP_SL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_SL different shapes")
  ENDIF
ENDIF


IF (ALLOCATED(SELF%SP_RR) .AND. ALLOCATED(RHS%SP_RR)) THEN
  IF (ALL(SHAPE(SELF%SP_RR) == SHAPE(RHS%SP_RR))) THEN
    SELF%SP_RR = PA*SELF%SP_RR + PB*RHS%SP_RR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_RR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CL) .AND. ALLOCATED(RHS%SP_CL)) THEN
  IF (ALL(SHAPE(SELF%SP_CL) == SHAPE(RHS%SP_CL))) THEN
    SELF%SP_CL = PA*SELF%SP_CL + PB*RHS%SP_CL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_CL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_OM) .AND. ALLOCATED(RHS%SP_OM)) THEN
  IF (ALL(SHAPE(SELF%SP_OM) == SHAPE(RHS%SP_OM))) THEN
    SELF%SP_OM = PA*SELF%SP_OM + PB*RHS%SP_OM
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_OM different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_EP) .AND. ALLOCATED(RHS%SP_EP)) THEN
  IF (ALL(SHAPE(SELF%SP_EP) == SHAPE(RHS%SP_EP))) THEN
    SELF%SP_EP = PA*SELF%SP_EP + PB*RHS%SP_EP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_EP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_X2) .AND. ALLOCATED(RHS%SP_X2)) THEN
  IF (ALL(SHAPE(SELF%SP_X2) == SHAPE(RHS%SP_X2))) THEN
    SELF%SP_X2 = PA*SELF%SP_X2 + PB*RHS%SP_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CI) .AND. ALLOCATED(RHS%SP_CI)) THEN
  IF (ALL(SHAPE(SELF%SP_CI) == SHAPE(RHS%SP_CI))) THEN
    SELF%SP_CI = PA*SELF%SP_CI + PB*RHS%SP_CI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SP_CI different shapes")
  ENDIF
ENDIF

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF) .AND. ALLOCATED(RHS%SD_VF)) THEN
  IF (ALL(SHAPE(SELF%SD_VF) == SHAPE(RHS%SD_VF))) THEN
    SELF%SD_VF = PA*SELF%SD_VF + PB*RHS%SD_VF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VF different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VP) .AND. ALLOCATED(RHS%SD_VP)) THEN
  IF (ALL(SHAPE(SELF%SD_VP) == SHAPE(RHS%SD_VP))) THEN
    SELF%SD_VP = PA*SELF%SD_VP + PB*RHS%SD_VP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VV) .AND. ALLOCATED(RHS%SD_VV)) THEN
  IF (ALL(SHAPE(SELF%SD_VV) == SHAPE(RHS%SD_VV))) THEN
    SELF%SD_VV = PA*SELF%SD_VV + PB*RHS%SD_VV
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VV different shapes")
  ENDIF
ENDIF
  
IF (ALLOCATED(SELF%SD_VN) .AND. ALLOCATED(RHS%SD_VN)) THEN
  IF (ALL(SHAPE(SELF%SD_VN) == SHAPE(RHS%SD_VN))) THEN
    SELF%SD_VN = PA*SELF%SD_VN + PB*RHS%SD_VN
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VN different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VH) .AND. ALLOCATED(RHS%SD_VH)) THEN
  IF (ALL(SHAPE(SELF%SD_VH) == SHAPE(RHS%SD_VH))) THEN
    SELF%SD_VH = PA*SELF%SD_VH + PB*RHS%SD_VH
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VH different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VA) .AND. ALLOCATED(RHS%SD_VA)) THEN
  IF (ALL(SHAPE(SELF%SD_VA) == SHAPE(RHS%SD_VA))) THEN
    SELF%SD_VA = PA*SELF%SD_VA + PB*RHS%SD_VA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VC) .AND. ALLOCATED(RHS%SD_VC)) THEN
  IF (ALL(SHAPE(SELF%SD_VC) == SHAPE(RHS%SD_VC))) THEN
    SELF%SD_VC = PA*SELF%SD_VC + PB*RHS%SD_VC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V2) .AND. ALLOCATED(RHS%SD_V2)) THEN
  IF (ALL(SHAPE(SELF%SD_V2) == SHAPE(RHS%SD_V2))) THEN
    SELF%SD_V2 = PA*SELF%SD_V2 + PB*RHS%SD_V2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_V2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V3) .AND. ALLOCATED(RHS%SD_V3)) THEN
  IF (ALL(SHAPE(SELF%SD_V3) == SHAPE(RHS%SD_V3))) THEN
    SELF%SD_V3 = PA*SELF%SD_V3 + PB*RHS%SD_V3
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_V3 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VD) .AND. ALLOCATED(RHS%SD_VD)) THEN
  IF (ALL(SHAPE(SELF%SD_VD) == SHAPE(RHS%SD_VD))) THEN
    SELF%SD_VD = PA*SELF%SD_VD + PB*RHS%SD_VD
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VD different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WS) .AND. ALLOCATED(RHS%SD_WS)) THEN
  IF (ALL(SHAPE(SELF%SD_WS) == SHAPE(RHS%SD_WS))) THEN
    SELF%SD_WS = PA*SELF%SD_WS + PB*RHS%SD_WS
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_WS different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WW) .AND. ALLOCATED(RHS%SD_WW)) THEN
  IF (ALL(SHAPE(SELF%SD_WW) == SHAPE(RHS%SD_WW))) THEN
    SELF%SD_WW = PA*SELF%SD_WW + PB*RHS%SD_WW
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_WW different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VX) .AND. ALLOCATED(RHS%SD_VX)) THEN
  IF (ALL(SHAPE(SELF%SD_VX) == SHAPE(RHS%SD_VX))) THEN
    SELF%SD_VX = PA*SELF%SD_VX + PB*RHS%SD_VX
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_VX different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XA) .AND. ALLOCATED(RHS%SD_XA)) THEN
  IF (ALL(SHAPE(SELF%SD_XA) == SHAPE(RHS%SD_XA))) THEN
    SELF%SD_XA = PA*SELF%SD_XA + PB*RHS%SD_XA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_XA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XR) .AND. ALLOCATED(RHS%SD_XR)) THEN
  IF (ALL(SHAPE(SELF%SD_XR) == SHAPE(RHS%SD_XR))) THEN
    SELF%SD_XR = PA*SELF%SD_XR + PB*RHS%SD_XR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_XR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_X2) .AND. ALLOCATED(RHS%SD_X2)) THEN
  IF (ALL(SHAPE(SELF%SD_X2) == SHAPE(RHS%SD_X2))) THEN
    SELF%SD_X2 = PA*SELF%SD_X2 + PB*RHS%SD_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFL) .AND. ALLOCATED(RHS%SD_SFL)) THEN
  IF (ALL(SHAPE(SELF%SD_SFL) == SHAPE(RHS%SD_SFL))) THEN
    SELF%SD_SFL = PA*SELF%SD_SFL + PB*RHS%SD_SFL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_SFL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFO) .AND. ALLOCATED(RHS%SD_SFO)) THEN
  IF (ALL(SHAPE(SELF%SD_SFO) == SHAPE(RHS%SD_SFO))) THEN
    SELF%SD_SFO = PA*SELF%SD_SFO + PB*RHS%SD_SFO
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_SFO different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_PF) .AND. ALLOCATED(RHS%SD_PF)) THEN
  IF (ALL(SHAPE(SELF%SD_PF) == SHAPE(RHS%SD_PF))) THEN
    SELF%SD_PF = PA*SELF%SD_PF + PB*RHS%SD_PF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:AXPBY SD_PF different shapes")
  ENDIF
ENDIF


!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:AXPBY_SURF',1,ZHOOK_HANDLE)
  
END SUBROUTINE AXPBY_SURF

!=========================================================================

SUBROUTINE DOT_PROD_SURF(FLD1,FLD2,PPROD)

IMPLICIT NONE
TYPE(TSURF),     INTENT(IN)  :: FLD1
TYPE(TSURF),     INTENT(IN)  :: FLD2
REAL(KIND=JPRB), INTENT(OUT) :: PPROD
  
REAL(KIND=JPRB)    :: ZTMP
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

#include "dotprod2.intfb.h"
#include "dotprod3.intfb.h"

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DOT_PROD_SURF',0,ZHOOK_HANDLE)

PPROD = 0.0_JPRB

! pronostic groups:
IF (ALLOCATED(FLD1%SP_SB) .AND. ALLOCATED(FLD2%SP_SB)) THEN
  CALL DOTPROD3(FLD1%SP_SB,FLD2%SP_SB,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_SG) .AND. ALLOCATED(FLD2%SP_SG)) THEN
  CALL DOTPROD2(FLD1%SP_SG,FLD2%SP_SG,ZTMP)
  PPROD = PPROD + ZTMP 
ENDIF

IF (ALLOCATED(FLD1%SP_SL) .AND. ALLOCATED(FLD2%SP_SL)) THEN
  CALL DOTPROD2(FLD1%SP_SL,FLD2%SP_SL,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_RR) .AND. ALLOCATED(FLD2%SP_RR)) THEN
  CALL DOTPROD2(FLD1%SP_RR,FLD2%SP_RR,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_CL) .AND. ALLOCATED(FLD2%SP_CL)) THEN
  CALL DOTPROD2(FLD1%SP_CL,FLD2%SP_CL,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_OM) .AND. ALLOCATED(FLD2%SP_OM)) THEN
  CALL DOTPROD3(FLD1%SP_OM,FLD2%SP_OM,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_EP) .AND. ALLOCATED(FLD2%SP_EP)) THEN
  CALL DOTPROD3(FLD1%SP_EP,FLD2%SP_EP,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_X2) .AND. ALLOCATED(FLD2%SP_X2)) THEN
  CALL DOTPROD2(FLD1%SP_X2,FLD2%SP_X2,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SP_CI) .AND. ALLOCATED(FLD2%SP_CI)) THEN
  CALL DOTPROD2(FLD1%SP_CI,FLD2%SP_CI,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

! diagnostic groups:
IF (ALLOCATED(FLD1%SD_VF) .AND. ALLOCATED(FLD2%SD_VF)) THEN
  CALL DOTPROD2(FLD1%SD_VF,FLD2%SD_VF,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VP) .AND. ALLOCATED(FLD2%SD_VP)) THEN
  CALL DOTPROD2(FLD1%SD_VP,FLD2%SD_VP,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VV) .AND. ALLOCATED(FLD2%SD_VV)) THEN
  CALL DOTPROD2(FLD1%SD_VV,FLD2%SD_VV,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VN) .AND. ALLOCATED(FLD2%SD_VN)) THEN
  CALL DOTPROD2(FLD1%SD_VN,FLD2%SD_VN,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VH) .AND. ALLOCATED(FLD2%SD_VH)) THEN
  CALL DOTPROD2(FLD1%SD_VH,FLD2%SD_VH,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VA) .AND. ALLOCATED(FLD2%SD_VA)) THEN
  CALL DOTPROD2(FLD1%SD_VA,FLD2%SD_VA,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VC) .AND. ALLOCATED(FLD2%SD_VC)) THEN
  CALL DOTPROD2(FLD1%SD_VC,FLD2%SD_VC,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_V2) .AND. ALLOCATED(FLD2%SD_V2)) THEN
  CALL DOTPROD2(FLD1%SD_V2,FLD2%SD_V2,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF
  
IF (ALLOCATED(FLD1%SD_V3) .AND. ALLOCATED(FLD2%SD_V3)) THEN
  CALL DOTPROD3(FLD1%SD_V3,FLD2%SD_V3,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VD) .AND. ALLOCATED(FLD2%SD_VD)) THEN
  CALL DOTPROD2(FLD1%SD_VD,FLD2%SD_VD,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_WS) .AND. ALLOCATED(FLD2%SD_WS)) THEN
  CALL DOTPROD2(FLD1%SD_WS,FLD2%SD_WS,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_WW) .AND. ALLOCATED(FLD2%SD_WW)) THEN
  CALL DOTPROD2(FLD1%SD_WW,FLD2%SD_WW,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_VX) .AND. ALLOCATED(FLD2%SD_VX)) THEN
  CALL DOTPROD2(FLD1%SD_VX,FLD2%SD_VX,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_XA) .AND. ALLOCATED(FLD2%SD_XA)) THEN
  CALL DOTPROD3(FLD1%SD_XA,FLD2%SD_XA,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_XR) .AND. ALLOCATED(FLD2%SD_XR)) THEN
  CALL DOTPROD3(FLD1%SD_XR,FLD2%SD_XR,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_X2) .AND. ALLOCATED(FLD2%SD_X2)) THEN
  CALL DOTPROD2(FLD1%SD_X2,FLD2%SD_X2,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_SFL) .AND. ALLOCATED(FLD2%SD_SFL)) THEN
  CALL DOTPROD2(FLD1%SD_SFL,FLD2%SD_SFL,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_SFO) .AND. ALLOCATED(FLD2%SD_SFO)) THEN
  CALL DOTPROD2(FLD1%SD_SFO,FLD2%SD_SFO,ZTMP)
  PPROD = PPROD + ZTMP
ENDIF

IF (ALLOCATED(FLD1%SD_PF) .AND. ALLOCATED(FLD2%SD_PF)) THEN
 CALL DOTPROD3(FLD1%SD_PF,FLD2%SD_PF,ZTMP)
 PPROD = PPROD + ZTMP
ENDIF


!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DOT_PROD_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE DOT_PROD_SURF

!=========================================================================

SUBROUTINE RANDOM_SURF(SELF)
 
IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF

REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:RANDOM_SURF',0,ZHOOK_HANDLE)

! pronostic groups:
IF (ALLOCATED(SELF%SP_SB)) THEN
  CALL RANDOM_NUMBER(SELF%SP_SB)
  SELF%SP_SB = 2.0_JPRB*SELF%SP_SB - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_SG)) THEN
  CALL RANDOM_NUMBER(SELF%SP_SG)
  SELF%SP_SG = 2.0_JPRB*SELF%SP_SG - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_SL)) THEN
  CALL RANDOM_NUMBER(SELF%SP_SL)
  SELF%SP_SL = 2.0_JPRB*SELF%SP_SL - 1.0_JPRB
ENDIF


IF (ALLOCATED(SELF%SP_RR)) THEN
  CALL RANDOM_NUMBER(SELF%SP_RR)
  SELF%SP_RR = 2.0_JPRB*SELF%SP_RR - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_CL)) THEN
  CALL RANDOM_NUMBER(SELF%SP_CL)
  SELF%SP_CL = 2.0_JPRB*SELF%SP_CL - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_OM)) THEN
  CALL RANDOM_NUMBER(SELF%SP_OM)
  SELF%SP_OM = 2.0_JPRB*SELF%SP_OM - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_EP)) THEN
  CALL RANDOM_NUMBER(SELF%SP_EP)
  SELF%SP_EP = 2.0_JPRB*SELF%SP_EP - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_X2)) THEN
  CALL RANDOM_NUMBER(SELF%SP_X2)
  SELF%SP_X2 = 2.0_JPRB*SELF%SP_X2 - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SP_CI)) THEN
  CALL RANDOM_NUMBER(SELF%SP_CI)
  SELF%SP_CI = 2.0_JPRB*SELF%SP_CI - 1.0_JPRB
ENDIF

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VF)
  SELF%SD_VF = 2.0_JPRB*SELF%SD_VF - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VP)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VP)
  SELF%SD_VP = 2.0_JPRB*SELF%SD_VP - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VV)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VV)
  SELF%SD_VV = 2.0_JPRB*SELF%SD_VV - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VN)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VN)
  SELF%SD_VN = 2.0_JPRB*SELF%SD_VN - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VH)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VH)
  SELF%SD_VH = 2.0_JPRB*SELF%SD_VH - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VA)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VA)
  SELF%SD_VA = 2.0_JPRB*SELF%SD_VA - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VC)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VC)
  SELF%SD_VC = 2.0_JPRB*SELF%SD_VC - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_V2)) THEN
  CALL RANDOM_NUMBER(SELF%SD_V2)
  SELF%SD_V2 = 2.0_JPRB*SELF%SD_V2 - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_V3)) THEN
  CALL RANDOM_NUMBER(SELF%SD_V3)
  SELF%SD_V3 = 2.0_JPRB*SELF%SD_V3 - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VD)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VD)
  SELF%SD_VD = 2.0_JPRB*SELF%SD_VD - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_WS)) THEN
  CALL RANDOM_NUMBER(SELF%SD_WS)
  SELF%SD_WS = 2.0_JPRB*SELF%SD_WS - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_WW)) THEN
  CALL RANDOM_NUMBER(SELF%SD_WW)
  SELF%SD_WW = 2.0_JPRB*SELF%SD_WW - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_VX)) THEN
  CALL RANDOM_NUMBER(SELF%SD_VX)
  SELF%SD_VX = 2.0_JPRB*SELF%SD_VX - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_XA)) THEN
  CALL RANDOM_NUMBER(SELF%SD_XA)
  SELF%SD_XA = 2.0_JPRB*SELF%SD_XA - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_XR)) THEN
  CALL RANDOM_NUMBER(SELF%SD_XR)
  SELF%SD_XR = 2.0_JPRB*SELF%SD_XR - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_X2)) THEN
  CALL RANDOM_NUMBER(SELF%SD_X2)
  SELF%SD_X2 = 2.0_JPRB*SELF%SD_X2 - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_SFL)) THEN
  CALL RANDOM_NUMBER(SELF%SD_SFL)
  SELF%SD_SFL = 2.0_JPRB*SELF%SD_SFL - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_SFO)) THEN
  CALL RANDOM_NUMBER(SELF%SD_SFO)
  SELF%SD_SFO = 2.0_JPRB*SELF%SD_SFO - 1.0_JPRB
ENDIF

IF (ALLOCATED(SELF%SD_PF)) THEN
  CALL RANDOM_NUMBER(SELF%SD_PF)
  SELF%SD_PF = 2.0_JPRB*SELF%SD_PF - 1.0_JPRB
ENDIF

!IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:RANDOM_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE RANDOM_SURF

END MODULE SURFACE_FIELDS_MIX
