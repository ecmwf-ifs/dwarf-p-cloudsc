! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_GLOBAL_ATLAS_STATE_MOD
  ! Driver module to manage the setup and teardown of the dwarf memory state
  USE PARKIND1,  ONLY : JPIM, JPRB
  USE YOMPHYDER, ONLY : STATE_TYPE
  USE YOECLDP,   ONLY : NCLV, YRECLDP, YRECLDP_LOAD_PARAMETERS
  USE YOMCST,    ONLY : YOMCST_LOAD_PARAMETERS
  USE YOETHF,    ONLY : YOETHF_LOAD_PARAMETERS
  USE YOEPHLI  , ONLY : YREPHLI, YREPHLI_LOAD_PARAMETERS

  USE FILE_IO_MOD, ONLY: INPUT_INITIALIZE, INPUT_FINALIZE, LOAD_SCALAR, LOAD_ARRAY
  USE EXPAND_ATLAS_MOD, ONLY: LOADVAR_ATLAS, LOADSTATE_ATLAS
  USE VALIDATE_ATLAS_MOD, ONLY: VALIDATEVAR_ATLAS, VALIDATESTATE_ATLAS
  USE CLOUDSC_MPI_MOD, ONLY: IRANK

  USE ATLAS_MODULE
  USE, INTRINSIC :: ISO_C_BINDING
  USE ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS_MODULE

  IMPLICIT NONE

  TYPE VAR3D_PTR
      REAL(C_DOUBLE), POINTER :: PTR(:,:,:)
  END TYPE
  TYPE VAR2D_PTR
      REAL(C_DOUBLE), POINTER :: PTR(:,:)
  END TYPE

  !INTEGER, PARAMETER PLCRIT_AER = 1
  !INTEGER, PARAMETER PLCRIT_AER = 2 
  !IN_VAR_NAMES(PLCRIT_AER)
  CHARACTER(LEN=10), PARAMETER, DIMENSION(30) :: IN_VAR_NAMES = (/ &
      "PLCRIT_AER", "PICRIT_AER", "PRE_ICE   ", "PCCN      ", "PNICE     ", "PT        ", "PQ        ", &
      "PVFA      ", "PVFL      ", "PVFI      ", "PDYNA     ", "PDYNL     ", "PDYNI     ", "PHRSW     ", &
      "PHRLW     ", "PVERVEL   ", "PAP       ", "PLU       ", "PLUDE     ", "PSNDE     ", "PMFU      ", &
      "PMFD      ", "PA        ", "PSUPSAT   ", &
      "PLSM      ", "LDCUM     ", "KTYPE     ", "PAPH      ", "PEXTRA    ", "PCLV      " /)
  CHARACTER(LEN=16), PARAMETER, DIMENSION(16) :: OUT_VAR_NAMES = (/ &
      "PFSQLF          ", "PFSQIF          ", "PFCQLNG         ", "PFCQNNG         ", "PFSQRF          ", &
      "PFSQSF          ", "PFCQRNG         ", "PFCQSNG         ", "PFSQLTUR        ", "PFSQITUR        ", &
      "PFPLSL          ", "PFPLSN          ", "PFHPSL          ", "PFHPSN          ", "PCOVPTOT        ", &
      "PRAINFRAC_TOPRFZ" /) 

  TYPE CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW
    TYPE(VAR2D_PTR), DIMENSION(24) :: IN_VARS_2D_REAL64
    TYPE(VAR2D_PTR), DIMENSION(15) :: OUT_VARS_2D_REAL64

    ! Input field variables and tendencies
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PLCRIT_AER(:,:)
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PICRIT_AER(:,:) 
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PRE_ICE(:,:) 
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PCCN(:,:)     ! liquid cloud condensation nuclei
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PNICE(:,:)    ! ice number concentration (cf. CCN)
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PT(:,:)       ! T at start of callpar
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PQ(:,:)       ! Q at start of callpar
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PVFA(:,:)     ! CC from VDF scheme
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PVFL(:,:)     ! Liq from VDF scheme
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PVFI(:,:)     ! Ice from VDF scheme
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PDYNA(:,:)    ! CC from Dynamics
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PDYNL(:,:)    ! Liq from Dynamics
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PDYNI(:,:)    ! Liq from Dynamics
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PHRSW(:,:)    ! Short-wave heating rate
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PHRLW(:,:)    ! Long-wave heating rate
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PVERVEL(:,:)  ! Vertical velocity
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PAP(:,:)      ! Pressure on full levels
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PLU(:,:)      ! Conv. condensate
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PLUDE(:,:)    ! Conv. detrained water 
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PSNDE(:,:)    ! Conv. detrained snow
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PMFU(:,:)     ! Conv. mass flux up
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PMFD(:,:)     ! Conv. mass flux down
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PA(:,:)       ! Original Cloud fraction (t)
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PSUPSAT(:,:)

    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PLSM(:)       ! Land fraction (0-1) 
    LOGICAL,        POINTER, CONTIGUOUS :: LDCUM(:)      ! Convection active
    INTEGER,        POINTER, CONTIGUOUS :: KTYPE(:)      ! Convection type 0,1,2
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PAPH(:,:)     ! Pressure on half levels
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PEXTRA(:,:,:) ! extra fields
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PCLV(:,:,:) 

    TYPE(STATE_TYPE) :: TENDENCY_CML ! cumulative tendency used for final output
    TYPE(STATE_TYPE) :: TENDENCY_TMP ! cumulative tendency used as input
    TYPE(STATE_TYPE) :: TENDENCY_LOC ! local tendency from cloud scheme

    ! Output fields used for validation
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQLF(:,:)   ! Flux of liquid
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQIF(:,:)   ! Flux of ice
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFCQLNG(:,:)  ! -ve corr for liq
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFCQNNG(:,:)  ! -ve corr for ice
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQRF(:,:)   ! Flux diagnostics
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQSF(:,:)   ! for DDH, generic
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFCQRNG(:,:)  ! rain
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFCQSNG(:,:)  ! snow
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQLTUR(:,:) ! liquid flux due to VDF
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFSQITUR(:,:) ! ice flux due to VDF
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFPLSL(:,:)   ! liq+rain sedim flux
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFPLSN(:,:)   ! ice+snow sedim flux
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFHPSL(:,:)   ! Enthalpy flux for liq
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PFHPSN(:,:)   ! Enthalpy flux for ice
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PCOVPTOT(:,:) ! Precip fraction
    REAL(C_DOUBLE), POINTER, CONTIGUOUS :: PRAINFRAC_TOPRFZ(:) 

    CONTAINS
    PROCEDURE :: GET_BLOCK => CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK
  END TYPE CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW

  TYPE CLOUDSC_GLOBAL_ATLAS_STATE_FIELDS
    !TYPE(atlas_field), allocatlabe :: fields(:)
  
      ! Input field variables and tendencies
    TYPE(atlas_field) :: fPLCRIT_AER
    TYPE(atlas_field) :: fPICRIT_AER 
    TYPE(atlas_field) :: fPRE_ICE 
    TYPE(atlas_field) :: fPCCN     ! liquid cloud condensation nuclei
    TYPE(atlas_field) :: fPNICE    ! ice number concentration (cf. CCN)
    TYPE(atlas_field) :: fPT       ! T at start of callpar
    TYPE(atlas_field) :: fPQ       ! Q at start of callpar
    TYPE(atlas_field) :: fPVFA     ! CC from VDF scheme
    TYPE(atlas_field) :: fPVFL     ! Liq from VDF scheme
    TYPE(atlas_field) :: fPVFI     ! Ice from VDF scheme
    TYPE(atlas_field) :: fPDYNA    ! CC from Dynamics
    TYPE(atlas_field) :: fPDYNL    ! Liq from Dynamics
    TYPE(atlas_field) :: fPDYNI    ! Liq from Dynamics
    TYPE(atlas_field) :: fPHRSW    ! Short-wave heating rate
    TYPE(atlas_field) :: fPHRLW    ! Long-wave heating rate
    TYPE(atlas_field) :: fPVERVEL  ! Vertical velocity
    TYPE(atlas_field) :: fPAP      ! Pressure on full levels
    TYPE(atlas_field) :: fPLU      ! Conv. condensate
    TYPE(atlas_field) :: fPLUDE    ! Conv. detrained water 
    TYPE(atlas_field) :: fPSNDE    ! Conv. detrained snow
    TYPE(atlas_field) :: fPMFU     ! Conv. mass flux up
    TYPE(atlas_field) :: fPMFD     ! Conv. mass flux down
    TYPE(atlas_field) :: fPA       ! Original Cloud fraction (t)
    TYPE(atlas_field) :: fPSUPSAT

    TYPE(atlas_field) :: fPLSM       ! Land fraction (0-1) 
    TYPE(atlas_field) :: fLDCUM      ! Convection active
    TYPE(atlas_field) :: fKTYPE      ! Convection type 0,1,2
    TYPE(atlas_field) :: fPAPH     ! Pressure on half levels
    TYPE(atlas_field) :: fPEXTRA ! extra fields
    TYPE(atlas_field) :: fPCLV 

    TYPE(atlas_field)  :: fTENDENCY_CML ! cumulative tendency used for final output
    TYPE(atlas_field)  :: fTENDENCY_TMP ! cumulative tendency used as input
    TYPE(atlas_field)  :: fTENDENCY_LOC ! local tendency from cloud scheme

    ! Output fields used for validation
    TYPE(atlas_field) :: fPFSQLF   ! Flux of liquid
    TYPE(atlas_field) :: fPFSQIF   ! Flux of ice
    TYPE(atlas_field) :: fPFCQLNG  ! -ve corr for liq
    TYPE(atlas_field) :: fPFCQNNG  ! -ve corr for ice
    TYPE(atlas_field) :: fPFSQRF   ! Flux diagnostics
    TYPE(atlas_field) :: fPFSQSF   ! for DDH, generic
    TYPE(atlas_field) :: fPFCQRNG  ! rain
    TYPE(atlas_field) :: fPFCQSNG  ! snow
    TYPE(atlas_field) :: fPFSQLTUR ! liquid flux due to VDF
    TYPE(atlas_field) :: fPFSQITUR ! ice flux due to VDF
    TYPE(atlas_field) :: fPFPLSL   ! liq+rain sedim flux
    TYPE(atlas_field) :: fPFPLSN   ! ice+snow sedim flux
    TYPE(atlas_field) :: fPFHPSL   ! Enthalpy flux for liq
    TYPE(atlas_field) :: fPFHPSN   ! Enthalpy flux for ice
    TYPE(atlas_field) :: fPCOVPTOT ! Precip fraction
    TYPE(atlas_field) :: fPRAINFRAC_TOPRFZ 
  CONTAINS
    PROCEDURE :: SETUP => CLOUDSC_GLOBAL_ATLAS_SETUP_BLOCK
  END TYPE CLOUDSC_GLOBAL_ATLAS_STATE_FIELDS

  TYPE CLOUDSC_GLOBAL_ATLAS_STATE
    ! Memory state containing raw fields annd tendencies for CLOUDSC dwarf
    !
    ! Note that the global state has an additional outermost block
    ! dimension allocated for each field variable.
    INTEGER(KIND=JPIM)                   :: NPROMA, KLEV    ! Grid points and vertical levels per block
    INTEGER(KIND=JPIM)                   :: NGPTOT, NBLOCKS ! Total number of grid points and blocks
    INTEGER(KIND=JPIM)                   :: KFLDX 
    LOGICAL                              :: LDSLPHY 
    LOGICAL                              :: LDMAINCALL      ! T if main call to cloudsc
    REAL(KIND=JPRB)                      :: PTSPHY          ! Physics timestep

  CONTAINS
    PROCEDURE :: LOAD => CLOUDSC_GLOBAL_ATLAS_STATE_LOAD
    PROCEDURE :: VALIDATE => CLOUDSC_GLOBAL_ATLAS_STATE_VALIDATE
  END TYPE CLOUDSC_GLOBAL_ATLAS_STATE

CONTAINS

  SUBROUTINE CLOUDSC_GLOBAL_ATLAS_SETUP_BLOCK(SELF, FSET)
    CLASS(CLOUDSC_GLOBAL_ATLAS_STATE_FIELDS), INTENT(INOUT) :: SELF
    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET

    SELF%fPLCRIT_AER = FSET%FIELD("PLCRIT_AER")
    SELF%fPICRIT_AER = FSET%FIELD("PICRIT_AER")
    SELF%fPRE_ICE = FSET%FIELD("PRE_ICE")
    SELF%fPCCN = FSET%FIELD("PCCN")
    SELF%fPNICE = FSET%FIELD("PNICE")
    SELF%fPT = FSET%FIELD("PT")
    SELF%fPQ = FSET%FIELD("PQ")
    SELF%fPVFA = FSET%FIELD("PVFA")
    SELF%fPVFL = FSET%FIELD("PVFL")
    SELF%fPVFI = FSET%FIELD("PVFI")
    SELF%fPDYNA = FSET%FIELD("PDYNA")
    SELF%fPDYNL = FSET%FIELD("PDYNL")
    SELF%fPDYNI = FSET%FIELD("PDYNI")
    SELF%fPHRSW = FSET%FIELD("PHRSW")
    SELF%fPHRLW = FSET%FIELD("PHRLW")
    SELF%fPVERVEL = FSET%FIELD("PVERVEL")
    SELF%fPAP = FSET%FIELD("PAP")
    SELF%fPLU = FSET%FIELD("PLU")
    SELF%fPLUDE = FSET%FIELD("PLUDE")
    SELF%fPSNDE = FSET%FIELD("PSNDE")
    SELF%fPMFU = FSET%FIELD("PMFU")
    SELF%fPMFD = FSET%FIELD("PMFD")
    SELF%fPA = FSET%FIELD("PA")
    SELF%fPSUPSAT = FSET%FIELD("PSUPSAT")
    SELF%fPLSM = FSET%FIELD("PLSM")
    SELF%fLDCUM = FSET%FIELD("LDCUM")
    SELF%fKTYPE = FSET%FIELD("KTYPE")
    SELF%fPAPH = FSET%FIELD("PAPH")
    SELF%fPEXTRA = FSET%FIELD("PEXTRA")
    SELF%fPCLV = FSET%FIELD("PCLV")

    SELF%fTENDENCY_CML = FSET%FIELD('TENDENCY_CML')
    SELF%fTENDENCY_TMP = FSET%FIELD('TENDENCY_TMP')
    SELF%fTENDENCY_LOC = FSET%FIELD('TENDENCY_LOC')

    SELF%fPFSQLF = FSET%FIELD("PFSQLF")
    SELF%fPFSQIF = FSET%FIELD("PFSQIF")
    SELF%fPFCQLNG = FSET%FIELD("PFCQLNG")
    SELF%fPFCQNNG = FSET%FIELD("PFCQNNG")
    SELF%fPFSQRF = FSET%FIELD("PFSQRF")
    SELF%fPFSQSF = FSET%FIELD("PFSQSF")
    SELF%fPFCQRNG = FSET%FIELD("PFCQRNG")
    SELF%fPFCQSNG = FSET%FIELD("PFCQSNG")
    SELF%fPFSQLTUR = FSET%FIELD("PFSQLTUR")
    SELF%fPFSQITUR = FSET%FIELD("PFSQITUR")
    SELF%fPFPLSL = FSET%FIELD("PFPLSL")
    SELF%fPFPLSN = FSET%FIELD("PFPLSN")
    SELF%fPFHPSL = FSET%FIELD("PFHPSL")
    SELF%fPFHPSN = FSET%FIELD("PFHPSN")
    SELF%fPCOVPTOT = FSET%FIELD("PCOVPTOT")
    SELF%fPRAINFRAC_TOPRFZ = FSET%FIELD("PRAINFRAC_TOPRFZ")
  END SUBROUTINE  

  SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK(SELF, FIELDS, IBLK)
    CLASS(CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK_VIEW), INTENT(INOUT) :: SELF
    CLASS(CLOUDSC_GLOBAL_ATLAS_STATE_FIELDS), INTENT(INOUT) :: FIELDS
    INTEGER, INTENT(IN) :: IBLK

    REAL(C_DOUBLE), POINTER :: TMP3D(:,:,:)

    ! NOTE the last six input variables need special treatment - different types
    CALL FIELDS%fPLCRIT_AER%DATA(SELF%PLCRIT_AER, IBLK)
    CALL FIELDS%fPICRIT_AER%DATA(SELF%PICRIT_AER, IBLK)
    CALL FIELDS%fPRE_ICE%DATA(SELF%PRE_ICE, IBLK)
    CALL FIELDS%fPCCN%DATA(SELF%PCCN, IBLK)
    CALL FIELDS%fPNICE%DATA(SELF%PNICE, IBLK)
    CALL FIELDS%fPT%DATA(SELF%PT, IBLK)
    CALL FIELDS%fPQ%DATA(SELF%PQ, IBLK)
    CALL FIELDS%fPVFA%DATA(SELF%PVFA, IBLK)
    CALL FIELDS%fPVFL%DATA(SELF%PVFL, IBLK)
    CALL FIELDS%fPVFI%DATA(SELF%PVFI, IBLK)
    CALL FIELDS%fPDYNA%DATA(SELF%PDYNA, IBLK)
    CALL FIELDS%fPDYNL%DATA(SELF%PDYNL, IBLK)
    CALL FIELDS%fPDYNI%DATA(SELF%PDYNI, IBLK)
    CALL FIELDS%fPHRSW%DATA(SELF%PHRSW, IBLK)
    CALL FIELDS%fPHRLW%DATA(SELF%PHRLW, IBLK)
    CALL FIELDS%fPVERVEL%DATA(SELF%PVERVEL, IBLK)
    CALL FIELDS%fPAP%DATA(SELF%PAP, IBLK)
    CALL FIELDS%fPLU%DATA(SELF%PLU, IBLK)
    CALL FIELDS%fPLUDE%DATA(SELF%PLUDE, IBLK)
    CALL FIELDS%fPSNDE%DATA(SELF%PSNDE, IBLK)
    CALL FIELDS%fPMFU%DATA(SELF%PMFU, IBLK)
    CALL FIELDS%fPMFD%DATA(SELF%PMFD, IBLK)
    CALL FIELDS%fPA%DATA(SELF%PA, IBLK)
    CALL FIELDS%fPSUPSAT%DATA(SELF%PSUPSAT, IBLK)
    CALL FIELDS%fPLSM%DATA(SELF%PLSM, IBLK)
    CALL FIELDS%fLDCUM%DATA(SELF%LDCUM, IBLK)
    CALL FIELDS%fKTYPE%DATA(SELF%KTYPE, IBLK)
    CALL FIELDS%fPAPH%DATA(SELF%PAPH, IBLK)
    CALL FIELDS%fPEXTRA%DATA(SELF%PEXTRA, IBLK)
    CALL FIELDS%fPCLV%DATA(SELF%PCLV, IBLK)

    CALL FIELDS%fTENDENCY_CML%DATA(TMP3D, IBLK)
    SELF%TENDENCY_CML%T   => TMP3D(:,:,1)
    SELF%TENDENCY_CML%A   => TMP3D(:,:,2)
    SELF%TENDENCY_CML%Q   => TMP3D(:,:,3)
    SELF%TENDENCY_CML%CLD => TMP3D(:,:,4:)
    CALL FIELDS%fTENDENCY_TMP%DATA(TMP3D, IBLK)
    SELF%TENDENCY_TMP%T   => TMP3D(:,:,1)
    SELF%TENDENCY_TMP%A   => TMP3D(:,:,2)
    SELF%TENDENCY_TMP%Q   => TMP3D(:,:,3)
    SELF%TENDENCY_TMP%CLD => TMP3D(:,:,4:)
    CALL FIELDS%fTENDENCY_LOC%DATA(TMP3D, IBLK)
    SELF%TENDENCY_LOC%T   => TMP3D(:,:,1)
    SELF%TENDENCY_LOC%A   => TMP3D(:,:,2)
    SELF%TENDENCY_LOC%Q   => TMP3D(:,:,3)
    SELF%TENDENCY_LOC%CLD => TMP3D(:,:,4:)

    CALL FIELDS%fPFSQLF%DATA(SELF%PFSQLF, IBLK)
    CALL FIELDS%fPFSQIF%DATA(SELF%PFSQIF, IBLK)
    CALL FIELDS%fPFCQLNG%DATA(SELF%PFCQLNG, IBLK)
    CALL FIELDS%fPFCQNNG%DATA(SELF%PFCQNNG, IBLK)
    CALL FIELDS%fPFSQRF%DATA(SELF%PFSQRF, IBLK)
    CALL FIELDS%fPFSQSF%DATA(SELF%PFSQSF, IBLK)
    CALL FIELDS%fPFCQRNG%DATA(SELF%PFCQRNG, IBLK)
    CALL FIELDS%fPFCQSNG%DATA(SELF%PFCQSNG, IBLK)
    CALL FIELDS%fPFSQLTUR%DATA(SELF%PFSQLTUR, IBLK)
    CALL FIELDS%fPFSQITUR%DATA(SELF%PFSQITUR, IBLK)
    CALL FIELDS%fPFPLSL%DATA(SELF%PFPLSL, IBLK)
    CALL FIELDS%fPFPLSN%DATA(SELF%PFPLSN, IBLK)
    CALL FIELDS%fPFHPSL%DATA(SELF%PFHPSL, IBLK)
    CALL FIELDS%fPFHPSN%DATA(SELF%PFHPSN, IBLK)
    CALL FIELDS%fPCOVPTOT%DATA(SELF%PCOVPTOT, IBLK)
    CALL FIELDS%fPRAINFRAC_TOPRFZ%DATA(SELF%PRAINFRAC_TOPRFZ, IBLK)
  END SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_BLOCK

  SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_LOAD(SELF, FSET, NPROMA, NGPTOT, NGPTOTG)
    ! Load reference input data via serialbox
    CLASS(CLOUDSC_GLOBAL_ATLAS_STATE) :: SELF
    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN) :: NPROMA, NGPTOT
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG

    TYPE(ATLAS_STRUCTUREDGRID) :: GRID
    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    INTEGER(KIND=JPIM) :: KLON, IVAR, B
    TYPE(VAR3D_PTR), DIMENSION(24) :: IN_VARS_3D_REAL64
    TYPE(VAR3D_PTR), DIMENSION(15) :: OUT_VARS_3D_REAL64
    REAL(C_DOUBLE), POINTER :: TMP3D(:,:,:)
    REAL(C_DOUBLE), POINTER :: TMP2D(:,:)
    TYPE(ATLAS_FIELD) :: FIELD

    CALL INPUT_INITIALIZE(NAME='input')

    CALL LOAD_SCALAR('KLON', KLON)
    CALL LOAD_SCALAR('KLEV', SELF%KLEV)
    CALL LOAD_SCALAR('KFLDX', SELF%KFLDX)

    GRID = ATLAS_REGULARLONLATGRID(NGPTOT, 1)
    FSPACE = ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS(GRID, LEVELS=SELF%KLEV, NPROMA=NPROMA, HALO=0)
    SELF%NBLOCKS = FSPACE%NBLKS()
    FSET = ATLAS_FIELDSET()

    DO IVAR = 1, SIZE(IN_VAR_NAMES) - 6 ! last six variables are special
        CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME=TRIM(IN_VAR_NAMES(IVAR)), KIND=ATLAS_REAL(JPRB)))
    ENDDO
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PLSM",   KIND=ATLAS_REAL(JPRB),    LEVELS=0))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="LDCUM",  KIND=ATLAS_LOGICAL(),     LEVELS=0))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="KTYPE",  KIND=ATLAS_INTEGER(JPIM), LEVELS=0))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PAPH",   KIND=ATLAS_REAL(JPRB),    LEVELS=SELF%KLEV+1))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PEXTRA", KIND=ATLAS_REAL(JPRB),    VARIABLES=MAX(1,SELF%KFLDX)))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PCLV",   KIND=ATLAS_REAL(JPRB),    VARIABLES=MAX(1,NCLV)))

    DO IVAR = 1, SIZE(IN_VAR_NAMES)
        CALL LOADVAR_ATLAS(FSET, TRIM(IN_VAR_NAMES(IVAR)), KLON, NGPTOTG)
    ENDDO

    FIELD = FSPACE%CREATE_FIELD(NAME='TENDENCY_CML', KIND=ATLAS_REAL(JPRB), VARIABLES=3+NCLV)
    CALL FSET%ADD(FIELD)
    FIELD = FSPACE%CREATE_FIELD(NAME='TENDENCY_TMP', KIND=ATLAS_REAL(JPRB), VARIABLES=3+NCLV)
    CALL FSET%ADD(FIELD)
    FIELD = FSPACE%CREATE_FIELD(NAME='TENDENCY_LOC', KIND=ATLAS_REAL(JPRB), VARIABLES=3+NCLV)
    CALL FSET%ADD(FIELD)

    ! The STATE_TYPE arrays are tricky, as the AOSOA layout needs to be expictly
    ! unrolled at every step, and we rely on dirty hackery to do this.
    CALL LOADSTATE_ATLAS(FSET, 'TENDENCY_CML', KLON, NGPTOTG)
    CALL LOADSTATE_ATLAS(FSET, 'TENDENCY_TMP', KLON, NGPTOTG)
    ! Output fields are simply allocated and zero'd
    DO IVAR = 1, SIZE(OUT_VAR_NAMES) - 2
        CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME=TRIM(OUT_VAR_NAMES(IVAR)), KIND=ATLAS_REAL(JPRB), LEVELS=SELF%KLEV+1))
    ENDDO
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PCOVPTOT", KIND=ATLAS_REAL(JPRB)))
    CALL FSET%ADD(FSPACE%CREATE_FIELD(NAME="PRAINFRAC_TOPRFZ", KIND=ATLAS_REAL(JPRB), LEVELS=0))

    DO IVAR = 1, SIZE(OUT_VAR_NAMES) - 1
        FIELD = FSET%FIELD(TRIM(OUT_VAR_NAMES(IVAR)))
        CALL FIELD%DATA(TMP3D)
        !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B) schedule(runtime)
        DO B=1, SELF%NBLOCKS
           TMP3D(:,:,B) = 0.0_JPRB
        END DO
        !$omp end parallel do
    ENDDO
    ! DEBUG
    !FIELD = FSET%FIELD("PAP")
    !call field%data(tmp3d)
    !print *, MINVAL(tmp3d), MAXVAL(tmp3d)

    FIELD = FSET%FIELD("PRAINFRAC_TOPRFZ")
    CALL FIELD%DATA(TMP2D)
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B) schedule(runtime)
    DO B=1, SELF%NBLOCKS
       TMP2D(:,B) = 0.0_JPRB
    END DO
    !$OMP END PARALLEL DO

    ! Initialize global parameters from the input file
    CALL LOAD_SCALAR('PTSPHY', SELF%PTSPHY)
    CALL LOAD_SCALAR('LDSLPHY', SELF%LDSLPHY)
    CALL LOAD_SCALAR('LDMAINCALL', SELF%LDMAINCALL)
    CALL YOMCST_LOAD_PARAMETERS()
    CALL YOETHF_LOAD_PARAMETERS()
    CALL YRECLDP_LOAD_PARAMETERS()
    CALL YREPHLI_LOAD_PARAMETERS()

    CALL INPUT_FINALIZE()
  END SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_LOAD

  SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_VALIDATE(SELF, FSET, NGPTOT, NGPTOTG)
    ! Validate the correctness of output against reference data
    CLASS(CLOUDSC_GLOBAL_ATLAS_STATE) :: SELF
    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    INTEGER(KIND=JPIM), INTENT(IN) :: NGPTOT
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG

    INTEGER(KIND=JPIM) :: KLON, IVAR

    CALL INPUT_INITIALIZE(NAME='reference')
    CALL LOAD_SCALAR('KLON', KLON)
    print *, "KLON = ", KLON
    CALL INPUT_FINALIZE()

    ! Write variable validation header
    IF (IRANK == 0) THEN
      print '(1X,A20,1X,A3,5(1X,A20))', &
           & 'Variable','Dim', 'MinValue','MaxValue','AbsMaxErr','AvgAbsErr/GP','MaxRelErr-%'
    END IF


    ! Actual variable validation
    CALL VALIDATEVAR_ATLAS(FSET, 'PLUDE', KLON, NGPTOTG)
    DO IVAR = 1, SIZE(OUT_VAR_NAMES)
        CALL VALIDATEVAR_ATLAS(FSET, OUT_VAR_NAMES(IVAR), KLON, NGPTOTG)
    ENDDO
    CALL VALIDATESTATE_ATLAS(FSET, 'TENDENCY_LOC', KLON, NGPTOTG)

  END SUBROUTINE CLOUDSC_GLOBAL_ATLAS_STATE_VALIDATE

END MODULE CLOUDSC_GLOBAL_ATLAS_STATE_MOD
