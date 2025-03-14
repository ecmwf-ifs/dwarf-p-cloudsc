! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM DWARF_CLOUDSC

USE PARKIND1, ONLY: JPIM, JPIB
USE CLOUDSC_MPI_MOD, ONLY: CLOUDSC_MPI_INIT, CLOUDSC_MPI_END, NUMPROC, IRANK

USE YOECLDP  , ONLY : YRECLDP
USE YOMCST   , ONLY : YRCST 
USE YOETHF   , ONLY : YRTHF

#ifdef _OPENMP
USE OMP_LIB
#endif

#ifdef CLOUDSC_FIELD
USE CLOUDSC_FIELD_STATE_MOD, ONLY: CLOUDSC_FIELD_STATE
USE CLOUDSC_DRIVER_FIELD_MOD, ONLY: CLOUDSC_DRIVER_FIELD
#else
USE CLOUDSC_GLOBAL_STATE_MOD, ONLY: CLOUDSC_GLOBAL_STATE
USE CLOUDSC_DRIVER_MOD, ONLY: CLOUDSC_DRIVER
#endif

IMPLICIT NONE

CHARACTER(LEN=20) :: CLARG
INTEGER(KIND=JPIM) :: IARGS, LENARG, JARG, I
CHARACTER(LEN=20) :: PACKED_STORAGE
LOGICAL :: USE_PACKED

INTEGER(KIND=JPIM) :: NUMOMP   = 1     ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOTG  = 16384 ! Number of grid points (as read from command line)
INTEGER(KIND=JPIM) :: NPROMA   = 32    ! NPROMA blocking factor (currently active)
INTEGER(KIND=JPIM) :: NGPTOT           ! Local number of grid points

#ifdef CLOUDSC_FIELD
TYPE(CLOUDSC_FIELD_STATE) :: GLOBAL_STATE
#else
TYPE(CLOUDSC_GLOBAL_STATE) :: GLOBAL_STATE
#endif

IARGS = COMMAND_ARGUMENT_COUNT()

! Get the number of OpenMP threads to use for the benchmark
if (IARGS >= 1) then
  CALL GET_COMMAND_ARGUMENT(1, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NUMOMP
  if (NUMOMP <= 0) then
#ifdef _OPENMP
    NUMOMP = OMP_GET_MAX_THREADS()
#else
    ! if arg is 0 or negative, and OpenMP disabled; defaults to 1
    NUMOMP = 1
#endif
  end if
end if

! Initialize MPI environment
CALL CLOUDSC_MPI_INIT(NUMOMP)

! Get total number of grid points (NGPTOT) with which to run the benchmark
IF (IARGS >= 2) THEN
  CALL GET_COMMAND_ARGUMENT(2, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NGPTOTG
END IF

! Determine local number of grid points
NGPTOT = (NGPTOTG - 1) / NUMPROC + 1
if (IRANK == NUMPROC - 1) then
  NGPTOT = NGPTOTG - (NUMPROC - 1) * NGPTOT
end if

! Get the block size (NPROMA) for which to run the benchmark  
IF (IARGS >= 3) THEN
  CALL GET_COMMAND_ARGUMENT(3, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NPROMA
ENDIF

#ifdef CLOUDSC_FIELD
! FIXME: this needs to be fixed
! USE_PACKED = .TRUE.
CALL GET_ENVIRONMENT_VARIABLE('CLOUDSC_PACKED_STORAGE', PACKED_STORAGE)
USE_PACKED = TRIM(PACKED_STORAGE) == 'ON' .OR. TRIM(PACKED_STORAGE) == '1'
! Create a global memory state using FIELD objects from serialized input data
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG, USE_PACKED=USE_PACKED)
#else
! Create a global memory state from serialized input data
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)
#endif

! Call the driver to perform the parallel loop over our kernel
#ifdef CLOUDSC_FIELD
CALL CLOUDSC_DRIVER_FIELD( NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, NGPTOTG, &
     & GLOBAL_STATE%KFLDX, GLOBAL_STATE%PTSPHY, &
     & GLOBAL_STATE%AUX, GLOBAL_STATE%FLUX, &
     & GLOBAL_STATE%TENDENCY_TMP, GLOBAL_STATE%TENDENCY_LOC, &
     & YRCST, YRTHF, YRECLDP)
#else
CALL CLOUDSC_DRIVER( NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, NGPTOTG, &
     & GLOBAL_STATE%KFLDX, GLOBAL_STATE%PTSPHY, &
     & GLOBAL_STATE%PT, GLOBAL_STATE%PQ, &
     & GLOBAL_STATE%TENDENCY_CML, GLOBAL_STATE%TENDENCY_TMP, GLOBAL_STATE%TENDENCY_LOC, &
     & GLOBAL_STATE%PVFA,    GLOBAL_STATE%PVFL,  GLOBAL_STATE%PVFI, &
     & GLOBAL_STATE%PDYNA,   GLOBAL_STATE%PDYNL, GLOBAL_STATE%PDYNI, &
     & GLOBAL_STATE%PHRSW,   GLOBAL_STATE%PHRLW, &
     & GLOBAL_STATE%PVERVEL, GLOBAL_STATE%PAP,   GLOBAL_STATE%PAPH, &
     & GLOBAL_STATE%PLSM,    GLOBAL_STATE%LDCUM, GLOBAL_STATE%KTYPE, &
     & GLOBAL_STATE%PLU,     GLOBAL_STATE%PLUDE, GLOBAL_STATE%PSNDE, &
     & GLOBAL_STATE%PMFU,    GLOBAL_STATE%PMFD, &
     & GLOBAL_STATE%PA,      GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,&
     & GLOBAL_STATE%PLCRIT_AER, GLOBAL_STATE%PICRIT_AER, GLOBAL_STATE%PRE_ICE, &
     & GLOBAL_STATE%PCCN,     GLOBAL_STATE%PNICE,&
     & GLOBAL_STATE%PCOVPTOT, GLOBAL_STATE%PRAINFRAC_TOPRFZ, &
     & GLOBAL_STATE%PFSQLF,   GLOBAL_STATE%PFSQIF ,  GLOBAL_STATE%PFCQNNG,  GLOBAL_STATE%PFCQLNG, &
     & GLOBAL_STATE%PFSQRF,   GLOBAL_STATE%PFSQSF ,  GLOBAL_STATE%PFCQRNG,  GLOBAL_STATE%PFCQSNG, &
     & GLOBAL_STATE%PFSQLTUR, GLOBAL_STATE%PFSQITUR, &
     & GLOBAL_STATE%PFPLSL,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL,   GLOBAL_STATE%PFHPSN, &
     & YRCST, YRTHF, YRECLDP)
#endif

! Validate the output against serialized reference data
CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
#ifdef CLOUDSC_FIELD
CALL GLOBAL_STATE%FINALIZE()
#endif

! Tear down MPI environment
CALL CLOUDSC_MPI_END()

END PROGRAM DWARF_CLOUDSC
