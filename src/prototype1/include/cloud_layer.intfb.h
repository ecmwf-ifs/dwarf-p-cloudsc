! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTERFACE
SUBROUTINE CLOUD_LAYER(&
 & KDIM, LDSLPHY, LDMAINCALL, PAUX, state, tendency_cml, tendency_tmp, tendency_dyn, tendency_vdf, PRAD,&
 & PSAVTENDSAT, PSURF, LLKEYS,&
 & AUXL, FLUX, PDIAG,&
 & tendency_loc) 
USE PARKIND1 ,ONLY : JPIM, JPRB
USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE,&
 & FLUX_TYPE, AUX_DIAG_LOCAL_TYPE, AUX_RAD_TYPE, KEYS_LOCAL_TYPE,&
 & SURF_AND_MORE_TYPE, AUX_DIAG_TYPE 
TYPE (DIMENSION_TYPE) , INTENT (IN) :: KDIM
LOGICAL , INTENT (IN) :: LDSLPHY
LOGICAL , INTENT (IN) :: LDMAINCALL
TYPE (AUX_TYPE) , INTENT (IN) :: PAUX
TYPE (STATE_TYPE) , INTENT (IN) :: state
TYPE (STATE_TYPE) , INTENT (IN) :: tendency_cml
TYPE (STATE_TYPE) , INTENT (IN) :: tendency_tmp
TYPE (STATE_TYPE) , INTENT (IN) :: tendency_dyn
TYPE (STATE_TYPE) , INTENT (IN) :: tendency_vdf
TYPE (AUX_RAD_TYPE) , INTENT (IN) :: PRAD
REAL(KIND=JPRB) , INTENT (IN) :: PSAVTENDSAT(KDIM%KLON,KDIM%KLEV)
TYPE (SURF_AND_MORE_TYPE) , INTENT(INOUT) :: PSURF
TYPE (KEYS_LOCAL_TYPE) , INTENT (IN) :: LLKEYS
TYPE (AUX_DIAG_LOCAL_TYPE) , INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE) , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE) , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE) , INTENT (OUT) :: tendency_loc
END SUBROUTINE CLOUD_LAYER
END INTERFACE
