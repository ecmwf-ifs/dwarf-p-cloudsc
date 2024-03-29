! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMPHY2

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*
!     ------------------------------------------------------------------
!     CONSTANTES DEFINISSANT LE CONTEXTE DE L'EXPERIENCE POUR LA
!     PHYSIQUE : DISCRETISATION TEMPORELLE, STRUCTURE VERTICALE,
!     SORTIES.
!       TSPHY     : PAS DE TEMPS DE LA PHYSIQUE.
!                 : PHYSICS TIME STEP.
!       XMUCVPP   : FACTEUR DE CVPP POUR "ANTI-FIBRIL.". 
!                 : "ANTI-FIBRIL." FACTOR FOR CVPP.
!       XMULAF    : FACTEUR "ANTI-FIBRIL." (1. RECOMMANDE SI ACTIF).
!                 : "ANTI-FIBRIL." FACTOR (1. RECOMMENDED IF ACTIVATED).
!       XDAMP     : DAMPING FACTOR USED IN THE NEW TREATMENT OF
!                   SHALLOW CONVECTION.
!                   IF ZERO => OLD TREATMENT
!                   IF /= 0 => XMUCVPP SHOULD BE ZERO
!       LMULAF    : CONTROL DE l'ANTI-FIBRIL. SUR LA VERTICALE
!                 : CONTROL "ANTI-FIBRIL." ON VERTICAL
!       HCLP      : HAUTEUR MOYENNE DE LA CLP ( EN GENERAL 1500M ).
!                 : MEAN PBL DEPTH ( IN GENERAL 1500M ).
!       HTCLS     : HAUTEUR METEO POUR T ET Q ( EN GENERAL 2M ).
!                 : SCREEN HEIGHT FOR T AND Q ( IN GENERAL 2M ).
!       HVCLS     : HAUTEUR METEO POUR U ET V ( EN GENERAL 10M ).
!                 : MEASURING HEIGHT FOR U AND V (IN GENERAL 10M ).

!       HTSHM     : HAUTEUR DE TRANSITION EN "S" NUAGES HAUTS/MOYENS.
!                 : TRANSITION HEIGHT IN "S" COORDINATE H/M CLOUDS.
!       HTSML     : HAUTEUR DE TRANSITION EN "S" NUAGES MOYENS/BAS.
!                 : TRANSITION HEIGHT IN "S" COORDINATE M/L CLOUDS.
!       LWMOCLOUD : ACTIVATION DU DIAGNOSTIC DE L'HAUTEUR NUAGES SELON L'OMM.
!                 : ACTIVATE CLOUD HEIGHT DIAGNOSTICS ACCORDING TO WMO.
!       HWMOHIGH  : HAUTEUR DE NUAGES HAUTS SELON L'OMM.
!                 : HEIGHT OF HIGH CLOUDS ACCORDING TO WMO.
!       HWMOLOW   : HAUTEUR DE NUAGES BAS SELON L'OMM.
!                 : HEIGHT OF LOW CLOUDS ACCORDING TO WMO.
!       NTSHM     : INDICE DU NIVEAU DE TRANSITION NUAGES HAUTS/MOYENS.
!                 : TRANSITION LEVEL BETWEEN HIGH/MEDIUM CLOUDS.
!       NTSML     : INDICE DU NIVEAU DE TRANSITION NUAGES MOYENS/BAS.
!                 : TRANSITION LEVEL BETWEEN MEDIUM/LOW CLOUDS.
!       LRAFTUR   : ACTIVATION DU DIAGNOSTIC DES RAFALES TURBULENTES
!                   ACTIVATE DIAGNOSTIC OF TURBULENT GUSTS
!       LRAFTKE   : ACTIVATION DU DIAGNOSTIC DES RAFALES TURBULENTES AVEC LA TKE
!                   ACTIVATE DIAGNOSTIC OF TURBULENT GUSTS FROM TKE
!       GZ0RAF    : Z0 FOIS G UTILISE POUR LE CALCUL DES RAFALES TURBULENTES
!                   Z0 TIMES G USED TO COMPUTE TURBULENT GUSTS
!       FACRAF    : COEFFICIENT DE CALCUL DES RAFALES TURBULENTES
!                   COEFFICIENT FOR THE COMPUTATION OF TURBULENT GUSTS
!       HTKERAF   : HAUTEUR DE CALCUL DES RAFALES TURBULENTES
!                   HEIGHT FOR THE COMPUTATION OF TURBULENT GUSTS

INTEGER(KIND=JPIM) :: NTSHM
INTEGER(KIND=JPIM) :: NTSML
REAL(KIND=JPRB) :: TSPHY
REAL(KIND=JPRB) :: XMUCVPP
REAL(KIND=JPRB) :: XMULAF
REAL(KIND=JPRB) :: XDAMP
REAL(KIND=JPRB) :: HCLP
REAL(KIND=JPRB) :: HTCLS
REAL(KIND=JPRB) :: HVCLS
REAL(KIND=JPRB) :: HTSHM
REAL(KIND=JPRB) :: HTSML
REAL(KIND=JPRB) :: HWMOHIGH
REAL(KIND=JPRB) :: HWMOLOW
REAL(KIND=JPRB) :: GZ0RAF
REAL(KIND=JPRB) :: FACRAF
REAL(KIND=JPRB) :: HTKERAF
LOGICAL :: LRAFTUR
LOGICAL :: LMULAF
LOGICAL :: LWMOCLOUD
LOGICAL :: LRAFTKE
!     ------------------------------------------------------------------
END MODULE YOMPHY2
