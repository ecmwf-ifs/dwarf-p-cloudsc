! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCST_initialise

use yomcst, only : rg, rcpd, rd, retv, rlvtt, rlstt, rlmlt, rtt, rv, &
       & rg_h, rcpd_h, rd_h, retv_h, rlvtt_h, rlstt_h, rlmlt_h, rtt_h, rv_h

IMPLICIT NONE


CONTAINS

SUBROUTINE YOMCST_LOAD_PARAMETERS()
USE FILE_IO_MOD, ONLY : LOAD_SCALAR_real
    CALL LOAD_SCALAR_real('RG', RG_h)
    CALL LOAD_SCALAR_real('RD', RD_h)
    CALL LOAD_SCALAR_real('RCPD', RCPD_h)
    CALL LOAD_SCALAR_real('RETV', RETV_h)
    CALL LOAD_SCALAR_real('RLVTT', RLVTT_h)
    CALL LOAD_SCALAR_real('RLSTT', RLSTT_h)
    CALL LOAD_SCALAR_real('RLMLT', RLMLT_h)
    CALL LOAD_SCALAR_real('RTT', RTT_h)
    CALL LOAD_SCALAR_real('RV', RV_h)

    rg=rg_h; rd=rd_h; rcpd=rcpd_h; retv=retv_h; rlvtt=rlvtt_h; rlstt=rlstt_h; rlmlt=rlmlt_h; rtt=rtt_h; rv=rv_h
  END SUBROUTINE YOMCST_LOAD_PARAMETERS

END MODULE YOMCST_initialise
