# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)

set(ECBUILD_Fortran_FLAGS "-fpe0")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -convert big_endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume noold_maxminloc")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -diag-disable=10441")

set(ECBUILD_Fortran_FLAGS_BIT "-g -O2 -traceback")
set(ECBUILD_C_FLAGS_BIT "-g -O2 -diag-disable=10441")
set(ECBUILD_CXX_FLAGS_BIT "-g -O2 -diag-disable=10441")
