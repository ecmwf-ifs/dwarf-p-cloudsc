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

set(ECBUILD_Fortran_FLAGS "-g")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -qopenmp-threadprivate compat")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume byterecl")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -convert big_endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -traceback")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -align array64byte")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -warn nounused,nouncalled")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -march=core-avx2")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-functions")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-limit=1500")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Winline")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -no-fma")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume realloc_lhs")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fp-model precise")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ftz")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fp-speculation=safe")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fast-transcendentals")

####################################################################
# LINK FLAGS
####################################################################

set( ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr -Ktrap=fp" )
set( ECBUILD_MODULE_LINKER_FLAGS "-Wl,--eh-frame-hdr -Ktrap=fp -Wl,-Map,loadmap" )
set( ECBUILD_EXE_LINKER_FLAGS    "-Wl,--eh-frame-hdr -Ktrap=fp -Wl,-Map,loadmap -Wl,--as-needed" )
set( ECBUILD_CXX_IMPLICIT_LINK_LIBRARIES "${LIBCRAY_CXX_RTS}" CACHE STRING "" )
