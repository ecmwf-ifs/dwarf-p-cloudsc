# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI OFF )
set( ENABLE_ACC OFF )

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS   "-homp -hnoacc" CACHE STRING "" )
set( OpenMP_C_LIB_NAMES   "sci_cray_mpi_mp;sci_cray_mp;craymp;atomic" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-homp -hnoacc" CACHE STRING "" )


####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")

set(ECBUILD_Fortran_FLAGS_BIT "-O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")

####################################################################
# LINK FLAGS
####################################################################

if( EXISTS "$ENV{CC_X86_64}/lib/x86-64/libcray-c++-rts.so" )
  set( LIBCRAY_CXX_RTS "$ENV{CC_X86_64}/lib/x86-64/libcray-c++-rts.so" )
elseif( EXISTS "$ENV{CC_X86_64}/lib/libcray-c++-rts.so" )
  set( LIBCRAY_CXX_RTS "$ENV{CC_X86_64}/lib/libcray-c++-rts.so" )
endif()

set( ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr -Ktrap=fp" )
set( ECBUILD_MODULE_LINKER_FLAGS "-Wl,--eh-frame-hdr -Ktrap=fp -Wl,-Map,loadmap" )
set( ECBUILD_EXE_LINKER_FLAGS    "-Wl,--eh-frame-hdr -Ktrap=fp -Wl,-Map,loadmap -Wl,--as-needed" )
set( ECBUILD_CXX_IMPLICIT_LINK_LIBRARIES "${LIBCRAY_CXX_RTS}" CACHE STRING "" )

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
