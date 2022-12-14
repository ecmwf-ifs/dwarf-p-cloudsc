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

set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if( ENABLE_MPI )
set(CMAKE_C_COMPILER "mpifcc")
set(CMAKE_CXX_COMPILER "mpiFCC")
set(CMAKE_Fortran_COMPILER "mpifrt")
else()
set(CMAKE_C_COMPILER "fcc")
set(CMAKE_CXX_COMPILER "FCC")
set(CMAKE_Fortran_COMPILER "frt")
endif()

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS   "-Kopenmp -Nfjomplib" CACHE STRING "" )
set( OpenMP_C_LIB_NAMES   "" CACHE STRING "" )
set( OpenMP_CXX_FLAGS   "-Kopenmp -Nfjomplib" CACHE STRING "" )
set( OpenMP_CXX_LIB_NAMES   "" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-Kopenmp -Nfjomplib" CACHE STRING "" )

####################################################################
# COMMON FLAGS
####################################################################

foreach(LANG IN ITEMS C CXX Fortran LINKER)
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -fpic")
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -Klto")
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -SSL2")
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -Kopenmp")
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -Nfjomplib")
endforeach()

set(CMAKE_EXE_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")

set( ECBUILD_Fortran_FLAGS_BIT "-Kfast -O3 -KA64FX -KSVE -KARMV8_3_A -Ksimd=2 -Kassume=notime_saving_compilation -DNDEBUG" )

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
