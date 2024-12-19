# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Disable MPI in Github runner with NVHPC
set( ENABLE_MPI OFF CACHE STRING "" )

####################################################################
# OpenMP FLAGS
####################################################################

# Note: OpenMP_Fortran_FLAGS gets overwritten by the FindOpenMP module
# unless its stored as a cache variable
set( OpenMP_Fortran_FLAGS   "-mp -mp=bind,allcores,numa" CACHE STRING "" )
set( OpenMP_C_FLAGS         "-mp -mp=bind,allcores,numa" CACHE STRING "" )

####################################################################
# OpenACC FLAGS
####################################################################

# NB: We have to add `-mp` again to avoid undefined symbols during linking
# (smells like an Nvidia bug)
set( OpenACC_Fortran_FLAGS "-acc=gpu -mp" CACHE STRING "" )
# Enable this to get more detailed compiler output
# set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Minfo" )

####################################################################
# COMMON FLAGS
####################################################################

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

set(ECBUILD_Fortran_FLAGS "-fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mframe")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz")

set( ECBUILD_Fortran_FLAGS_BIT "-O2 -gopt" )

set( ECBUILD_C_FLAGS "-O2 -gopt -traceback" )

set( ECBUILD_CXX_FLAGS "-O2 -gopt" )
