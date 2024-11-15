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
set( ENABLE_USE_STMT_FUNC ON CACHE STRING "" )

####################################################################
# OpenMP FLAGS
####################################################################

# set( OpenMP_C_FLAGS           "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" CACHE STRING "" )
# set( OpenMP_CXX_FLAGS         "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" CACHE STRING "" )
# set( OpenMP_Fortran_FLAGS     "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" CACHE STRING "" )

# set( OpenMP_C_FLAGS           "-fopenmp" CACHE STRING "")
# set( OpenMP_CXX_FLAGS         "-fopenmp" CACHE STRING "")
# set( OpenMP_Fortran_FLAGS     "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn --offload-arch=gfx90a" CACHE STRING "")

# set( OpenMP_C_LIB_NAMES       "craymp" CACHE STRING "" )
# set( OpenMP_CXX_LIB_NAMES     "craymp" CACHE STRING "" )
# set( OpenMP_Fortran_LIB_NAMES "craymp" CACHE STRING "" )
# set( OpenMP_craymp_LIBRARY    "/opt/cray/pe/cce/16.0.1/cce/x86_64/lib/libcraymp.so" CACHE STRING "" )

####################################################################
# OpenACC FLAGS
####################################################################

# set( OpenACC_C_FLAGS "-hacc" CACHE STRING "" )
# set( OpenACC_CXX_FLAGS "-hacc" CACHE STRING "" )
# set( OpenACC_Fortran_FLAGS "-hacc" CACHE STRING "" )

####################################################################
# Compiler FLAGS
####################################################################

# set(ECBUILD_C_FLAGS "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn --offload-arch=gfx90a") #  -march=gfx90a")
# set(ECBUILD_CXX_FLAGS "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn --offload-arch=gfx90a") #  -march=gfx90a")

# General Flags (add to default)
# set(ECBUILD_Fortran_FLAGS "-v -g -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn --offload-arch=gfx90a") #  -march=gfx90a")
set(ECBUILD_Fortran_FLAGS "-fopenmp --offload-arch=gfx90a -lFortranRuntimeHostDevice") #  -module-dir '/scratch/project_465000527/staneker/cloudsc-old-PRs/module_files'")
# set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fhermetic-module-files")
# set(ECBUILD_Fortran_FLAGS "-hcontiguous")
# set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
# set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")

# set(ECBUILD_Fortran_FLAGS_BIT "-O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")

if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
  set(CMAKE_HIP_ARCHITECTURES gfx90a)
endif()

# list(APPEND CMAKE_PREFIX_PATH "/home/rip/Qt/5.12.1/gcc_64")

# select OpenMP pragma to be used
set( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL OFF CACHE BOOL "" )
