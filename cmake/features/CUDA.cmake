# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_CUDA )

  enable_language( CUDA )

  if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    set( CLOUDSC_CUDA_FLAGS "-lineinfo -maxrregcount=128" )
  else()
    foreach(cuda_arch IN LISTS CMAKE_CUDA_ARCHITECTURES)
       if(cuda_arch LESS 70)
           ecbuild_warn("Specified CUDA architecture ${cuda_arch} is unsupported and will not compile. 
                            Minimum required >= 70. I am attempting to use native CUDA architecture instead.")
           if( ${CMAKE_CUDA_ARCHITECTURES_NATIVE} LESS 70 )
             ecbuild_error("Native CUDA architecture ${CMAKE_CUDA_ARCHITECTURES_NATIVE} is also unsupported.
                                   Please provide proper arch file or CMAKE_CUDA_ARCHITECTURES at build time or deactivate CUDA variant.")
           else()
             set( CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_NATIVE})
             ecbuild_warn("CUDA architecture is now set to ${CMAKE_CUDA_ARCHITECTURES_NATIVE} (native version)." )
             break()
           endif()
       endif()
    endforeach()
    set( CLOUDSC_CUDA_FLAGS "-lineinfo -maxrregcount=128 -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}" )
  endif()

  if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    set( CLOUDSC_CUDA_OPT_FLAGS "-O0 -g -G" )
  else()
    set( CLOUDSC_CUDA_OPT_FLAGS  "--ptxas-options=-O3 -O3 -use_fast_math" )
  endif()

endif()
