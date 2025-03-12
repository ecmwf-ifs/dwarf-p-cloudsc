# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/common/module" )
set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/common/include" )

# Available precisions
set( CLOUDSC_PRECISIONS )
if( HAVE_DOUBLE_PRECISION )
  list( APPEND CLOUDSC_PRECISIONS dp )
endif()
if( HAVE_SINGLE_PRECISION )
  list( APPEND CLOUDSC_PRECISIONS sp )
endif()

foreach( prec IN LISTS CLOUDSC_PRECISIONS )

  if( "${prec}" STREQUAL "sp" )
    list( APPEND CLOUDSC_DEFINITIONS SINGLE )
  endif()

  # Make sure to build the common libs first
  add_subdirectory( common common.${prec} )

  # A few variants only support double precision
  if( "${prec}" STREQUAL "dp" )
    add_subdirectory( prototype1 cloudsc_prototype1.${prec} )
    add_subdirectory( cloudsc_c cloudsc_c.${prec} )
    if( HAVE_PYIFACE )
      add_subdirectory( cloudsc_pyiface cloudsc_pyiface.${prec} )
    endif()
    if( HAVE_PYTHON_F2PY )
      add_subdirectory( cloudsc_python cloudsc_python.${prec} )
    endif()
  endif()

  # Add all remaining subdirectories with precision-specific binary directories
  add_subdirectory( cloudsc_fortran cloudsc_fortran.${prec} )
  add_subdirectory( cloudsc_gpu cloudsc_gpu.${prec} )

  if( HAVE_CUDA )
    add_subdirectory( cloudsc_cuda cloudsc_cuda.${prec} )
  endif()
  if( HAVE_ATLAS )
    add_subdirectory( cloudsc_fortran_atlas cloudsc_fortran_atlas.${prec} )
  endif()
  if( HAVE_HIP )
    add_subdirectory( cloudsc_hip cloudsc_hip.${prec} )
  endif()
  if( HAVE_LOKI )
    add_subdirectory( cloudsc_loki cloudsc_loki.${prec} )
  endif()
  if( HAVE_SYCL )
    add_subdirectory( cloudsc_sycl cloudsc_sycl.${prec} )
  endif()

endforeach()
