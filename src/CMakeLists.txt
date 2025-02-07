# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/common/module" )
set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/common/include" )

add_subdirectory(prototype1)
add_subdirectory(common)
add_subdirectory(cloudsc_fortran)
add_subdirectory(cloudsc_fortran_atlas)
add_subdirectory(cloudsc_pyiface)
add_subdirectory(cloudsc_python)
add_subdirectory(cloudsc_c)
if( HAVE_CUDA )
add_subdirectory(cloudsc_cuda)
endif()
add_subdirectory(cloudsc_hip)
add_subdirectory(cloudsc_sycl)
add_subdirectory(cloudsc_gpu)
if( HAVE_LOKI )
add_subdirectory(cloudsc_loki)
endif()
