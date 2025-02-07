# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Macro to create executable and accompanying tests

function( cloudsc_add_fortran_exe )

    set( options "" )
    set( single_value_args TARGET DRIVER COMMON_LIB )
    set( multi_value_args DEFINITIONS )

    cmake_parse_arguments( _CL "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

    if( _CL_UNPARSED_ARGUMENTS )
        ecbuild_critical( "Unknown keywords given to cloudsc_add_fortran_exe(): \"${_CL_UNPARSED_ARGUMENTS}\"")
    endif()

    ecbuild_add_executable(
        TARGET ${_CL_TARGET}-${prec}
        SOURCES
            dwarf_cloudsc.F90
            ${_CL_DRIVER}
            cloudsc.F90
        LIBS
            cloudsc-common-${prec}-${_CL_COMMON_LIB}
        DEFINITIONS
            ${_CL_DEFINITIONS}
    )

    set_target_properties( ${_CL_TARGET}-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${_CL_TARGET}"
    )

    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-serial
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-omp
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-mpi
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-mpi-omp
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )
endfunction()

### Definition of executable variants

# Default Fortran driver
cloudsc_add_fortran_exe(
    TARGET dwarf-cloudsc-fortran
    DRIVER cloudsc_driver_mod.F90
    COMMON_LIB lib
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)

# Default Fortran driver using FIELD_API
if( HAVE_FIELD_API )
    cloudsc_add_fortran_exe(
        TARGET dwarf-cloudsc-fortran-field
        DRIVER cloudsc_driver_field_mod.F90
        COMMON_LIB field_api
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            CLOUDSC_FIELD
    )
endif()
