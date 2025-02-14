# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Helper macro to create a target and associated tests

macro( cloudsc_add_hip_exe TARGET DRIVER KERNEL )

    ecbuild_add_library(
        TARGET ${TARGET}-${prec}-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/dtype.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
            ${DRIVER}
            ${KERNEL}
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
            $<${HAVE_OMP}:OpenMP::OpenMP_CXX>
            hip::device
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            $<$<STREQUAL:"${prec}","sp">:SINGLE>
    )

    set_source_files_properties( ${DRIVER} ${KERNEL}
        PROPERTIES LANGUAGE HIP )

    if( CMAKE_CXX_COMPILER_ID MATCHES "Cray" AND HAVE_OMP )
        # Workaround for Linker issue with Cray compiler, see
        # https://gitlab.kitware.com/cmake/cmake/-/issues/24402
        target_link_options( ${TARGET}-${prec}-lib PUBLIC "-fopenmp" )
    endif()

    ecbuild_add_executable(
        TARGET      ${TARGET}-${prec}
        SOURCES     dwarf_cloudsc.cpp
        LIBS
            ${TARGET}-${prec}-lib
    )


    ecbuild_add_test(
        TARGET      ${TARGET}-${prec}-serial
        COMMAND     ${TARGET}-${prec}
        ARGS        1 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

endmacro()

### Kernel variants

# SCC
cloudsc_add_hip_exe(
    dwarf-cloudsc-hip
    "cloudsc/cloudsc_driver.h;cloudsc/cloudsc_driver.cpp"
    "cloudsc/cloudsc_c.h;cloudsc/cloudsc_c.cpp"
)

# SCC-HOIST
cloudsc_add_hip_exe(
    dwarf-cloudsc-hip-hoist
    "cloudsc/cloudsc_driver_hoist.h;cloudsc/cloudsc_driver_hoist.cpp"
    "cloudsc/cloudsc_c_hoist.h;cloudsc/cloudsc_c_hoist.cpp"
)

# SCC-K-CACHING
cloudsc_add_hip_exe(
    dwarf-cloudsc-hip-k-caching
    "cloudsc/cloudsc_driver.h;cloudsc/cloudsc_driver.cpp"
    "cloudsc/cloudsc_c_k_caching.h;cloudsc/cloudsc_c_k_caching.cpp"
)
