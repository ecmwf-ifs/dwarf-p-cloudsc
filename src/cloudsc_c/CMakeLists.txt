# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_C
    DESCRIPTION "Build the C version CLOUDSC using Serialbox" DEFAULT ON
    CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_C )

    set( CMAKE_C_STANDARD 11 )
    set( CMAKE_C_STANDARD_REQUIRED ON )

    ecbuild_add_library(
        TARGET dwarf-cloudsc-c-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/yoecldp_c.c
            cloudsc/yoethf_c.h
            cloudsc/yoethf_c.c
            cloudsc/yomcst_c.h
            cloudsc/yomcst_c.c
            cloudsc/load_state.h
            cloudsc/load_state.c
            cloudsc/cloudsc_c.h
            cloudsc/cloudsc_c.c
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.c
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.c
            cloudsc/mycpu.h
            cloudsc/mycpu.c
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c
        SOURCES dwarf_cloudsc.c
        LIBS dwarf-cloudsc-c-lib
    )

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-serial
        COMMAND bin/dwarf-cloudsc-c
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-omp
        COMMAND bin/dwarf-cloudsc-c
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        # The OpenMP C runtime for some reason doesn't appreciate setting the number of threads
        # in the application when OMP_NUM_THREADS is set, therefore we disable this here
        # OMP 4
        ENVIRONMENT OMP_STACKSIZE=1G
        CONDITION HAVE_OMP
    )

else()
    ecbuild_info( "Serialbox not found, disabling C prototype" )
endif()
