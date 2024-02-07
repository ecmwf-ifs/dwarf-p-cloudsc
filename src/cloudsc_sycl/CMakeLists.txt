# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_SYCL
    DESCRIPTION "Build the SYCL version CLOUDSC using Serialbox" DEFAULT ON
    CONDITION Serialbox_FOUND AND HAVE_SYCL
)

if( HAVE_CLOUDSC_SYCL )

    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
	    cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c.kernel
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
	PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            Serialbox::Serialbox_C
    	    $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    add_sycl_to_target(
	    TARGET
	    	dwarf-cloudsc-scc-sycl-lib
	    SOURCES
	    	cloudsc/cloudsc_c.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-sycl-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-serial
        COMMAND bin/dwarf-cloudsc-scc-sycl
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    ######

    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-hoist-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_hoist.kernel
            cloudsc/cloudsc_driver_hoist.h
            cloudsc/cloudsc_driver_hoist.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
	    Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    add_sycl_to_target(
	    TARGET
	    	dwarf-cloudsc-scc-hoist-sycl-lib
	    SOURCES
	    	cloudsc/cloudsc_c_hoist.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-hoist-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-hoist-sycl-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-hoist-serial
        COMMAND bin/dwarf-cloudsc-scc-hoist-sycl
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    ######

    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-k-caching-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_k_caching.kernel
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver_k_caching.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
	    Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    add_sycl_to_target(
	    TARGET 
	    	dwarf-cloudsc-scc-k-caching-sycl-lib 
	    SOURCES 
	    	cloudsc/cloudsc_c_k_caching.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-k-caching-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-k-caching-sycl-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-k-caching-serial
        COMMAND bin/dwarf-cloudsc-scc-k-caching-sycl
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

else()
	ecbuild_info( "Serialbox not found, disabling SYCL version" )
endif()
