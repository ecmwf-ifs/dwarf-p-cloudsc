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
	    cloudsc/yoethf_c.h
	    cloudsc/yomcst_c.h
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
    )
    target_link_libraries( dwarf-cloudsc-c-lib PRIVATE Serialbox::Serialbox_C )
    if( TARGET OpenMP::OpenMP_C )
        target_link_libraries( dwarf-cloudsc-c-lib PRIVATE OpenMP::OpenMP_C )
    endif()
    target_include_directories( dwarf-cloudsc-c-lib PUBLIC $<INSTALL_INTERFACE:include> )
    target_include_directories( dwarf-cloudsc-c-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc> )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c
        SOURCES dwarf_cloudsc.c
    )
    target_link_libraries(
        dwarf-cloudsc-c PRIVATE
        dwarf-cloudsc-c-lib
    )

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

else()
    message(STATUS "Serialbox not found, disabling C prototype" )
endif()
