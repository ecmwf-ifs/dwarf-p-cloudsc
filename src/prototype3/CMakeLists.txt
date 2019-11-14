ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)


if( HAVE_SERIALBOX )

    ecbuild_add_library(
        TARGET dwarf-cloudsc-v3-lib
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
    target_link_libraries( dwarf-cloudsc-v3-lib PRIVATE
                           Serialbox::SerialboxCShared
                           OpenMP::OpenMP_Fortran )
    target_include_directories( dwarf-cloudsc-v3-lib PUBLIC $<INSTALL_INTERFACE:include> )
    target_include_directories( dwarf-cloudsc-v3-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc> )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-v3
        SOURCES
            dwarf_cloudsc.c
    )
    target_link_libraries( dwarf-cloudsc-v3 dwarf-cloudsc-v3-lib )


else()
    message(STATUS "Serialbox not found, disabling prototype3" )
endif()
