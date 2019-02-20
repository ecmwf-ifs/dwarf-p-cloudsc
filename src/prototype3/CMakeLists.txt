ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)


if( HAVE_SERIALBOX )

    ecbuild_add_library(
        TARGET dwarf-cloudsc-v3-lib
	INSTALL_HEADERS LISTED
        SOURCES
	    cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.c
	    cloudsc/cloudsc_c.h
            cloudsc/cloudsc_c.c
    )
    target_link_libraries( dwarf-cloudsc-v3-lib PRIVATE
                           Serialbox::SerialboxFortranShared
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
