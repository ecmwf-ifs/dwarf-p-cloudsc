ecbuild_add_library(
    TARGET  dwarf-cloudsc-v2-lib
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
    INCLUDES
        ${SERIALBOX_INCLUDE_DIRS}
    SOURCES 
        modules/parkind1.F90    	    
        cloudsc/cloudsc_driver.F90
    LIBS
        ${SERIALBOX_FORTRAN_LIBRARIES}
)

ecbuild_add_executable(
    TARGET  dwarf-cloudsc-v2
    SOURCES
        dwarf_cloudsc.F90 
    LIBS
        dwarf-cloudsc-v2-lib
)
