include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

ecbuild_add_library(
    TARGET  dwarf-cloudsc-v2-lib
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
    INCLUDES
        ${SERIALBOX_INCLUDE_DIRS}
    SOURCES 
        modules/parkind1.F90
        modules/yoecldp.F90
        modules/yommp0.F90
        modules/yomcst.F90
        modules/yoethf.F90
        modules/yoephli.F90
        modules/yomjfh.F90
        modules/routines.F90
        modules/abor1.F90
        cloudsc/load_state_mod.F90
        cloudsc/cloudsc_global_state_mod.F90
        cloudsc/cloudsc_driver_mod.F90
        cloudsc/cloudsc.F90
        cloudsc/cuadjtq.F90
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
