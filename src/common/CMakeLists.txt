ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

ecbuild_add_library( TARGET cloudsc-common-lib
    SOURCES
        module/parkind1.F90
        module/yoecldp.F90
        module/yomcst.F90
        module/yoethf.F90
        module/yoephli.F90
        module/yomphyder.F90
        module/fcttre_mod.F90
        module/fccld_mod.F90
        module/routines.F90
        module/abor1.F90
        module/timer_mod.F90
        module/mycpu.c
        module/ec_pmon_mod.F90
        module/load_array_mod.F90
        module/validate_mod.F90
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)
target_include_directories( cloudsc-common-lib PUBLIC $<INSTALL_INTERFACE:include> )
target_include_directories( cloudsc-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories( cloudsc-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module>)

if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
    target_link_libraries( cloudsc-common-lib PRIVATE OpenMP::OpenMP_Fortran )
endif()

if( TARGET Serialbox::Serialbox_Fortran )
    target_link_libraries( cloudsc-common-lib PRIVATE Serialbox::Serialbox_Fortran )
endif()
