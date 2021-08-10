ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_PROTOTYPE1
    DESCRIPTION "Build the original operational Fortran prototype" DEFAULT ON
)

list(APPEND dwarf_P_cloudMicrophysics_IFSScheme_srcs
cloudsc/cloud_layer.F90
cloudsc/cloudsc.F90
cloudsc/cloudsc_driver.F90
cloudsc/cloudsc_in.F90
cloudsc/cloudsc_out.F90
cloudsc/cuadjtq.F90
)

list(APPEND dwarf_P_cloudMicrophysics_IFSScheme_support_srcs
support/abor1.F90
support/diag_mod.F90
support/diff_mod.F90
support/expand_mod.F90
support/serialize_mod.F90
support/mycpu.c
support/parkind1.F90
support/routines.F90
support/surface_fields_mix.F90
support/timer_mod.F90
support/yoecldp.F90
support/yoephli.F90
support/yoethf.F90
support/yomcst.F90
support/yomct0.F90
support/yomct3.F90
support/yomdim.F90
support/yomdyn.F90
support/yomjfh.F90
support/yomlun.F90
support/yomlun_ifsaux.F90
support/yommp0.F90
support/yomphy2.F90
support/yomphyder.F90
)

if( HAVE_CLOUDSC_PROTOTYPE1 )

    # Define the binary build of the original CLOUDSC prototype
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    ecbuild_add_library(TARGET dwarf-P-cloudMicrophysics-IFSScheme-libs
        SOURCES
            ${dwarf_P_cloudMicrophysics_IFSScheme_srcs}
            ${dwarf_P_cloudMicrophysics_IFSScheme_support_srcs}
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    if( TARGET Serialbox::Serialbox_Fortran )
        target_link_libraries( dwarf-P-cloudMicrophysics-IFSScheme-libs PRIVATE Serialbox::Serialbox_Fortran )
    endif()

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-P-cloudMicrophysics-IFSScheme-libs PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    ecbuild_add_executable(TARGET dwarf-P-cloudMicrophysics-IFSScheme
        SOURCES cloudsc_dwarf.F90
    )
    target_link_libraries( dwarf-P-cloudMicrophysics-IFSScheme dwarf-P-cloudMicrophysics-IFSScheme-libs )

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/cloudsc.bin ${CMAKE_CURRENT_BINARY_DIR}/../../../cloudsc.bin )
endif()
