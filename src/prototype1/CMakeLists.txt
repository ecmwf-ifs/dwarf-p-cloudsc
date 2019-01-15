ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

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

list(APPEND CLOUDSC_DEFINITIONS ECMWF)
list(APPEND CLOUDSC_DEFINITIONS )
if(${EC_OS_NAME} MATCHES "linux" )
  list( APPEND CLOUDSC_DEFINITIONS LINUX)
endif()

if(${EC_OS_NAME} MATCHES "macosx")
  list( APPEND CLOUDSC_DEFINITIONS DARWIN)
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Note, that we link statically here, since Serialbox only exposes static Fortran libs.
# Serailbox is used to generate the off-line input and reference data for alternative
# versions of this dwarf, eg. prototype2.
#
# We should return to dynamic linking once Serialbox::SerialboxFortranShared is available.

ecbuild_add_library(TARGET dwarf-P-cloudMicrophysics-IFSScheme-libs
                    TYPE STATIC
                    LINKER_LANGUAGE CXX
                    DEFINITIONS
                        ${CLOUDSC_DEFINITIONS}
                    SOURCES
                        ${dwarf_P_cloudMicrophysics_IFSScheme_srcs}
                        ${dwarf_P_cloudMicrophysics_IFSScheme_support_srcs}
)

if( Serialbox_FOUND )
  target_compile_definitions( dwarf-P-cloudMicrophysics-IFSScheme-libs PRIVATE HAVE_SERIALBOX )
  target_link_libraries( dwarf-P-cloudMicrophysics-IFSScheme-libs PRIVATE
                         Serialbox::SerialboxFortranStatic )
endif()

target_link_libraries( dwarf-P-cloudMicrophysics-IFSScheme-libs PRIVATE  OpenMP::OpenMP_Fortran )

ecbuild_add_executable(TARGET  dwarf-P-cloudMicrophysics-IFSScheme
                       SOURCES cloudsc_dwarf.F90 
                       LIBS
		           dwarf-P-cloudMicrophysics-IFSScheme-libs
)
