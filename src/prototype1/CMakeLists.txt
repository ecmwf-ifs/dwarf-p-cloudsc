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


# Default install dir if not given by the user
if(NOT SERIALBOX_ROOT)
  set(SERIALBOX_ROOT "/home/rd/naml/apps/serialbox2/install")
  message(STATUS "Set serialbox2 install directory to: ${SERIALBOX_ROOT}")
endif()
set( SERIALBOX_USE_SHARED_LIBS ON )
# Find installed version of serialbox
find_package(Serialbox REQUIRED COMPONENTS Fortran HINTS
  "${SERIALBOX_ROOT}/cmake/"
)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


ecbuild_add_library(TARGET dwarf-P-cloudMicrophysics-IFSScheme-libs
                    DEFINITIONS
                        ${CLOUDSC_DEFINITIONS}
		    INCLUDES
		        ${SERIALBOX_INCLUDE_DIRS}
                    SOURCES
                        ${dwarf_P_cloudMicrophysics_IFSScheme_srcs}
                        ${dwarf_P_cloudMicrophysics_IFSScheme_support_srcs}
                    LIBS
		        ${SERIALBOX_FORTRAN_LIBRARIES}
)

ecbuild_add_executable(TARGET  dwarf-P-cloudMicrophysics-IFSScheme
                       SOURCES cloudsc_dwarf.F90 
                       LIBS
		           dwarf-P-cloudMicrophysics-IFSScheme-libs
)
