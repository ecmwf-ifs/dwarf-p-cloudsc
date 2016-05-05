list(APPEND dwarf3_cloudsc_srcs
cloudsc/cloud_layer.F90
cloudsc/cloudsc.F90
cloudsc/cloudsc_driver.F90
cloudsc/cloudsc_in.F90
cloudsc/cloudsc_out.F90
cloudsc/cuadjtq.F90
)

list(APPEND dwarf3_support_srcs
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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


ecbuild_add_library(TARGET dwarf3-cloudsc-libs
                    DEFINITIONS
                        ${CLOUDSC_DEFINITIONS}
                    SOURCES
                        ${dwarf3_cloudsc_srcs}
                        ${dwarf3_support_srcs}
)


ecbuild_add_executable(TARGET  dwarf3-cloudsc
                       SOURCES cloudsc_dwarf.F90 
                       LIBS    dwarf3-cloudsc-libs                                
)

