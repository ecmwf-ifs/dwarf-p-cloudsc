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

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Clang")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    ecbuild_add_fortran_flags("-O2"                              )
    ecbuild_add_fortran_flags("-g"                               )
    ecbuild_add_fortran_flags("-fPIC"                            )
    ecbuild_add_fortran_flags("-fopenmp"                         )
    ecbuild_add_fortran_flags("-ffpe-trap=invalid,zero,overflow" )
    ecbuild_add_fortran_flags("-fstack-arrays"                   )
    ecbuild_add_fortran_flags("-fconvert=big-endian"             )
    ecbuild_add_fortran_flags("-fbacktrace"                      )
    ecbuild_add_fortran_flags("-fno-second-underscore"           )
    ecbuild_add_fortran_flags("-ffree-form"                      )
    ecbuild_add_fortran_flags("-ffast-math"                      )
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    ecbuild_add_fortran_flags("-O2"                           )
    ecbuild_add_fortran_flags("-g"                            )
    ecbuild_add_fortran_flags("-fpic"                         )
    ecbuild_add_fortran_flags("-qopenmp"                      )
    ecbuild_add_fortran_flags("-qopenmp-threadprivate compat" )
    ecbuild_add_fortran_flags("-assume byterecl"              )
    ecbuild_add_fortran_flags("-fstack-arrays"                )
    ecbuild_add_fortran_flags("-convert big-endian"           )
    ecbuild_add_fortran_flags("-traceback"                    )
    ecbuild_add_fortran_flags("-align array64byte"            )
    ecbuild_add_fortran_flags("-warn nounused,nouncalled"     )
    ecbuild_add_fortran_flags("-xHost"                        )
    ecbuild_add_fortran_flags("-finline-functions"            )
    ecbuild_add_fortran_flags("-finline-limit=500"            )
    ecbuild_add_fortran_flags("-Winline"                      )
    ecbuild_add_fortran_flags("-Wl, --as-needed"              )
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    ecbuild_add_fortran_flags("-O2"                    )
    ecbuild_add_fortran_flags("-g"                     )
    ecbuild_add_fortran_flags("-fpic"                  )
    ecbuild_add_fortran_flags("-mp"                    )
    ecbuild_add_fortran_flags("-mp=bind,allcores,numa" )
    ecbuild_add_fortran_flags("-Mbyteswapio"           )     
    ecbuild_add_fortran_flags("-Mstack_arrays"         )
    ecbuild_add_fortran_flags("-Mrecursive"            )
    ecbuild_add_fortran_flags("-Ktrap=fp"              )
    ecbuild_add_fortran_flags("-Kieee"                 )
    ecbuild_add_fortran_flags("-Mdaz"                  )
    ecbuild_add_fortran_flags("-Wl, --as-needed"       )
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
    ecbuild_add_fortran_flags("-O0"                    )
    ecbuild_add_fortran_flags("-hflex_mp=conservative" )
    ecbuild_add_fortran_flags("--hfp1"                 )
    ecbuild_add_fortran_flags("-hadd_paren"            )
    ecbuild_add_fortran_flags("-hcontiguous"           )
    ecbuild_add_fortran_flags("-hbyteswapio"           )
    ecbuild_add_fortran_flags("-Ktrap=fp"              )
    ecbuild_add_fortran_flags("-Wl, --as-needed"       )
endif()

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

