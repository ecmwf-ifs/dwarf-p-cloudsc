### use of cloudsc-run for tests
ecbuild_add_option( FEATURE CLOUDSC_RUN
                    DEFAULT ON
                    DESCRIPTION "Use cloudsc/tools/cloudsc-run to run cloudsc tests" )

if( HAVE_CLOUDSC_RUN )
    # set( CMAKE_CROSSCOMPILING_EMULATOR ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
    # set( MPIEXEC_EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
    # set( MPIEXEC_NUMPROC_FLAG='-n' )
    set( CMAKE_CROSSCOMPILING_EMULATOR ${PROJECT_SOURCE_DIR}/tools/cloudsc-run )
endif()

