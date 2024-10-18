# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI
    DESCRIPTION "Use Loki source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

# Define the CLAW-based Loki dwarf variants as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI_CLAW
    DESCRIPTION "Build optimized CPU and GPU version of CLOUDSC derived from Loki+CLAW" DEFAULT OFF
    CONDITION HAVE_CLOUDSC_LOKI
)

if( HAVE_CLOUDSC_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    # OFP frontend cannot deal with statement functions, so we toggle them here
    set( CLOUDSC_DEFINE_STMT_FUNC "" )
    if(NOT "${LOKI_FRONTEND}" STREQUAL "ofp")
    set( CLOUDSC_DEFINE_STMT_FUNC CLOUDSC_STMT_FUNC )
    endif()

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem/cloudsc.idem.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE idem
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
            loki-idem/cloudsc.idem.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-serial
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )


    ####################################################
    ##  "Single Column Coalesced" (SCC) mode          ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc/cloudsc.scc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DATA_OFFLOAD
        REMOVE_OPENMP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc
        SOURCES
            dwarf_cloudsc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
            loki-scc/cloudsc.scc.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-serial
        COMMAND bin/dwarf-cloudsc-loki-scc
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )


    ####################################################
    ##  "Single Column Coalesced Field" (SCC-FIELD)   ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ##   * Uses Field API for device data offloading  ##
    ####################################################

    if ( HAVE_FIELD_API AND field_api_HAVE_ACC )
        set( CLOUDSC_FORTRAN_SRC_PATH ${PROJECT_SOURCE_DIR}/src/cloudsc_fortran )
    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-field/cloudsc.scc_field.F90
            loki-scc-field/cloudsc_driver_field_mod.scc_field.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-field
        DEPENDS
            cloudsc.F90
            ${CLOUDSC_FORTRAN_SRC_PATH}/cloudsc_driver_field_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-field
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki_field_offload.config
        CPP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
            ${COMMON_MODULE}/cloudsc_field_state_mod.F90
            ${COMMON_MODULE}/cloudsc_flux_type_mod.F90
            ${COMMON_MODULE}/cloudsc_aux_type_mod.F90
            ${COMMON_MODULE}/cloudsc_state_type_mod.F90
        SOURCES
            ${CLOUDSC_FORTRAN_SRC_PATH}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-field
        SOURCES
            ${PROJECT_SOURCE_DIR}/src/cloudsc_gpu/dwarf_cloudsc_gpu.F90
            loki-scc-field/cloudsc.scc_field.F90
            loki-scc-field/cloudsc_driver_field_mod.scc_field.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-field
    	PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-field
    )


    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-field
        COMMAND bin/dwarf-cloudsc-loki-scc-field
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )

endif()

    # Create symlink for the input data
    if( HAVE_SERIALBOX )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
    endif()

    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()
endif()
