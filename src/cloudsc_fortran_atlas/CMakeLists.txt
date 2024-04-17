# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_ATLAS
    DESCRIPTION "Build the Fortran version CLOUDSC using Atlas and Serialbox" DEFAULT ON
    CONDITION atlas_FOUND AND (Serialbox_FOUND OR HDF5_FOUND)
)

if( HAVE_CLOUDSC_FORTRAN_ATLAS )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-fortran-atlas
        SOURCES
            cloudsc_global_atlas_state_mod.F90
            expand_atlas_mod.F90
            validate_atlas_mod.F90
            cloudsc_driver_mod.F90
            cloudsc.F90
            dwarf_cloudsc_atlas.F90
        LIBS
            cloudsc-common-lib
            atlas_f
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

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

    ecbuild_add_test(
        TARGET dwarf-cloudsc-fortran-atlas-serial
        COMMAND bin/dwarf-cloudsc-fortran-atlas
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-fortran-atlas-omp
        COMMAND bin/dwarf-cloudsc-fortran-atlas
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-fortran-atlas-mpi
        COMMAND bin/dwarf-cloudsc-fortran-atlas
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-fortran-mpi-atlas-omp
        COMMAND bin/dwarf-cloudsc-fortran-atlas
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )

endif()
