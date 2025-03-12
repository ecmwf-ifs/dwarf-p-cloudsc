# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

ecbuild_add_executable(
    TARGET dwarf-cloudsc-fortran-atlas-${prec}
    SOURCES
        cloudsc_global_atlas_state_mod.F90
        expand_atlas_mod.F90
        validate_atlas_mod.F90
        cloudsc_driver_mod.F90
        cloudsc.F90
        dwarf_cloudsc_atlas.F90
    LIBS
        cloudsc-common-${prec}-atlas
        atlas_f
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)

set_target_properties( dwarf-cloudsc-fortran-atlas-${prec}
    PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module"
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-fortran-atlas-${prec}-serial
    COMMAND dwarf-cloudsc-fortran-atlas-${prec}
    ARGS 1 100 16
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
)
ecbuild_add_test(
    TARGET dwarf-cloudsc-fortran-atlas-${prec}-omp
    COMMAND dwarf-cloudsc-fortran-atlas-${prec}
    ARGS 4 100 16
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 4
    CONDITION HAVE_OMP
)
ecbuild_add_test(
    TARGET dwarf-cloudsc-fortran-atlas-${prec}-mpi
    COMMAND dwarf-cloudsc-fortran-atlas-${prec}
    ARGS 1 100 16
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    MPI 2
    OMP 1
    CONDITION HAVE_MPI
)
ecbuild_add_test(
    TARGET dwarf-cloudsc-fortran-mpi-atlas-${prec}-omp
    COMMAND dwarf-cloudsc-fortran-atlas-${prec}
    ARGS 4 100 16
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    MPI 2
    OMP 4
    CONDITION HAVE_OMP AND HAVE_MPI
)
