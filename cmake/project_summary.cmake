# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( HAVE_ACC )
    ecbuild_info("")
    ecbuild_info("OpenACC: ")
    if( OpenACC_C_FOUND )
        ecbuild_info("  - C libraries:           ${OpenACC_C_LIBRARIES}" )
        ecbuild_info("  - C flags:               ${OpenACC_C_FLAGS}" )
    else()
        ecbuild_info("  - OpenACC_C_FOUND:       FALSE")
    endif()
    if( OpenACC_Fortran_FOUND )
        ecbuild_info("  - Fortran libraries:     ${OpenACC_Fortran_LIBRARIES}" )
        ecbuild_info("  - Fortran flags:         ${OpenACC_Fortran_FLAGS}" )
    else()
        ecbuild_info("  - OpenACC_Fortran_FOUND: FALSE")
    endif()
endif()

if( HAVE_OMP )
    ecbuild_info("")
    ecbuild_info("OpenMP: ")
    if( OpenMP_C_FOUND )
        ecbuild_info("  - C libraries:           ${OpenMP_C_LIBRARIES}" )
        ecbuild_info("  - C flags:               ${OpenMP_C_FLAGS}" )
    else()
        ecbuild_info("  - OpenMP_C_FOUND:        FALSE")
    endif()
    if( OpenMP_Fortran_FOUND )
        ecbuild_info("  - Fortran libraries:     ${OpenMP_Fortran_LIBRARIES}" )
        ecbuild_info("  - Fortran flags:         ${OpenMP_Fortran_FLAGS}" )
    else()
        ecbuild_info("  - OpenMP_Fortran_FOUND: FALSE")
    endif()
endif()

ecbuild_info("")
