# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( ${CMAKE_VERSION} VERSION_LESS "3.25" )
  if ( DWARF_P_CLOUDSC_ENABLE_ACC OR (NOT DEFINED DWARF_P_CLOUDSC_ENABLE_ACC AND (ENABLE_ACC OR NOT DEFINED ENABLE_ACC)) )
    # See https://gitlab.kitware.com/cmake/cmake/-/issues/23691, fixed in CMake 3.25
    # (TL;DR: FindOpenACC sets OpenACC_<LANG>_FOUND correctly but does not set
    #  OpenACC_FOUND unless all three C, CXX, and Fortran have been found - even if
    #  only one language has been requested via COMPONENTS)
    find_package( OpenACC COMPONENTS Fortran )
    if( OpenACC_Fortran_FOUND )
      set( OpenACC_FOUND ON )
    endif()
  endif()
endif()
