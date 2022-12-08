# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.



if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  ecbuild_add_fortran_flags("-Mlarge_arrays")

  # should really be part of configuration, or ecbuild default?
  ecbuild_add_fortran_flags("-traceback"      BUILD DEBUG )
  ecbuild_add_fortran_flags("-fast"           BUILD RELEASE )
  ecbuild_add_fortran_flags("-gopt -fast"     BUILD RELWITHDEBINFO )
endif()



macro( cloudsc_add_compile_options )
  set( options  )
  set( single_value_args FLAGS )
  set( multi_value_args SOURCES )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to cloudsc_add_compile_flags(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()
  if(NOT _PAR_SOURCES)
    ecbuild_critical("SOURCES keyword missing to cloudsc_add_compile_flags()")
  endif()
  if(NOT _PAR_FLAGS)
    ecbuild_critical("FLAGS keyword missing to cloudsc_add_compile_flags()")
  endif()
  foreach( _file ${_PAR_SOURCES} )
    ecbuild_warn("Adding custom compile flags for file ${_file} : [${_PAR_FLAGS}]")
    if( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file} )
        ecbuild_error("${_file} does not exist")
    endif()
    set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS "${_PAR_FLAGS}" )
  endforeach()
endmacro()

