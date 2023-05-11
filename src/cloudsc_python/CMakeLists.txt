# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_PYTHON_F2PY
    DESCRIPTION "Build the pure Python variant from Loki transpilation" DEFAULT OFF
    CONDITION HDF5_FOUND
)

if( HAVE_CLOUDSC_PYTHON_F2PY )

    # Utilities to manage Python virtual environments
    include( python_venv )

    # Set up a custom venv for this variant and install the necessary dependencies
    set( cloudsc_VENV_PATH ${CMAKE_BINARY_DIR}/venv_cloudsc )
    setup_python_venv( ${cloudsc_VENV_PATH} )

    if( NOT Python3_EXECUTABLE )
        ecbuild_error("[PyIface] Could not find Python3 executable in virtualenv")
    endif()

    # Update to latest pip versionxs
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade pip)

    # Install the "cloudscf2py" Python package and runner in editable mode
    add_custom_command( OUTPUT ${Python3_VENV_BIN}/cloudsc_f2py.py
        COMMAND ${Python3_EXECUTABLE} -m pip install -e ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "[CLOUDSC-Python] Installing cloudscf2py into virtualenv [${cloudsc_VENV_PATH}]"
    )

    # Copy the CLI driver script into the bin directory for execution
    add_custom_command( OUTPUT ${CMAKE_BINARY_DIR}/bin/cloudsc_f2py.py
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${Python3_VENV_BIN}/cloudsc_f2py.py ${CMAKE_BINARY_DIR}/bin/cloudsc_f2py.py
        DEPENDS ${Python3_VENV_BIN}/cloudsc_f2py.py
        COMMENT "[CLOUDSC-Python] Linking Python driver scripts from virtualenv [${cloudsc_VENV_PATH}]"
    )

    # Add runner script as a custom executable target
    add_custom_target( cloudsc-f2py ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/cloudsc_f2py.py )

    ecbuild_add_test(
        COMMAND bin/cloudsc_f2py.py
        ARGS --ngptot=100 --nproma=16
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
    )


    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5
	  ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5
      )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5
	  ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5
      )
    endif()

endif()
