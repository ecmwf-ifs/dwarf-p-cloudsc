# (C) Copyright 2018- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

##############################################################################
#.rst:
#
# cloudsc_find_python_venv
# ========================
#
# Find Python 3 inside a virtual environment. ::
#
#   cloudsc_find_python_venv(VENV_PATH)
#
# It finds the Python3 Interpreter from a virtual environment at
# the given location (`VENV_PATH`)
#
# Options
# -------
#
# :VENV_PATH: The path to the virtual environment
#
# Output variables
# ----------------
# :Python3_FOUND:       Exported into parent scope from FindPython3
# :Python3_EXECUTABLE:  Exported into parent scope from FindPython3
# :Python3_VENV_BIN:    The path to the virtual environment's `bin` directory
# :ENV{VIRTUAL_ENV}:    Environment variable with the virtual environment directory,
#                       emulating the activate script
#
##############################################################################

function( cloudsc_find_python_venv VENV_PATH )

    # Update the environment with VIRTUAL_ENV variable (mimic the activate script)
    set( ENV{VIRTUAL_ENV} ${VENV_PATH} )

    # Change the context of the search to only find the venv
    set( Python3_FIND_VIRTUALENV ONLY )

    # Unset Python3_EXECUTABLE because it is also an input variable
    #  (see documentation, Artifacts Specification section)
    unset( Python3_EXECUTABLE )

    # Launch a new search
    find_package( Python3 COMPONENTS Interpreter REQUIRED QUIET )

    # Find the binary directory of the virtual environment
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; import os.path; print(os.path.dirname(sys.executable), end='')"
        OUTPUT_VARIABLE Python3_VENV_BIN
    )

    # Forward variables to parent scope
    foreach ( _VAR_NAME Python3_FOUND Python3_EXECUTABLE Python3_VENV_BIN )
        set( ${_VAR_NAME} ${${_VAR_NAME}} PARENT_SCOPE )
    endforeach()

endfunction()

##############################################################################
#.rst:
#
# cloudsc_create_python_venv
# ==========================
#
# Find Python 3 and create a virtual environment. ::
#
#   cloudsc_create_python_venv(VENV_PATH)
#
# Installation procedure
# ----------------------
#
# It creates a virtual environment at the given location (`VENV_PATH`)
#
# Options
# -------
#
# :VENV_PATH: The path to use for the virtual environment
#
##############################################################################

function( cloudsc_create_python_venv VENV_PATH )

    # Discover only system install Python 3
    set( Python3_FIND_VIRTUALENV STANDARD )
    find_package( Python3 COMPONENTS Interpreter REQUIRED )

    # Create a loki virtualenv
    message( STATUS "Create Python virtual environment ${VENV_PATH}" )
    execute_process( COMMAND ${Python3_EXECUTABLE} -m venv "${VENV_PATH}" )

endfunction()

##############################################################################
#.rst:
#
# cloudsc_setup_python_venv
# =================
#
# Find Python 3, create a virtual environment and make it available. ::
#
#   cloudsc_setup_python_venv(VENV_PATH)
#
# It combines calls to `cloudsc_create_python_venv` and `cloudsc_find_python_venv`
#
# Options
# -------
#
# :VENV_PATH: The path to use for the virtual environment
#
# Output variables
# ----------------
# :Python3_FOUND:       Exported into parent scope from FindPython3
# :Python3_EXECUTABLE:  Exported into parent scope from FindPython3
# :Python3_VENV_BIN:    The path to the virtual environment's `bin` directory
# :ENV{VIRTUAL_ENV}:    Environment variable with the virtual environment directory,
#                       emulating the activate script
#
##############################################################################

macro( cloudsc_setup_python_venv VENV_PATH )

    # Create the virtual environment
    cloudsc_create_python_venv( ${VENV_PATH} )

    # Discover Python in the virtual environment and set-up variables
    cloudsc_find_python_venv( ${VENV_PATH} )

endmacro()
