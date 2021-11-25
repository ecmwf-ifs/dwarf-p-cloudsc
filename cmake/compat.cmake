# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( CMAKE_VERSION VERSION_LESS 3.16 )
    set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/compat/3.16;${CMAKE_MODULE_PATH} )
endif()

