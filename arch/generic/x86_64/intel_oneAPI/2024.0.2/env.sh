# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Source me to get the correct configure/build/run environment

# As this is a generic arch, it won't load all necessary dependencies, e.g.
# python, intel oneAPI compilers, hdf5, cmake, etc

export CC=icx
export CXX=icpx
export F90=ifx
export FC=ifx
export F77=ifx

ulimit -s unlimited

set -x
echo foo
# Restore tracing to stored setting
#{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
set +x
echo bar
export ECBUILD_TOOLCHAIN="./toolchain.cmake"
