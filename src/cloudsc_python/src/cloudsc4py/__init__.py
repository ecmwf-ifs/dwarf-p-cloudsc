# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import os

import gt4py.cartesian.config as gt_config

import cloudsc4py.physics


# customize compilation/linking of GT4Py generated code
cxxflags = os.environ.get("CXXFLAGS", "")
if cxxflags != "":
    gt_config.build_settings["extra_compile_args"]["cxx"] += cxxflags.split(" ")

lflags = os.environ.get("LFLAGS", "")
if lflags != "":
    gt_config.build_settings["extra_link_args"] += lflags.split(" ")
