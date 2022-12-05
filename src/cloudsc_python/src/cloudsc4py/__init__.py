# -*- coding: utf-8 -*-
import os

import gt4py.config as gt_config

import cloudsc4py.physics


# customize compilation/linking of GT4Py generated code
cxxflags = os.environ.get("CXXFLAGS", "")
if cxxflags != "":
    gt_config.build_settings["extra_compile_args"]["cxx"] += cxxflags.split(" ")

lflags = os.environ.get("LFLAGS", "")
if lflags != "":
    gt_config.build_settings["extra_link_args"] += lflags.split(" ")
