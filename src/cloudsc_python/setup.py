# -*- coding: utf-8 -*-
from setuptools import setup
import sys


if sys.version_info.major < 3:
    print("Python 3.x is required.")
    sys.exit(1)


setup(use_scm_version=False)
