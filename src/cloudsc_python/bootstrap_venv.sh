# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#!/bin/bash

PYTHON=$(which python3)
PIP_UPGRADE=${PIP_UPGRADE:-1}
VENV=${VENV:-venv}
FRESH_INSTALL=${FRESH_INSTALL:-1}
INSTALL_PRE_COMMIT=${INSTALL_PRE_COMMIT:-1}
INSTALL_CUPY=${INSTALL_CUPY:-0}
CUPY_VERSION=${CUPY_VERSION:-cupy}


function install()
{
  # activate environment
  source "$VENV"/bin/activate

  # upgrade pip and setuptools
  if [ "$PIP_UPGRADE" -ne 0 ]; then
    pip install --upgrade pip setuptools wheel
  fi

  # install cloudsc4py
  pip install -e .

  # install gt sources
  python -m gt4py.gt_src_manager install

  # setup gt4py cache
  mkdir -p gt_cache
  echo -e "\nexport GT_CACHE_ROOT=$PWD/gt_cache" >> "$VENV"/bin/activate

  # install cupy
  if [ "$INSTALL_CUPY" -eq 1 ]; then
    pip install "$CUPY_VERSION"
  fi

  # install development packages
  pip install -r requirements_dev.txt

  # install pre-commit
  if [ "$INSTALL_PRE_COMMIT" -eq 1 ]; then
    pre-commit install
  fi

  # deactivate environment
  deactivate
}


if [ "$FRESH_INSTALL" -eq 1 ]; then
  echo -e "Creating new environment..."
  rm -rf "$VENV"
  $PYTHON -m venv "$VENV"
fi


install || deactivate


echo -e ""
echo -e "Command to activate environment:"
echo -e "\t\$ source $VENV/bin/activate"
echo -e ""
echo -e "Command to deactivate environment:"
echo -e "\t\$ deactivate"
echo -e ""
