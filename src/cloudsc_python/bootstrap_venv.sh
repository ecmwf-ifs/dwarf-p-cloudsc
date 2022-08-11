#!/bin/bash

PYTHON=$(which python3)
PIP_UPGRADE=${PIP_UPGRADE:-1}
VENV=${VENV:-venv}
FRESH_INSTALL=${FRESH_INSTALL:-1}
EXTRAS=${EXTRAS:-}


function install()
{
  # activate environment
  source $VENV/bin/activate

  # upgrade pip and setuptools
  if [ "$PIP_UPGRADE" -ne 0 ]; then
    pip install --upgrade pip setuptools wheel
  fi

  # install cloudsc4py
  if [ -z "$EXTRAS" ]; then
    pip install -e .
  else
    pip install -e .["$EXTRAS"]
  fi

  # install gt sources
  python -m gt4py.gt_src_manager install

  # install development packages
  pip install -r requirements_dev.txt

  # deactivate environment
  deactivate
}


if [ "$FRESH_INSTALL" -eq 1 ]; then
  echo -e "Creating new environment..."
  rm -rf $VENV
  $PYTHON -m venv $VENV
fi


install || deactivate


echo -e ""
echo -e "Command to activate environment:"
echo -e "\t\$ source $VENV/bin/activate"
echo -e ""
echo -e "Command to deactivate environment:"
echo -e "\t\$ deactivate"
echo -e ""
