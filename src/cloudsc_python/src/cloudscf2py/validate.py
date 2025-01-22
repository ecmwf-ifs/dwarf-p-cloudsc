# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import numpy as np


def validate(fields, ref_fields, kidia, kfdia, transpose=False):
    """
    Compare norms for given output fields against reference fields.

    By default, this assumes all data in row-major C layout. The
    :data:`transpose` flag can be used to change result data arrays
    into column-major Fortran layout.

    Parameters
    ----------
    fields : dict
        Dict of result fields from the run
    ref_fields : dict
        Dict of reference fields from file against which to compare
    kidia : int
        Start index for horizontal column dimension
    kfdia : int
        End index for horizontal column dimension
    transpose : bool
        Flag to apply a transpose step to Fortran-layout inputs
    """

    _field_names = [
        'plude', 'pcovptot', 'prainfrac_toprfz', 'pfsqlf', 'pfsqif',
        'pfcqlng', 'pfcqnng', 'pfsqrf', 'pfsqsf', 'pfcqrng', 'pfcqsng',
        'pfsqltur', 'pfsqitur', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t',
        'tendency_loc_cld'
    ]
    ngptot = kfdia - kidia + 1

    print(
        "             Variable Dim             MinValue             MaxValue"
        "            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%"
    )
    for name in _field_names:
        if transpose:
            # Transpose input Fortran arrays back to C-layout
            fields[name] = fields[name].transpose()

        if len(fields[name].shape) == 1:
            f = fields[name][kidia-1:kfdia]
            ref = ref_fields[name][kidia-1:kfdia]
        elif len(fields[name].shape) == 2:
            f = fields[name][:,kidia-1:kfdia]
            ref = ref_fields[name][:,kidia-1:kfdia]
        elif len(fields[name].shape) == 3:
            f = fields[name][:,:,kidia-1:kfdia]
            ref = ref_fields[name][:,:,kidia-1:kfdia]
        else:
            f = fields[name]
            ref = ref_fields[name]
        zsum = np.sum(np.absolute(ref))
        zerrsum = np.sum(np.absolute(f - ref))
        zeps = np.finfo(np.float64).eps
        print(
            ' {fname:>20}     {fmin:20.13e}  {fmax:20.13e}  {absmax:20.13e} '
            ' {absavg:20.13e}  {maxrel:20.13e}'.format(
                fname=name.upper(), fmin=f.min(), fmax=f.max(),
                absmax=np.absolute(f - ref).max(),
                absavg=np.sum(np.absolute(f - ref)) / ngptot,
                maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum) if zsum < zeps else zerrsum/zsum)
            )
        )
