#!/usr/bin/python3

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from pathlib import Path
import serialbox as ser
import h5py
import numpy as np

class Config:

    def __init__(self, ser_dir, ser_prefix, h5_dir, h5_file):
        self.ser_dir = ser_dir
        self.ser_prefix = ser_prefix
        self.h5_dir = Path(h5_dir)
        self.h5_file = h5_file


def serialbox2hdf5(config):

    with h5py.File(config.h5_dir/config.h5_file, 'w') as output_file:

        input_data = ser.Serializer(ser.OpenModeKind.Read, config.ser_dir, config.ser_prefix)
        savepoint = input_data.savepoint_list()[0]

        # for field in input_data.fieldnames():
        #     print('{:12} {:12}'.format(field, ','.join(str(d) for d in input_data.get_field_metainfo(field).dims)))

        for field in input_data.fieldnames():
            data = input_data.read(field, savepoint)
            if len(data.shape) > 1:
                data = data.transpose()
            output_file.create_dataset(field, data=data)

        for field, value in input_data.global_metainfo:
            if isinstance(value, np.ndarray):
                data = value
                if len(data.shape) > 1:
                    data = data.transpose()
            else:
                data = np.array([value])
            output_file.create_dataset(field, shape=(1,), data=data)

        default_metainfo = {
                'KLON': np.array([100]),
                'KLEV': np.array([137])
            }
        for field, value in default_metainfo.items():
            if field not in input_data.global_metainfo.to_dict():
                output_file.create_dataset(field, shape=(1,), data=value)


def verify(config):
    with h5py.File(config.h5_dir/config.h5_file, 'r') as h5f:
        sf = ser.Serializer(ser.OpenModeKind.Read, config.ser_dir, config.ser_prefix)
        savepoint = sf.savepoint_list()[0]

        h5f_fields = set(list(h5f))
        sf_fields = set(sf.fieldnames()) | set(f for f, _ in sf.global_metainfo)

        assert(h5f_fields == sf_fields)

        for field in sf.fieldnames():
            data = sf.read(field, savepoint)
            if len(data.shape) > 1:
                assert(np.all(data.transpose() == h5f[field]))
            else:
                assert(np.all(data == h5f[field]))

        for field, value in sf.global_metainfo:
            assert(h5f[field][0] == value)


if __name__ == '__main__':
    configs = [
        Config(ser_dir='../data', ser_prefix='input', h5_dir='../config-files', h5_file='input.h5'),
        Config(ser_dir='../data', ser_prefix='reference', h5_dir='../config-files', h5_file='reference.h5')
    ]
    for config in configs:
        serialbox2hdf5(config)
        verify(config)
