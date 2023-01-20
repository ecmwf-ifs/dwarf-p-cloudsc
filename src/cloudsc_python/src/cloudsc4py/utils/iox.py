# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
from datetime import timedelta
from functools import lru_cache
import h5py
import numpy as np
from pydantic import BaseModel
from typing import TYPE_CHECKING

from cloudsc4py.utils.f2py import ported_method

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Optional, Type

    from cloudsc4py.framework.config import DataTypes


class YoecldpParameters(BaseModel):
    NCLDQI: int
    NCLDQL: int
    NCLDQR: int
    NCLDQS: int
    NCLDQV: int
    NCLV: int


class YoethfParameters(BaseModel):
    R2ES: float
    R3IES: float
    R3LES: float
    R4IES: float
    R4LES: float
    R5ALSCP: float
    R5ALVCP: float
    R5IES: float
    R5LES: float
    RALFDCP: float
    RALSDCP: float
    RALVDCP: float
    RKOOP1: float
    RKOOP2: float
    RTICE: float
    RTICECU: float
    RTWAT: float
    RTWAT_RTICECU_R: float
    RTWAT_RTICE_R: float


class YomcstParameters(BaseModel):
    RCPD: float
    RD: float
    RETV: float
    RG: float
    RLMLT: float
    RLSTT: float
    RLVTT: float
    RTT: float
    RV: float


class YrecldpParameters(BaseModel):
    LAERICEAUTO: bool
    LAERICESED: bool
    LAERLIQAUTOCP: bool
    LAERLIQAUTOCPB: bool
    LAERLIQAUTOLSP: bool
    LAERLIQCOLL: bool
    LCLDBUDGET: bool
    LCLDEXTRA: bool
    NAECLBC: int
    NAECLDU: int
    NAECLOM: int
    NAECLSS: int
    NAECLSU: int
    NAERCLD: int
    NBETA: int
    NCLDDIAG: int
    NCLDTOP: int
    NSHAPEP: int
    NSHAPEQ: int
    NSSOPT: int
    RAMID: float
    RAMIN: float
    RCCN: float
    RCCNOM: float
    RCCNSS: float
    RCCNSU: float
    RCLCRIT: float
    RCLCRIT_LAND: float
    RCLCRIT_SEA: float
    RCLDIFF: float
    RCLDIFF_CONVI: float
    RCLDMAX: float
    RCLDTOPCF: float
    RCLDTOPP: float
    RCL_AI: float
    RCL_APB1: float
    RCL_APB2: float
    RCL_APB3: float
    RCL_AR: float
    RCL_AS: float
    RCL_BI: float
    RCL_BR: float
    RCL_BS: float
    RCL_CDENOM1: float
    RCL_CDENOM2: float
    RCL_CDENOM3: float
    RCL_CI: float
    RCL_CONST1I: float
    RCL_CONST1R: float
    RCL_CONST1S: float
    RCL_CONST2I: float
    RCL_CONST2R: float
    RCL_CONST2S: float
    RCL_CONST3I: float
    RCL_CONST3R: float
    RCL_CONST3S: float
    RCL_CONST4I: float
    RCL_CONST4R: float
    RCL_CONST4S: float
    RCL_CONST5I: float
    RCL_CONST5R: float
    RCL_CONST5S: float
    RCL_CONST6I: float
    RCL_CONST6R: float
    RCL_CONST6S: float
    RCL_CONST7S: float
    RCL_CONST8S: float
    RCL_CR: float
    RCL_CS: float
    RCL_DI: float
    RCL_DR: float
    RCL_DS: float
    RCL_DYNVISC: float
    RCL_FAC1: float
    RCL_FAC2: float
    RCL_FZRAB: float
    RCL_FZRBB: float
    RCL_KA273: float
    RCL_KKAac: float
    RCL_KKAau: float
    RCL_KKBac: float
    RCL_KKBaun: float
    RCL_KKBauq: float
    RCL_KK_cloud_num_land: float
    RCL_KK_cloud_num_sea: float
    RCL_SCHMIDT: float
    RCL_X1I: float
    RCL_X1R: float
    RCL_X1S: float
    RCL_X2I: float
    RCL_X2R: float
    RCL_X2S: float
    RCL_X3I: float
    RCL_X3S: float
    RCL_X41: float
    RCL_X4R: float
    RCL_X4S: float
    RCOVPMIN: float
    RDENSREF: float
    RDENSWAT: float
    RDEPLIQREFDEPTH: float
    RDEPLIQREFRATE: float
    RICEHI1: float
    RICEHI2: float
    RICEINIT: float
    RKCONV: float
    RKOOPTAU: float
    RLCRITSNOW: float
    RLMIN: float
    RNICE: float
    RPECONS: float
    RPRC1: float
    RPRC2: float
    RPRECRHMAX: float
    RSNOWLIN1: float
    RSNOWLIN2: float
    RTAUMEL: float
    RTHOMO: float
    RVICE: float
    RVRAIN: float
    RVRFACTOR: float
    RVSNOW: float


class HDF5Reader:
    f: h5py.File
    data_types: DataTypes

    def __init__(self, filename: str, data_types: DataTypes) -> None:
        self.f = h5py.File(filename)
        self.data_types = data_types

    def __del__(self) -> None:
        self.f.close()

    def get_field(self, name: str) -> np.ndarray:
        ds = self.f.get(name, None)
        if ds is None:
            raise RuntimeError(f"Unknown field `{name}`.")

        if ds.ndim == 1:
            return self._get_field_1d(ds, name)
        elif ds.ndim == 2:
            return self._get_field_2d(ds, name)
        elif ds.ndim == 3:
            return self._get_field_3d(ds, name)
        else:
            raise RuntimeError(f"The field `{name}` has unexpected shape {ds.shape}.")

    @lru_cache
    def get_nlev(self) -> int:
        return self.f["KLEV"][0]

    @lru_cache
    def get_nlon(self) -> int:
        return self.f["KLON"][0]

    def get_timestep(self) -> timedelta:
        return timedelta(seconds=self._get_parameter_f("PTSPHY"))

    @ported_method(from_file="common/module/yoecldp.F90", from_line=86, to_line=91)
    def get_yoecldp_parameters(self) -> YoecldpParameters:
        return YoecldpParameters(
            **{"NCLV": 5, "NCLDQL": 1, "NCLDQI": 2, "NCLDQR": 3, "NCLDQS": 4, "NCLDQV": 5}
        )

    @ported_method(from_file="common/module/yoethf.F90", from_line=79, to_line=99)
    def get_yoethf_parameters(self) -> YoethfParameters:
        return self._initialize_parameters(YoethfParameters)

    @ported_method(from_file="common/module/yomcst.F90", from_line=167, to_line=177)
    def get_yomcst_parameters(self) -> YomcstParameters:
        return self._initialize_parameters(YomcstParameters)

    @ported_method(from_file="common/module/yoecldp.F90", from_line=242, to_line=370)
    def get_yrecldp_parameters(self) -> YrecldpParameters:
        return self._initialize_parameters(
            YrecldpParameters, get_parameter_name=lambda attr_name: "YRECLDP_" + attr_name
        )

    def _get_field_1d(self, ds: h5py.Dataset, name: str) -> np.ndarray:
        nlon = self.get_nlon()
        nlev = self.get_nlev()
        if nlon <= ds.shape[0] <= nlon + 1 or nlev <= ds.shape[0] <= nlev + 1:
            return ds[:]
        else:
            raise RuntimeError(
                f"The field `{name}` is expected to have shape ({nlon}(+1),) or "
                f"({nlev}(+1),), but has shape {ds.shape}."
            )

    def _get_field_2d(self, ds, name):
        nlon = self.get_nlon()
        nlev = self.get_nlev()
        if nlon <= ds.shape[0] <= nlon + 1 and nlev <= ds.shape[1] <= nlev + 1:
            return ds[...]
        elif nlon <= ds.shape[1] <= nlon + 1 and nlev <= ds.shape[0] <= nlev + 1:
            return np.transpose(ds[...])
        else:
            raise RuntimeError(
                f"The field `{name}` is expected to have shape "
                f"({nlon}(+1), {nlev}(+1)) or ({nlev}(+1), {nlon}(+1)), "
                f"but has shape {ds.shape}."
            )

    def _get_field_3d(self, ds, name):
        nlon = self.get_nlon()
        nlev = self.get_nlev()

        if nlon in ds.shape:
            axes = [ds.shape.index(nlon)]
        elif nlon + 1 in ds.shape:
            axes = [ds.shape.index(nlon + 1)]
        else:
            raise RuntimeError(f"The field `{name}` has unexpected shape {ds.shape}.")

        if nlev in ds.shape:
            axes += [ds.shape.index(nlev)]
        elif nlev + 1 in ds.shape:
            axes += [ds.shape.index(nlev + 1)]
        else:
            raise RuntimeError(f"The field `{name}` has unexpected shape {ds.shape}.")

        axes += tuple({0, 1, 2} - set(axes))

        return np.transpose(ds[...], axes=axes)

    def _initialize_parameters(
        self,
        parameter_cls: Type[BaseModel],
        get_parameter_name: Optional[Callable[[str], str]] = None,
    ):
        init_dict = {}
        for attr_name, metadata in parameter_cls.schema()["properties"].items():
            param_name = (
                get_parameter_name(attr_name) if get_parameter_name is not None else attr_name
            )
            param_type = metadata["type"]
            if param_type == "boolean":
                init_dict[attr_name] = self._get_parameter_b(param_name)
            elif param_type == "number":
                init_dict[attr_name] = self._get_parameter_f(param_name)
            elif param_type == "integer":
                init_dict[attr_name] = self._get_parameter_i(param_name)
            else:
                raise ValueError(f"Invalid parameter type `{param_type}`.")
        return parameter_cls(**init_dict)

    def _get_parameter_b(self, name: str) -> bool:
        return self.data_types.bool(self.f.get(name, [True])[0])

    def _get_parameter_f(self, name: str) -> float:
        return self.data_types.float(self.f.get(name, [0.0])[0])

    def _get_parameter_i(self, name: str) -> int:
        return self.data_types.int(self.f.get(name, [0])[0])
