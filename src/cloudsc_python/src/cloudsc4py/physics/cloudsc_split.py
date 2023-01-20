# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from __future__ import annotations
from functools import cached_property
from itertools import repeat
import numpy as np
import sys
from typing import TYPE_CHECKING

from cloudsc4py.framework.components import ImplicitTendencyComponent
from cloudsc4py.framework.grid import I, J, K
from cloudsc4py.framework.storage import managed_temporary_storage
from cloudsc4py.utils.numpyx import assign

if TYPE_CHECKING:
    from datetime import timedelta
    from typing import Dict

    from sympl._core.typingx import PropertyDict

    from cloudsc4py.framework.config import GT4PyConfig
    from cloudsc4py.framework.grid import ComputationalGrid
    from cloudsc4py.utils.iox import (
        YoecldpParameters,
        YoethfParameters,
        YomcstParameters,
        YrecldpParameters,
    )
    from cloudsc4py.utils.typingx import StorageDict


class Cloudsc(ImplicitTendencyComponent):
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        yoecldp_parameters: YoecldpParameters,
        yoethf_parameters: YoethfParameters,
        yomcst_parameters: YomcstParameters,
        yrecldp_parameters: YrecldpParameters,
        *,
        enable_checks: bool = True,
        gt4py_config: GT4PyConfig,
    ) -> None:
        super().__init__(computational_grid, enable_checks=enable_checks, gt4py_config=gt4py_config)

        self.nlev = self.computational_grid.grids[I, J, K].shape[2]
        externals = {}
        externals.update(yoecldp_parameters.dict())
        externals.update(yoethf_parameters.dict())
        externals.update(yomcst_parameters.dict())
        externals.update(yrecldp_parameters.dict())
        externals.update(
            {
                "DEPICE": 1,
                "EPSEC": 1e-14,
                "EPSILON": 100 * sys.float_info.epsilon,
                "EVAPRAIN": 2,
                "EVAPSNOW": 1,
                "FALLQV": False,
                "FALLQL": False,
                "FALLQI": False,
                "FALLQR": True,
                "FALLQS": True,
                "MELTQV": -99,
                "MELTQL": yoecldp_parameters.NCLDQI,
                "MELTQI": yoecldp_parameters.NCLDQR,
                "MELTQR": yoecldp_parameters.NCLDQS,
                "MELTQS": yoecldp_parameters.NCLDQR,
                "NLEV": self.nlev,
                "PHASEQV": 0,
                "PHASEQL": 1,
                "PHASEQI": 2,
                "PHASEQR": 1,
                "PHASEQS": 2,
                "RDCP": yomcst_parameters.RD / yomcst_parameters.RCPD,
                "RLDCP": 1 / (yoethf_parameters.RALSDCP - yoethf_parameters.RALVDCP),
                "TW1": 1329.31,
                "TW2": 0.0074615,
                "TW3": 0.85e5,
                "TW4": 40.637,
                "TW5": 275.0,
                "VQV": 0.0,
                "VQL": 0.0,
                "VQI": yrecldp_parameters.RVICE,
                "VQR": yrecldp_parameters.RVRAIN,
                "VQS": yrecldp_parameters.RVSNOW,
                "WARMRAIN": 2,
            }
        )

        self.cloudsc_tendencies = self.compile_stencil("cloudsc_tendencies", externals)
        self.cloudsc_fluxes = self.compile_stencil("cloudsc_fluxes", externals)

    @cached_property
    def _input_properties(self) -> PropertyDict:
        # todo(stubbiali): sort out units
        return {
            "b_convection_on": {"grid": (I, J), "units": ""},
            "f_a": {"grid": (I, J, K), "units": ""},
            "f_ap": {"grid": (I, J, K), "units": ""},
            "f_aph": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_ccn": {"grid": (I, J, K), "units": ""},
            "f_hrlw": {"grid": (I, J, K), "units": ""},
            "f_hrsw": {"grid": (I, J, K), "units": ""},
            "f_icrit_aer": {"grid": (I, J, K), "units": ""},
            "f_lcrit_aer": {"grid": (I, J, K), "units": ""},
            "f_lsm": {"grid": (I, J), "units": ""},
            "f_lu": {"grid": (I, J, K), "units": ""},
            "f_lude": {"grid": (I, J, K), "units": ""},
            "f_mfd": {"grid": (I, J, K), "units": ""},
            "f_mfu": {"grid": (I, J, K), "units": ""},
            "f_nice": {"grid": (I, J, K), "units": ""},
            "f_qi": {"grid": (I, J, K), "units": ""},
            "f_ql": {"grid": (I, J, K), "units": ""},
            "f_qr": {"grid": (I, J, K), "units": ""},
            "f_qs": {"grid": (I, J, K), "units": ""},
            "f_qv": {"grid": (I, J, K), "units": ""},
            "f_re_ice": {"grid": (I, J, K), "units": ""},
            "f_snde": {"grid": (I, J, K), "units": ""},
            "f_supsat": {"grid": (I, J, K), "units": ""},
            "f_t": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_a": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_qi": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_ql": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_qr": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_qs": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_qv": {"grid": (I, J, K), "units": ""},
            "f_tnd_tmp_t": {"grid": (I, J, K), "units": ""},
            "f_vfi": {"grid": (I, J, K), "units": ""},
            "f_vfl": {"grid": (I, J, K), "units": ""},
            "f_w": {"grid": (I, J, K), "units": ""},
            "i_convection_type": {"grid": (I, J), "units": ""},
        }

    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        # todo(stubbiali): sort out units
        return {
            "f_a": {"grid": (I, J, K), "units": "s^-1"},
            "f_t": {"grid": (I, J, K), "units": "s^-1"},
            "f_qv": {"grid": (I, J, K), "units": "s^-1"},
            "f_ql": {"grid": (I, J, K), "units": "s^-1"},
            "f_qi": {"grid": (I, J, K), "units": "s^-1"},
            "f_qr": {"grid": (I, J, K), "units": "s^-1"},
            "f_qs": {"grid": (I, J, K), "units": "s^-1"},
        }

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        # todo(stubbiali): sort out units
        return {
            "f_covptot": {"grid": (I, J, K), "units": ""},
            "f_fcqlng": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fcqnng": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fcqrng": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fcqsng": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fhpsl": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fhpsn": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fplsl": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fplsn": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqif": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqitur": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqlf": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqltur": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqrf": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_fsqsf": {"grid": (I, J, K - 1 / 2), "units": ""},
            "f_rainfrac_toprfz": {"grid": (I, J), "units": ""},
        }

    def array_call(
        self,
        state: StorageDict,
        timestep: timedelta,
        out_tendencies: StorageDict,
        out_diagnostics: StorageDict,
        overwrite_tendencies: Dict[str, bool],
    ) -> None:
        with managed_temporary_storage(
            self.computational_grid,
            *repeat(((I, J), "float"), 6),
            ((I, J), "bool"),
            ((K,), "int"),
            *repeat(((I, J, K), "float"), 18),
            gt4py_config=self.gt4py_config,
        ) as (
            aph_s,
            cldtopdist,
            covpmax,
            covptot,
            paphd,
            trpaus,
            rainliq,
            klevel,
            foealfa,
            lneg_qi,
            lneg_ql,
            lneg_qr,
            lneg_qs,
            lude,
            pfplsi,
            pfplsl,
            pfplsr,
            pfplss,
            qi0,
            qin,
            ql0,
            qln,
            qr0,
            qrn,
            qs0,
            qsn,
        ):
            inputs = {
                "in_" + name.split("_", maxsplit=1)[1]: state[name]
                for name in self.input_properties
            }
            tendencies = {
                "out_tnd_loc_" + name.split("_", maxsplit=1)[1]: out_tendencies[name]
                for name in self.tendency_properties
            }
            diagnostics = {
                "out_" + name.split("_", maxsplit=1)[1]: out_diagnostics[name]
                for name in self.diagnostic_properties
            }
            temporaries = {
                "tmp_aph_s": aph_s,
                "tmp_cldtopdist": cldtopdist,
                "tmp_covpmax": covpmax,
                "tmp_covptot": covptot,
                "tmp_klevel": klevel,
                "tmp_paphd": paphd,
                "tmp_rainliq": rainliq,
                "tmp_trpaus": trpaus,
            }
            aph_s[...] = state["f_aph"][..., self.nlev]
            assign(klevel, np.arange(self.nlev + 1))

            inputs1 = inputs.copy()
            vfi = inputs1.pop("in_vfi")
            vfl = inputs1.pop("in_vfl")
            diagnostics1 = {
                "out_covptot": diagnostics["out_covptot"],
                "out_foealfa": foealfa,
                "out_lneg_qi": lneg_qi,
                "out_lneg_ql": lneg_ql,
                "out_lneg_qr": lneg_qr,
                "out_lneg_qs": lneg_qs,
                "out_lude": lude,
                "out_pfplsi": pfplsi,
                "out_pfplsl": pfplsl,
                "out_pfplsr": pfplsr,
                "out_pfplss": pfplss,
                "out_qi0": qi0,
                "out_qin": qin,
                "out_ql0": ql0,
                "out_qln": qln,
                "out_qr0": qr0,
                "out_qrn": qrn,
                "out_qs0": qs0,
                "out_qsn": qsn,
                "out_rainfrac_toprfz": diagnostics["out_rainfrac_toprfz"],
            }
            self.cloudsc_tendencies(
                **inputs1,
                **tendencies,
                **diagnostics1,
                **temporaries,
                dt=timestep.total_seconds(),
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )

            inputs2 = {
                "in_aph": inputs["in_aph"],
                "in_foealfa": foealfa,
                "in_lneg_qi": lneg_qi,
                "in_lneg_ql": lneg_ql,
                "in_lneg_qr": lneg_qr,
                "in_lneg_qs": lneg_qs,
                "in_lude": lude,
                "in_pfplsi": pfplsi,
                "in_pfplsl": pfplsl,
                "in_pfplsr": pfplsr,
                "in_pfplss": pfplss,
                "in_qi0": qi0,
                "in_qin": qin,
                "in_ql0": ql0,
                "in_qln": qln,
                "in_qr0": qr0,
                "in_qrn": qrn,
                "in_qs0": qs0,
                "in_qsn": qsn,
                "in_vfi": vfi,
                "in_vfl": vfl,
            }
            outputs2 = diagnostics.copy()
            outputs2.pop("out_covptot")
            outputs2.pop("out_rainfrac_toprfz")
            self.cloudsc_fluxes(
                **inputs2,
                **outputs2,
                dt=timestep.total_seconds(),
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K - 1 / 2].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )
