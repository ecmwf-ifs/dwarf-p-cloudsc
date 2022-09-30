# -*- coding: utf-8 -*-
from __future__ import annotations
import dataclasses
from functools import partial
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

from plot_utils import get_figure_and_axes, set_axes_properties, set_figure_properties

if TYPE_CHECKING:
    from typing import Optional


@dataclasses.dataclass
class Line:
    backend: str
    color: str
    linestyle: str
    marker: str
    linewidth: float
    markersize: float
    markerfacecolor: Optional[str] = None
    markeredgecolor: Optional[str] = None
    label: Optional[str] = None

    def __post_init__(self):
        self.markerfacecolor = self.markerfacecolor or self.color
        self.markeredgecolor = self.markeredgecolor or self.color
        self.label = self.label or self.backend


DEFAULT_LINEWIDTH = 1.5
DEFAULT_MARKERSIZE = 10
DefaultLine = partial(Line, linewidth=DEFAULT_LINEWIDTH, markersize=DEFAULT_MARKERSIZE, label=None)


@dataclasses.dataclass
class LinePool:
    output_file: str
    lines: list[Line]
    reference_backend: str = "fortran"


lines = [
    # DefaultLine("fortran", "black", "-", None),
    DefaultLine("gpu-scc", "magenta", "-", "p", label="FORTRAN: gpu-scc"),
    DefaultLine("gpu-scc-hoist", "pink", "-", "H", label="FORTRAN: gpu-scc-hoist"),
    DefaultLine("gt:cpu_kfirst", "cyan", "-", ">", label="Python: gt:cpu_kfirst"),
    DefaultLine("gt:cpu_ifirst", "blue", "-", "<", label="Python: gt:cpu_ifirst"),
    # DefaultLine("gt:gpu", "green", "-", "^", label="Python: gt:gpu"),
    DefaultLine("cuda", "purple", "-", "o", label="Python: cuda"),
    DefaultLine("dace:gpu", "orange", "-", "s", label="Python: dace:gpu"),
]
line_pool_00 = LinePool("../drivers/performance_cray.csv", lines)
line_pool_01 = LinePool("../drivers/performance_gnu.csv", lines)
line_pool_10 = LinePool("../drivers/performance_intel.csv", lines)
line_pool_11 = LinePool("../drivers/performance_nvidia.csv", lines)
nx_l = tuple(2**i for i in range(10, 18))


figure_properties = {
    "figsize": [12.5, 8],
    "fontsize": 16,
    "tight_layout": True,
    "figlegend_on": True,
    "figlegend_ax": 0,
    "figlegend_loc": "center right",
    "subplots_adjust_right": 0.7,
}
axes_properties_00 = {
    "fontsize": 16,
    "x_label": "",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": (),
    "y_label": "Speed-up w.r.t FORTRAN [-]",
    "y_scale": "symlog",
    "y_scale_kwargs": {"base": 2, "linthresh": 1},
    "y_lim": (0, 12),
    "y_ticks": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "y_tick_labels": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "legend_on": False,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
    "text": "$\\mathbf{a)}$ Cray",
    "text_loc": "upper left",
}
axes_properties_01 = {
    "fontsize": 16,
    "x_label": "",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": (),
    "y_label": "",
    "y_scale": "symlog",
    "y_scale_kwargs": {"base": 2, "linthresh": 1},
    "y_lim": (0, 12),
    "y_ticks": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "y_tick_labels": (),
    "legend_on": False,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
    "text": "$\\mathbf{b)}$ GNU",
    "text_loc": "upper left",
}
axes_properties_10 = {
    "fontsize": 16,
    "x_label": "Number of columns [-]",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": tuple(f"$2^{{{int(np.log2(nx))}}}$" for nx in nx_l),
    "y_label": "Speed-up w.r.t. FORTRAN [-]",
    "y_scale": "symlog",
    "y_scale_kwargs": {"base": 2, "linthresh": 1},
    "y_lim": (0, 12),
    "y_ticks": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "y_tick_labels": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "legend_on": False,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
    "text": "$\\mathbf{c)}$ Intel",
    "text_loc": "upper left",
}
axes_properties_11 = {
    "fontsize": 16,
    "x_label": "Number of columns [-]",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": tuple(f"$2^{{{int(np.log2(nx))}}}$" for nx in nx_l),
    "y_label": "",
    "y_scale": "symlog",
    "y_scale_kwargs": {"base": 2, "linthresh": 1},
    "y_lim": (0, 12),
    "y_ticks": (0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12),
    "y_tick_labels": (),
    "legend_on": False,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
    "text": "$\\mathbf{d)}$ NVIDIA",
    "text_loc": "upper left",
}


def plot_speedup(ax: plt.Axes, line_pool: LinePool) -> dict[str, list[float]]:
    speedups = {line.backend: [] for line in line_pool.lines}
    df = pd.read_csv(line_pool.output_file)
    df_ref = df[df.backend == line_pool.reference_backend]

    # collect data
    for nx in nx_l:
        df_ref_1 = df_ref[df_ref.nx == nx]
        runtime_ref = df_ref_1.loc[df_ref_1.index[0], "mean"]

        for line in line_pool.lines:
            df1 = df[df.backend == line.backend]
            df2 = df1[df1.nx == nx]
            if df2.empty:
                speedups[line.backend].append(math.nan)
            else:
                runtime = df2.loc[df2.index[0], "mean"]
                speedups[line.backend].append(runtime_ref / runtime)

    # plot
    ax.plot(nx_l, (1,) * len(nx_l), "k-", linewidth=1)
    for line in line_pool.lines:
        ax.plot(
            nx_l,
            speedups[line.backend],
            line.color,
            linestyle=line.linestyle,
            linewidth=line.linewidth,
            marker=line.marker,
            markersize=line.markersize,
            markerfacecolor=line.markerfacecolor,
            markeredgecolor=line.markeredgecolor,
            label=line.label,
        )

    return speedups


def plot_runtime(ax: plt.Axes, line_pool: LinePool) -> dict[str, list[float]]:
    runtimes = {line.backend: [] for line in line_pool.lines}
    df = pd.read_csv(line_pool.output_file)

    # collect data
    for nx in nx_l:
        for line in line_pool.lines:
            df1 = df[df.backend == line.backend]
            df2 = df1[df1.nx == nx]
            if df2.empty:
                runtimes[line.backend].append(math.nan)
            else:
                runtimes[line.backend].append(df2.loc[df2.index[0], "mean"])

    # plot
    for line in line_pool.lines:
        ax.plot(
            nx_l,
            runtimes[line.backend],
            line.color,
            linestyle=line.linestyle,
            linewidth=line.linewidth,
            marker=line.marker,
            markersize=line.markersize,
            markerfacecolor=line.markerfacecolor,
            markeredgecolor=line.markeredgecolor,
            label=line.label,
        )

    return runtimes


def main():
    plot_lines = plot_speedup

    fig, ax_00 = get_figure_and_axes(nrows=2, ncols=2, index=1, **figure_properties)
    plot_lines(ax_00, line_pool_00)
    set_axes_properties(ax_00, **axes_properties_00)

    _, ax_01 = get_figure_and_axes(fig=fig, nrows=2, ncols=2, index=2, **figure_properties)
    plot_lines(ax_01, line_pool_01)
    set_axes_properties(ax_01, **axes_properties_01)

    _, ax_10 = get_figure_and_axes(fig=fig, nrows=2, ncols=2, index=3, **figure_properties)
    plot_lines(ax_10, line_pool_10)
    set_axes_properties(ax_10, **axes_properties_10)

    _, ax_11 = get_figure_and_axes(fig=fig, nrows=2, ncols=2, index=4, **figure_properties)
    plot_lines(ax_11, line_pool_11)
    set_axes_properties(ax_11, **axes_properties_11)

    # _, ax_12 = get_figure_and_axes(fig=fig, nrows=2, ncols=2, index=3, **figure_properties)
    # handles, labels = ax_00.get_legend_handles_labels()
    # ax_01.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=ax_12.transAxes)
    # # ax_12.legend(handles, labels, borderaxespad=0)
    # ax_12.axis("off")

    set_figure_properties(fig, **figure_properties)
    plt.show()


if __name__ == "__main__":
    main()
