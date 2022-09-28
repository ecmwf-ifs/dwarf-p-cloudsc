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
    DefaultLine("gt:cpu_kfirst", "cyan", "-", ">"),
    DefaultLine("gt:cpu_ifirst", "blue", "-", "<"),
    # DefaultLine("gt:gpu", "green", "-", "^"),
    DefaultLine("cuda", "purple", "-", "o"),
    DefaultLine("dace:gpu", "orange", "-", "s"),
]
line_pool = LinePool("performance.csv", lines, reference_backend="fortran")
nx_l = tuple(2**i for i in range(10, 18))


figure_properties = {
    "figsize": [6, 6],
    "fontsize": 16,
    "tight_layout": True,
}
axes_properties = {
    "fontsize": 16,
    "x_label": "Number of columns [-]",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": tuple(f"$2^{{{int(np.log2(nx))}}}$" for nx in nx_l),
    "y_label": "Speed-up w.r.t. FORTRAN [-]",
    "y_lim": [0, 8],
    "y_ticks": range(0, 9),
    "y_tick_labels": None,
    "legend_on": True,
    "legend_fontsize": 14,
    "legend_framealpha": 1.0,
    "legend_loc": "upper left",
    "legend_ncol": 1,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
}


def plot_lines(ax: plt.Axes, line_pool: LinePool) -> dict[str, list[float]]:
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


def main():
    fig, ax = get_figure_and_axes(nrows=1, ncols=1, index=1, **figure_properties)
    plot_lines(ax, line_pool)
    set_axes_properties(ax, **axes_properties)
    set_figure_properties(fig, **figure_properties)
    plt.show()


if __name__ == "__main__":
    main()
