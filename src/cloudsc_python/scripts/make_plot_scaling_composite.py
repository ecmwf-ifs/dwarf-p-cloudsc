# -*- coding: utf-8 -*-
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np

from make_plot_scaling import DefaultLine, LinePool, plot_lines
from plot_utils import get_figure_and_axes, set_axes_properties, set_figure_properties


# >>> config: start
lines = [
    # DefaultLine("fortran", "black", "-", None),
    DefaultLine("gpu-scc", "magenta", "-", "p", label="FORTRAN: gpu-scc"),
    DefaultLine("gpu-scc-hoist", "pink", "-", "H", label="FORTRAN: gpu-scc-hoist"),
    DefaultLine("gt:cpu_kfirst", "cyan", "-", ">", label="Python: gt:cpu_kfirst"),
    DefaultLine("gt:cpu_ifirst", "blue", "-", "<", label="Python: gt:cpu_ifirst"),
    DefaultLine("gt:gpu", "green", "-", "^", label="Python: gt:gpu"),
    DefaultLine("cuda", "purple", "-", "o", label="Python: cuda"),
    DefaultLine("dace:gpu", "orange", "-", "s", label="Python: dace:gpu"),
]
line_pool_00 = LinePool("../data/performance_cray.csv", lines, host="dom")
line_pool_01 = LinePool("../data/performance_gnu.csv", lines, host="dom")
line_pool_10 = LinePool("../data/performance_intel.csv", lines, host="dom")
line_pool_11 = LinePool("../data/performance_nvidia.csv", lines, host="dom")
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
    "y_label": "Speed-up w.r.t. FORTRAN [-]",
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
# >>> config: end


def main():
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

    set_figure_properties(fig, **figure_properties)
    plt.show()


if __name__ == "__main__":
    main()
