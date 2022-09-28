# -*- coding: utf-8 -*-
from matplotlib import rcParams
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numbers
from typing import Optional


text_locations = {
    "upper right": 1,
    "upper left": 2,
    "lower left": 3,
    "lower right": 4,
    "right": 5,
    "center left": 6,
    "center right": 7,
    "lower center": 8,
    "upper center": 9,
    "center": 10,
}


linestyle_dict = {
    "-": "solid",
    "--": "dashed",
    ":": "dotted",
    "-.": "dashdot",
    "solid": "solid",
    "dashed": "dashed",
    "dotted": "dotted",
    "dashdot": "dashdot",
    "loosely dotted": (0, (1, 10)),
    "densely dotted": (0, (1, 1)),
    "loosely dashed": (0, (5, 10)),
    "densely dashed": (0, (5, 1)),
    "loosely dashdotted": (0, (3, 10, 1, 10)),
    "dashdotted": (0, (3, 5, 1, 5)),
    "densely dashdotted": (0, (3, 1, 1, 1)),
    "dashdotdotted": (0, (3, 5, 1, 5, 1, 5)),
    "loosely dashdotdotted": (0, (3, 10, 1, 10, 1, 10)),
    "densely dashdotdotted": (0, (3, 1, 1, 1, 1, 1)),
}


def get_figure_and_axes(
    fig: Optional[plt.Figure] = None,
    ax: Optional[plt.Axes] = None,
    default_fig: Optional[plt.Figure] = None,
    nrows: int = 1,
    ncols: int = 1,
    index: int = 1,
    **kwargs,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Get a :class:`matplotlib.figure.Figure` and a :class:`matplotlib.axes.Axes`, with the latter
    embedded in the former.
    """
    figsize = kwargs.get("figsize", (7, 7))
    fontsize = kwargs.get("fontsize", 12)
    projection = kwargs.get("projection", None)

    rcParams["font.size"] = fontsize

    if (fig is not None) and (ax is not None):
        try:
            if ax not in fig.get_axes():
                import warnings

                warnings.warn(
                    "Input axes do not belong to the input figure, "
                    "so the figure which the axes belong to is considered instead.",
                    RuntimeWarning,
                )

                out_fig, out_ax = ax.get_figure(), ax
            else:
                out_fig, out_ax = fig, ax
        except AttributeError:
            import warnings

            warnings.warn(
                "Input argument ``fig`` does not seem to be a matplotlib.figure.Figure, "
                "so the figure the axes belong to are considered.",
                RuntimeWarning,
            )

            out_fig, out_ax = ax.get_figure(), ax
    elif (fig is not None) and (ax is None):
        try:
            out_fig = fig
            out_ax = out_fig.add_subplot(nrows, ncols, index, projection=projection)
        except AttributeError:
            import warnings

            warnings.warn(
                "Input argument ``fig`` does not seem to be a matplotlib.figure.Figure, "
                "hence a proper matplotlib.figure.Figure object is created.",
                RuntimeWarning,
            )

            out_fig = plt.figure(figsize=figsize)
            out_ax = out_fig.add_subplot(nrows, ncols, index, projection=projection)
    elif (fig is None) and (ax is not None):
        out_fig, out_ax = ax.get_figure(), ax
    else:  # (fig is None) and (ax is None)
        if default_fig is None:
            out_fig = plt.figure(figsize=figsize)
            out_ax = out_fig.add_subplot(nrows, ncols, index, projection=projection)
        else:
            try:
                out_fig = default_fig
                out_ax = out_fig.add_subplot(nrows, ncols, index, projection=projection)
            except AttributeError:
                import warnings

                warnings.warn(
                    "Input argument ``default_fig`` does not seem to be a "
                    "matplotlib.figure.Figure, hence a proper matplotlib.figure.Figure "
                    "object is created.",
                    RuntimeWarning,
                )

                out_fig = plt.figure(figsize=figsize)
                out_ax = out_fig.add_subplot(nrows, ncols, index, projection=projection)

    return out_fig, out_ax


def set_figure_properties(fig: plt.Figure, **kwargs) -> None:
    """Ease the configuration of a :class:`matplotlib.figure.Figure`."""
    fontsize = kwargs.get("fontsize", 12)
    tight_layout = kwargs.get("tight_layout", True)
    tight_layout_rect = kwargs.get("tight_layout_rect", (0, 0, 1, 1))
    tight_layout_wpad = kwargs.get("tight_layout_wpad", None)
    tight_layout_hpad = kwargs.get("tight_layout_hpad", None)
    suptitle = kwargs.get("suptitle", "")
    x_label = kwargs.get("x_label", "")
    x_labelpad = kwargs.get("x_labelpad", 20)
    y_label = kwargs.get("y_label", "")
    y_labelpad = kwargs.get("y_labelpad", 20)
    figlegend_on = kwargs.get("figlegend_on", False)
    figlegend_multiple = kwargs.get("figlegend_multiple", False)
    figlegend_ax = kwargs.get("figlegend_ax", None)
    figlegend_framealpha = kwargs.get("figlegend_framealpha", 1.0)
    figlegend_loc = kwargs.get("figlegend_loc", "lower center")
    figlegend_ncol = kwargs.get("figlegend_ncol", 1)
    figlegend_title = kwargs.get("figlegend_title", None)
    right = kwargs.get("subplots_adjust_right", None)
    wspace = kwargs.get("subplots_adjust_wspace", None)
    hspace = kwargs.get("subplots_adjust_hspace", None)

    rcParams["font.size"] = fontsize

    if suptitle is not None and suptitle != "":
        fig.suptitle(suptitle, fontsize=fontsize + 1)

    if x_label != "" or y_label != "":
        ax = fig.add_subplot(111)
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_xticklabels([], visible=False)
        ax.set_yticks([])
        ax.set_yticklabels([], visible=False)

        if x_label != "":
            ax.set_xlabel(x_label, labelpad=x_labelpad)
        if y_label != "":
            ax.set_ylabel(y_label, labelpad=y_labelpad)

    if tight_layout:
        fig.tight_layout(rect=tight_layout_rect, w_pad=tight_layout_wpad, h_pad=tight_layout_hpad)

    fig.subplots_adjust(right=right, wspace=wspace, hspace=hspace)

    if figlegend_on:
        figlegend_ax = [figlegend_ax] if isinstance(figlegend_ax, int) else figlegend_ax
        axes = (
            [fig.get_axes()[i] for i in figlegend_ax]
            if figlegend_ax is not None
            else fig.get_axes()
        )

        title = (
            [figlegend_title] * len(axes)
            if figlegend_title is None or isinstance(figlegend_title, str)
            else figlegend_title
        )

        extra = patches.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0)

        if figlegend_multiple:
            framealpha = (
                [figlegend_framealpha] * len(axes)
                if isinstance(figlegend_framealpha, numbers.Number)
                else figlegend_framealpha
            )
            loc = [figlegend_loc] * len(axes) if isinstance(figlegend_loc, str) else figlegend_loc
            ncol = (
                [int(figlegend_ncol)] * len(axes)
                if isinstance(figlegend_ncol, numbers.Number)
                else figlegend_ncol
            )

            for i in range(len(axes)):
                handles, labels = axes[i].get_legend_handles_labels()
                fig.legend(
                    [extra, *handles] if title[i] is not None else handles,
                    [title[i], *labels] if title[i] is not None else labels,
                    framealpha=framealpha[i],
                    loc=loc[i],
                    ncol=ncol[i],
                )
        else:
            handles = []
            labels = []

            for i in range(len(axes)):
                _handles, _labels = axes[i].get_legend_handles_labels()

                if title[i] is not None:
                    handles = handles + [extra] + _handles
                    labels = labels + [title[i]] + _labels
                    if i < len(axes) - 1:
                        handles.append(extra)
                        labels.append("")
                else:
                    handles = handles + _handles
                    labels = labels + _labels

            fig.legend(
                handles,
                labels,
                framealpha=figlegend_framealpha,
                loc=figlegend_loc,
                ncol=figlegend_ncol,
            )


def set_axes_properties(ax: plt.Axes, **kwargs) -> None:
    """Ease the configuration of a :class:`matplotlib.axes.Axes` object."""
    fontsize = kwargs.get("fontsize", 12)
    # title
    title_center = kwargs.get("title_center", "")
    title_left = kwargs.get("title_left", "")
    title_right = kwargs.get("title_right", "")
    # x-axis
    x_label = kwargs.get("x_label", "")
    x_label_color = kwargs.get("x_label_color", "black")
    x_lim = kwargs.get("x_lim", None)
    x_invert = kwargs.get("x_invert", False)
    x_scale = kwargs.get("x_scale", "linear")
    x_scale_kwargs = kwargs.get("x_scale_kwargs", {})
    x_ticks = kwargs.get("x_ticks", None)
    x_tick_length = kwargs.get("x_tick_length", None)
    x_tick_labels = kwargs.get("x_tick_labels", None)
    x_tick_labels_color = kwargs.get("x_tick_labels_color", "black")
    x_tick_labels_rotation = kwargs.get("x_tick_labels_rotation", 0)
    x_tick_format = kwargs.get("x_tick_format", None)
    x_minor_ticks_visible = kwargs.get("x_minor_ticks_visible", False)
    x_visible = kwargs.get("x_visible", True)
    # y-axis
    y_label = kwargs.get("y_label", "")
    y_label_color = kwargs.get("y_label_color", "black")
    y_lim = kwargs.get("y_lim", None)
    y_invert = kwargs.get("y_invert", False)
    y_scale = kwargs.get("y_scale", "linear")
    y_scale_kwargs = kwargs.get("y_scale_kwargs", {})
    y_ticks = kwargs.get("y_ticks", None)
    y_tick_labels = kwargs.get("y_tick_labels", None)
    y_tick_labels_color = kwargs.get("y_tick_labels_color", "black")
    y_tick_labels_rotation = kwargs.get("y_tick_labels_rotation", 0)
    y_tick_format = kwargs.get("y_tick_format", None)
    y_minor_ticks_visible = kwargs.get("y_minor_ticks_visible", False)
    y_visible = kwargs.get("y_visible", True)
    # legend
    legend_on = kwargs.get("legend_on", False)
    legend_loc = kwargs.get("legend_loc", "best")
    legend_bbox_to_anchor = kwargs.get("legend_bbox_to_anchor", None)
    legend_framealpha = kwargs.get("legend_framealpha", 0.5)
    legend_ncol = kwargs.get("legend_ncol", 1)
    legend_fontsize = kwargs.get("legend_fontsize", fontsize)
    # textbox
    text = kwargs.get("text", None)
    text_loc = kwargs.get("text_loc", "")
    # grid
    grid_on = kwargs.get("grid_on", False)
    grid_properties = kwargs.get("grid_properties", None)
    # ax2 title
    ax2_title_center = kwargs.get("ax2_title_center", "")
    ax2_title_left = kwargs.get("ax2_title_left", "")
    ax2_title_right = kwargs.get("ax2_title_right", "")
    # x2-axis
    ax2_on = kwargs.get("ax2_on", False)
    x2_label = kwargs.get("x2_label", "")
    x2_label_color = kwargs.get("x2_label_color", "black")
    x2_lim = kwargs.get("x2_lim", None)
    x2_invert = kwargs.get("x2_invert", False)
    x2_scale = kwargs.get("x2_scale", None)
    x2_ticks = kwargs.get("x2_ticks", None)
    x2_tick_labels = kwargs.get("x2_tick_labels", None)
    x2_tick_labels_color = kwargs.get("x2_tick_labels_color", "black")
    x2_minor_ticks_visible = kwargs.get("x2_minor_ticks_visible", False)
    x2_visible = kwargs.get("x2_visible", True)
    # y2-axis
    y2_label = kwargs.get("y2_label", "")
    y2_label_color = kwargs.get("y2_label_color", "black")
    y2_lim = kwargs.get("y2_lim", None)
    y2_invert = kwargs.get("y2_invert", False)
    y2_scale = kwargs.get("y2_scale", None)
    y2_ticks = kwargs.get("y2_ticks", None)
    y2_tick_labels = kwargs.get("y2_tick_labels", None)
    y2_tick_labels_color = kwargs.get("y2_tick_labels_color", "black")
    y2_minor_ticks_visible = kwargs.get("y2_minor_ticks_visible", False)
    y2_visible = kwargs.get("y2_visible", True)

    rcParams["font.size"] = fontsize
    # rcParams['text.usetex'] = True

    # plot titles
    if ax.get_title(loc="center") == "":
        ax.set_title(
            title_center, loc="center", fontsize=rcParams["font.size"], pad=15 if ax2_on else 6
        )
    if ax.get_title(loc="left") == "":
        ax.set_title(
            title_left, loc="left", fontsize=rcParams["font.size"], pad=15 if ax2_on else 6
        )
    if ax.get_title(loc="right") == "":
        ax.set_title(
            title_right, loc="right", fontsize=rcParams["font.size"], pad=15 if ax2_on else 6
        )

    # axes labels
    if ax.get_xlabel() == "" and x_label is not None and x_label != "":
        ax.set(xlabel=x_label)
    if ax.get_ylabel() == "" and y_label is not None and y_label != "":
        ax.set(ylabel=y_label)

    # axes labelcolors
    if ax.get_xlabel() != "" and x_label_color != "":
        ax.xaxis.label.set_color(x_label_color)
    if ax.get_ylabel() != "" and y_label_color != "":
        ax.yaxis.label.set_color(y_label_color)

    # axes limits
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # invert the axes
    if x_invert:
        ax.invert_xaxis()
    if y_invert:
        ax.invert_yaxis()

    # axes scale
    if x_scale is not None:
        ax.set_xscale(x_scale, **x_scale_kwargs)
    if y_scale is not None:
        ax.set_yscale(y_scale, **y_scale_kwargs)

    # axes ticks
    if x_ticks is not None:
        ax.get_xaxis().set_ticks(x_ticks)
    if y_ticks is not None:
        ax.get_yaxis().set_ticks(y_ticks)

    # tick length
    if x_tick_length is not None:
        ax.get_xaxis().set_tick_params("both", length=x_tick_length)

    # axes tick labels
    if x_tick_labels is not None:
        ax.get_xaxis().set_ticklabels(x_tick_labels)
    if y_tick_labels is not None:
        ax.get_yaxis().set_ticklabels(y_tick_labels)

    # axes tick labels color
    if x_tick_labels_color != "":
        ax.tick_params(axis="x", colors=x_tick_labels_color)
    if y_tick_labels_color != "":
        ax.tick_params(axis="y", colors=y_tick_labels_color)

    # axes tick format
    if x_tick_format is not None:
        ax.xaxis.set_major_formatter(FormatStrFormatter(x_tick_format))
    if y_tick_format is not None:
        ax.yaxis.set_major_formatter(FormatStrFormatter(y_tick_format))

    # axes tick labels rotation
    plt.xticks(rotation=x_tick_labels_rotation)
    plt.yticks(rotation=y_tick_labels_rotation)

    # unlabelled axes ticks
    if not x_minor_ticks_visible:
        ax.get_xaxis().set_tick_params(which="minor", size=0)
        ax.get_xaxis().set_tick_params(which="minor", width=0)
    if not y_minor_ticks_visible:
        ax.get_yaxis().set_tick_params(which="minor", size=0)
        ax.get_yaxis().set_tick_params(which="minor", width=0)

    # axes visibility
    if not x_visible:
        ax.get_xaxis().set_visible(False)
    if not y_visible:
        ax.get_yaxis().set_visible(False)

    # legend
    if legend_on:
        if legend_bbox_to_anchor is None:
            ax.legend(
                loc=legend_loc,
                framealpha=legend_framealpha,
                ncol=legend_ncol,
                fontsize=legend_fontsize,
            )
        else:
            ax.legend(
                loc=legend_loc,
                framealpha=legend_framealpha,
                ncol=legend_ncol,
                fontsize=legend_fontsize,
                bbox_to_anchor=legend_bbox_to_anchor,
            )

    # text box
    if text is not None:
        ax.add_artist(AnchoredText(text, loc=text_locations[text_loc]))

    # plot grid
    if grid_on:
        gps = grid_properties if grid_properties is not None else {}
        ax.grid(True, **gps)

    if ax2_on:
        ax2 = ax.twiny()

        # x2-axis
        if x2_label != "":
            ax2.set_xlabel(x2_label, labelpad=8)
        if ax2.get_xlabel() != "" and x2_label_color != "":
            ax2.xaxis.label.set_color(x2_label_color)
        ax2.set_xlim(x2_lim if x2_lim is not None else ax.get_xlim())
        if x2_invert:
            ax2.invert_xaxis()
        if x2_scale is not None:
            ax2.set_xscale(x2_scale)
        if x2_ticks is not None:
            ax2.set_xticks(x2_ticks)
        if x2_tick_labels is not None:
            ax2.set_xticklabels(x2_tick_labels)
        if x2_tick_labels_color != "":
            ax2.tick_params(axis="x", colors=x2_tick_labels_color)
        if x2_minor_ticks_visible:
            ax2.get_xaxis().set_tick_params(which="minor", size=0)
            ax2.get_xaxis().set_tick_params(which="minor", width=0)
        if not x2_visible:
            ax2.get_xaxis().set_visible(False)
        # plt.ticklabel_format(axis="x", style="sci", scilimits=(-5, -5))

        # y2-axis
        if y2_label != "":
            ax2.set_ylabel(y2_label)
        if ax2.get_ylabel() != "" and y2_label_color != "":
            ax2.yaxis.label.set_color(y2_label_color)
        ax2.set_ylim(y2_lim if y2_lim is not None else ax.get_ylim())
        if y2_invert:
            ax2.invert_yaxis()
        if y2_scale is not None:
            ax2.set_yscale(y2_scale)
        if y2_ticks is not None:
            ax2.set_yticks(y2_ticks)
        if y2_tick_labels is not None:
            ax2.set_yticklabels(y2_tick_labels)
        if y2_tick_labels_color != "":
            ax2.tick_params(axis="y", colors=y2_tick_labels_color)
        if y2_minor_ticks_visible:
            ax2.get_yaxis().set_tick_params(which="minor", size=0)
            ax2.get_yaxis().set_tick_params(which="minor", width=0)
        if not y2_visible:
            ax2.get_yaxis().set_visible(False)

        # plot titles
        if ax2.get_title(loc="center") == "":
            ax2.set_title(ax2_title_center, loc="center", fontsize=rcParams["font.size"])
        if ax2.get_title(loc="left") == "":
            ax2.set_title(ax2_title_left, loc="left", fontsize=rcParams["font.size"] - 1)
        if ax2.get_title(loc="right") == "":
            ax2.set_title(ax2_title_right, loc="right", fontsize=rcParams["font.size"] - 1)


def add_annotation(ax: plt.Axes, **kwargs) -> None:
    """Add a text annotation to a plot."""
    # get keyword arguments
    fontsize = kwargs.get("fontsize", 16)
    text = kwargs.get("text", "")
    text_color = kwargs.get("text_color", "black")
    location = kwargs.get("location", (0, 0))
    horizontal_alignment = kwargs.get("horizontal_alignment", "left")
    vertical_alignment = kwargs.get("vertical_alignment", "center")

    # add annotation
    ax.annotate(
        text,
        location,
        horizontalalignment=horizontal_alignment,
        verticalalignment=vertical_alignment,
        fontsize=fontsize,
        color=text_color,
    )
