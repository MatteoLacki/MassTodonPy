# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet
from bokeh.resources import CDN
from bokeh.embed import file_html

from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Misc.io import create_folder_if_needed

def bokeh_spectrum(masstodon,
                   path="assigned_spectrum.html",
                   mode="inline",
                   show_plot=True,
                   width=800,
                   height=600,
                   **kwds):
    """Make a plot of the assigned spectrum.

    Parameters
    ==========
    path : string
        Path to where to save the output html file.
        If not provided, MassTodon will use the folder it is called from.
    mode : string
        The mode of plotting a bokeh plot.
    width : integer
        The width of the plot.
    height : integer
        The height of the plot.

    """
    create_folder_if_needed(path)
    output_file(path, mode=mode)
    D = masstodon.report.assigned_spectrum_data
    plot = figure(plot_width=width,
                  plot_height=height,
                  tools=D['tools'])
    plot.y_range.start = D['y_range_start']
    plot.xaxis.axis_label = D['x_label']
    plot.yaxis.axis_label = D['y_label']

    # The experimental data bars
    experimental_bars = plot.vbar(**D['exp_vbar'])

    # Horizontal threshold line
    if D['threshold_line']['intensity'] > 1.0:
        plot.line(*D['threshold_line']['args'],
                  **D['threshold_line']['kwds'])

    # Plotting rectangles
    source_rectangles = ColumnDataSource(D['rectangle_data'])
    fat_rectangles = plot.quad(source=source_rectangles,
                               **D['fat_rectangles'])

    slim_rectangles = plot.quad(source=source_rectangles,
                                **D['slim_rectangles'])

    hover_fat = HoverTool(renderers=[fat_rectangles],
                          tooltips=D['fat_rectangles_tooltips'])
    plot.add_tools(hover_fat)

    # plotting peak_group / local quality of peak fitting
    source_peak_groups = ColumnDataSource(D['peak_groups_data'])
    peak_group_intensities = plot.segment(source=source_peak_groups,
                                          **D['peak_groups'])

    hover_peak_groups = HoverTool(renderers=[peak_group_intensities],
                                  tooltips=D['peak_groups_tooltips'])
    plot.add_tools(hover_peak_groups)

    # Experimental Squares
    raw_spectrum = plot.square(**D['experimental_squares'])
    hover_squares = HoverTool(renderers=[raw_spectrum],
                              tooltips=D['experimental_squares_tooltips'])
    plot.add_tools(hover_squares)

    # Text labels
    source_clusters = ColumnDataSource(D['cluster_data'])
    labels = LabelSet(source=source_clusters,
                      **D['labels'])
    plot.add_layout(labels)
    if show_plot:
        show(plot)
    else:
        with open(path, 'w') as f:
            f.write(file_html(plot, CDN, path))
    return plot
