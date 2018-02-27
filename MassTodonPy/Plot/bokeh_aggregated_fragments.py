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

from bokeh.embed import file_html
from bokeh.models import HoverTool, LabelSet
from bokeh.plotting import ColumnDataSource, figure, output_file, show as show_plot
from bokeh.resources import CDN

from MassTodonPy.Misc.io import create_folder_if_needed
from MassTodonPy.Plot.Misc import aggregate_fragments


def bokeh_aggregated_fragments(masstodon,
                               path='aggregated_fragments.html',
                               mode='inline',
                               show=False,
                               width=0,
                               height=0,
                               _offset=.2,
                               **kwds):
    """Plot intensity of precursors, neglecting the quenched charge.

    Parameters
    ==========
    intensities : list
        List of estimated intensities. Indices correspond to charge - 1.
    path : string
        Path to where to save the output html file.
        If not provided, MassTodon will use the folder it is called from.
    mode : string
        The mode of plotting a bokeh plot.
    show : boolean
        Should we show the plot in your default browser?
    plot_width : integer
        The width of the plot.
    plot_height : integer
        The height of the plot.

    """
    data = aggregate_fragments(report=masstodon.report,
                               fasta=masstodon.precursor.fasta,
                               offset=_offset)
    source = ColumnDataSource(data=data)
    create_folder_if_needed(path)
    output_file(path, mode=mode)
    if not width or not height:
        p = figure()
    else:
        p = figure(plot_width=int(width),
                   plot_height=int(height))
    p.xaxis.visible = False
    p.yaxis.axis_label = "Intensity"
    left_bars = p.vbar(x='x_left',
                       width=2*_offset,
                       top='left',
                       color="orange",
                       source=source)
    right_bars= p.vbar(x='x_right',
                       width=2*_offset,
                       top='right',
                       color="navy",
                       source=source)
    left_labels = LabelSet(x='x_left', y='left', text='left_name',
                           level='glyph',
                           x_offset=-10,
                           source=source, render_mode='css')
    right_labels = LabelSet(x='x_right', y='right', text='right_name',
                            level='glyph',
                            x_offset=-5,
                            source=source, render_mode='css')
    aa_labels = LabelSet(x='x', y=0.0, text='aa',
                         level='glyph', y_offset=-20,
                         # x_offset=-10,
                         source=source, render_mode='css')
    p.add_layout(left_labels)
    p.add_layout(right_labels)
    p.add_layout(aa_labels)
    hover_bars_l = HoverTool(renderers=[left_bars],
                             mode='vline',
                             tooltips=[('name', '@left_name'),
                                       ('intensity', "@left{0.}")])
    p.add_tools(hover_bars_l)
    hover_bars_r = HoverTool(renderers=[right_bars],
                             mode='vline',
                             tooltips=[('fragment', '@right_name'),
                                       ('intensity', "@right{0.}")])
    p.add_tools(hover_bars_r)
    if show:
        show_plot(p)
    else:
        with open(path, 'w') as f:
            f.write(file_html(p, CDN, path))
    return p
