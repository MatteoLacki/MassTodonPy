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

from bokeh.plotting import ColumnDataSource, figure, output_file, show as show_plot
from bokeh.models import HoverTool, Span, LabelSet
from bokeh.resources import CDN
from bokeh.embed import file_html

from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Misc.io import create_folder_if_needed

def bokeh_fragments_intensity(masstodon,
                              path="fragments_intensity.html",
                              mode="inline",
                              show=True,
                              width=1000,
                              height=400,
                              **kwds):
    """Make a plot of the intensity of fragments.

    Parameters
    ==========
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
    create_folder_if_needed(path)
    output_file(path, mode='cdn')
    afi = masstodon.report.aggregeted_fragment_intensities()
    afi['z_minus'] = [-v for v in afi['z']]
    p = figure(plot_width=int(width),
               plot_height=int(height))
    p.xaxis.axis_label = 'Cleavage Site'
    p.yaxis.axis_label = 'Estimated Intensity'
    afi['x'] = list(range(len(afi['z'])))
    data = ColumnDataSource(afi)
    bars_c = p.vbar(source=data, x='x', top='c', width=.8)
    hover_bars_c = HoverTool(renderers=[bars_c],
                             tooltips=[('name', '@c_name'),
                                       ('intensity', "@c{0,0}")])
    p.add_tools(hover_bars_c)
    bars_z = p.vbar(source=data, x='x', top='z_minus', width=.8, color='red')
    hover_bars_z = HoverTool(renderers=[bars_z],
                             tooltips=[('name', '@z_name'),
                                       ('intensity', "@z{0,0}")])
    p.add_tools(hover_bars_z)
    # labels_c = LabelSet(source=data, x='x', y='c')
    # p.add_layout(labels_c)
    # labels_z = LabelSet(source=data, x='x', y='z')
    # p.add_layout([labels_z, labels_c])
    if show:
        show_plot(p)
    else:
        with open(path, 'w') as f:
            f.write(file_html(p, CDN, path))
    return p
