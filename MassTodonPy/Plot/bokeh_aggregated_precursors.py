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

from MassTodonPy.Misc.io import create_folder_if_needed

def bokeh_aggregated_precursors(masstodon,
                                path='aggregated_precusors.html',
                                mode='inline',
                                show=True,
                                width=400,
                                height=400,
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
    plot_width = kwds.get('')
    intensities = masstodon.report.get_aggregated_precursors()
    create_folder_if_needed(path)
    output_file(path, mode=mode)
    charges = list('q = {0}'.format(x) for x in range(1, 1+len(intensities)))
    p = figure(plot_width=int(width),
               plot_height=int(height),
               x_range=charges)
    p.xaxis.axis_label = 'Charge'
    p.yaxis.axis_label ='Estimated Intensity'
    data = ColumnDataSource({'intensity': intensities,
                             'charge': charges})
    bars = p.vbar(source=data, x='charge', top='intensity', width=1)
    hover_bars = HoverTool(renderers=[bars],
                           tooltips=[('charge', '@charge'),
                                     ('intensity', "@intensity{0,0}")],
                           mode='vline')
    p.add_tools(hover_bars)
    if show:
        show_plot(p)
    else:
        with open(path, 'w') as f:
            f.write(file_html(p, CDN, path))
    return p
