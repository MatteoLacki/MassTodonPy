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
from bokeh.models import HoverTool, Span

from MassTodonPy.Reporter.buffers import buffers
from MassTodonPy.Reporter.misc import Brick, Cluster, ColorGenerator, make_string_represenation


class Reporter(object):
    """Report MassTodon results."""

    def __init__(self, masstodon, max_buffer_len=0.5):
        self.masstodon = masstodon
        self.base_width = 10**(-self.masstodon.mz_digits)
        self.bricks = []
        self.clusters = []
        for name, thing in self.get_bricks_clusters():
            self.__dict__[name].append(thing)

        # order clusters by m/z ratios
        self.clusters.sort(key=lambda c: c.mz_L)

        # adding left and right buffers needed for plot
        mz_L, mz_R = list(zip(*((c.mz_L, c.mz_R) for c in self.clusters))) #TODO this is an overkill: smooth out the defintion of the buffers function.
        mz_lefts, mz_rights = buffers(mz_L, mz_R, max_length=max_buffer_len)
        for c, mz_left, mz_right in zip(self.clusters, mz_lefts, mz_rights):
            c.mz_left = mz_left
            c.mz_right = mz_right

        # lexicographically order the bricks
        self.bricks.sort(key=lambda b: (b.cluster.mz_L, -b.molecule.intensity))
        # adjust brick top/bottom/color
        self.color = ColorGenerator()
        prev_b = Brick(top=0.0, bottom=0.0, cluster=Cluster(mz_L=0.0)) # guardian
        for b in self.bricks:
            if b.cluster.mz_L is prev_b.cluster.mz_L:
                b.bottom = prev_b.top
                b.top = prev_b.top + b.intensity
            b.color = self.color(b.molecule)
            prev_b = b

    def get_bricks_clusters(self):
        """Generate an unorderd flow of bricks and clusters for the plot."""
        for sol in self.masstodon._solutions:
            for G in sol:
                if G[0] is 'G':
                    try:
                        min_mz = sol.node[G]['min_mz'] - self.base_width/2
                        max_mz = sol.node[G]['max_mz'] + self.base_width/2
                    except KeyError:
                        mz = sol.node[G]['mz']
                        min_mz = mz - self.base_width/2
                        max_mz = mz + self.base_width/2
                    intensity = sol.node[G]['intensity']
                    if intensity > 0:
                        intensity_d = intensity/(max_mz - min_mz)
                        estimate = sol.node[G]['estimate']
                        estimate_d = estimate/(max_mz - min_mz)
                        cluster = Cluster(mz_L=min_mz,
                                          mz_R=max_mz,
                                          mz_left=min_mz,     # updated later
                                          mz_right=max_mz,    # updated later
                                          intensity=intensity,
                                          intensity_d=intensity_d,
                                          estimate=estimate,
                                          estimate_d=estimate_d)
                        yield 'clusters', cluster
                        for I in sol[G]:
                            estimate = sol[G][I]['estimate']
                            for M in sol[I]:
                                if M[0] is 'M':
                                    mol = sol.node[M]['molecule']
                                    yield ('bricks',
                                            Brick(cluster=cluster,
                                                  top=estimate,
                                                  bottom=0.0,
                                                  color='#440154',# updated later
                                                  intensity=estimate,
                                                  molecule=mol))

    def aggregate(self):
        """Aggregate the results of MassTodon."""
        pass

    def write(self, file_path):
        """Write results to a file.

        Parameters
        ==========
        file_path : str
            Path to the output you want MassTodon to save its results.
            Supported path extensions are .txt and .csv.
        """
        pass

    def plot(self,
             file_path="assigned_spectrum.html",
             mode="inline",
             width=None,
             height=None,
             _mult=1):
        """Make a plot.

        Parameters
        ==========
        file_path : string
            Path to the output html file.
            If not provided, MassTodon will use the default temporary location
            used by the Bokeh module.
        mode : string
            The mode of plotting a bokeh plot.
        width : integer
            The width of the plot.
        height : integer
            The height of the plot.
        """
        output_file(file_path, mode=mode)

        if not width:
            width = 800 * _mult
        if not height:
            height = 400 * _mult

        # set up plot
        TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select save".split(" ")
        plot = figure(plot_width=width,
                      plot_height=height,
                      tools=TOOLS)
        plot.y_range.start = 0
        plot.xaxis.axis_label = 'mass/charge'
        plot.yaxis.axis_label = 'intensity'

        # add bars to experimentally observed peaks / overlay real spectrum
        mz_representation = make_string_represenation('mz',
                                                      self.masstodon.mz_digits)
        experimental_bars = plot.vbar(x=self.masstodon.spectrum.mz,
                                      top=self.masstodon.spectrum.intensity,
                                      width=self.base_width,
                                      color='black',
                                      alpha=.1)
        # Horizontal threshold line
        min_mz = max(min(self.masstodon.spectrum.mz) - 50, 0)
        max_mz = max(self.masstodon.spectrum.mz) + 50
        intensity_threshold = self.masstodon.spectrum.min_intensity

        plot.line([min_mz, max_mz],
                  [intensity_threshold,intensity_threshold],
                  line_width=2,
                  color = 'red')

        # Plotting the bricks / divisions of estimated peaks into constituents
        lists = zip(*((b.cluster.mz_L,
                       b.cluster.mz_R,
                       b.cluster.mz_left,
                       b.cluster.mz_right,
                       (b.cluster.mz_L + b.cluster.mz_R)/2.0,
                       b.bottom,
                       b.top,
                       b.color,
                       b.intensity,
                       b.molecule.name,
                       b.molecule.source,
                       str(b.molecule.formula),
                       b.molecule.q,
                       b.molecule.g,
                       b.molecule.intensity)
                      for b in self.bricks))

        data = dict(zip(('mz_L',
                         'mz_R',
                         'mz_left',
                         'mz_right',
                         'mz',
                         'bottom',
                         'top',
                         'color',
                         'intensity',
                         'molname',
                         'molsource',
                         'molformula',
                         'molq',
                         'molg',
                         'molintensity'),
                        lists))

        source = ColumnDataSource(data)

        fat_rectangles = plot.quad(top='top', bottom='bottom',
                                   left='mz_left', right='mz_right',
                                   color='color', source=source,
                                   alpha=.5)

        slim_rectangles = plot.quad(top='top', bottom='bottom',
                                   left='mz_L', right='mz_R',
                                   color='color', source=source)

        hover_fat = HoverTool(renderers=[fat_rectangles],
                              tooltips=[('peak intensity', "@intensity{0,0}"),
                                        ('m/z interval', mz_representation),
                                        ('name', '@molname'),
                                        ('formula','@molformula'),
                                        ('q', '@molq'),
                                        ('g', '@molg'),
                                        ('total intensity', "@molintensity{0,0}")])

        plot.add_tools(hover_fat)

        # plotting cluster / local quality of peak fitting
        cluster_lists = zip(*((c.mz_left,
                               c.mz_right,
                               (c.mz_left + c.mz_right)/2.0,
                               c.intensity,
                               c.estimate,
                               abs(c.intensity-c.estimate))
                              for c in self.clusters))

        data_clusters = dict(zip(('mz_left',
                                  'mz_right',
                                  'mz',
                                  'intensity',
                                  'estimate',
                                  'abserror'),
                                 cluster_lists))

        source_clusters = ColumnDataSource(data_clusters)

        cluster_intensities = plot.segment(x0='mz_left',
                                           x1='mz_right',
                                           y0='intensity',
                                           y1='intensity',
                                           color='red',
                                           line_width=3, source=source_clusters)

        hover_clusters = HoverTool(renderers=[cluster_intensities],
                                   tooltips=[('total intensity', "@intensity{0,0}"),
                                             ('total estimate', "@estimate{0,0}"),
                                             ('absolute error', "@abserror{0,0}"),
                                             ('m/z', mz_representation)])

        plot.add_tools(hover_clusters)

        # Experimental Squares
        raw_spectrum = plot.square(x=self.masstodon.spectrum.mz,
                                   y=self.masstodon.spectrum.intensity,
                                   size=5,
                                   color='black',
                                   alpha=.5)

        x_representation = make_string_represenation('x',
                                                     self.masstodon.mz_digits)

        hover_squares = HoverTool(renderers=[raw_spectrum],
                                  tooltips=[('intensity', "@y{0,0}"),
                                            ('m/z', x_representation)])
        plot.add_tools(hover_squares)

        show(plot)
        return plot
