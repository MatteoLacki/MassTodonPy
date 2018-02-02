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
import csv
import json
import os

from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Reporter.buffers import buffers
from MassTodonPy.Reporter.misc import Brick, Cluster, PeakGroup, ColorGenerator, make_string_represenation


class Reporter(object):
    """Report MassTodon results."""

    def __init__(self,
                 solutions,
                 molecules,
                 mz_digits,
                 spectrum,
                 max_buffer_len=0.5,
                 **kwds):
        """Reporting results of the MassTodon.

        Parameters
        ==========
        solutions :
            An iterable of solutions to the faced deconvolution problem.
        molecules :
            An iterable of generated molecules.
        mz_digits : int
            The number of significant digits used for m/z representation.
        max_buffer_len : float
            The maximal length of the visual buffer between peaks, i.e.
            the big rectangle width.
        kwds :
            Arguments to other methods.
        """
        self._bricks = []
        self._solutions = solutions
        self._molecules = molecules
        self._peak_groups = []
        self._mz_digits = int(mz_digits)
        self._min_interval_len = 10**(-self._mz_digits)
        self._spectrum = spectrum
        self._max_buffer_len = float(max_buffer_len)
        for name, thing in self._peakGroups_bricks_clusters():
            self.__dict__[name].append(thing)

        # order peak groups by m/z ratios
        self._peak_groups.sort(key=lambda c: c.mz_L)

        # adding left and right buffers needed for plot
        mz_L, mz_R = list(zip(*((c.mz_L, c.mz_R) for c in self._peak_groups))) #TODO this is an overkill: smooth out the defintion of the buffers function.
        mz_lefts, mz_rights = buffers(mz_L, mz_R,
                                      max_length=self._max_buffer_len)
        for c, mz_left, mz_right in zip(self._peak_groups, mz_lefts, mz_rights):
            c.mz_left = mz_left
            c.mz_right = mz_right

        # lexicographically order the bricks
        self._bricks.sort(key=lambda b: ( b.peak_group.mz_L,
                                         -b.molecule.intensity) )

        # adjust brick top/bottom/color
        self.color = ColorGenerator()
        prev_b = Brick(top=0.0, bottom=0.0, peak_group=PeakGroup(mz_L=0.0)) # guardian
        for b in self._bricks:
            if b.peak_group.mz_L is prev_b.peak_group.mz_L:
                b.bottom = prev_b.top
                b.top = prev_b.top + b.intensity
            b.color = self.color(b.molecule)
            prev_b = b

        # build up clusters: as many as solutions.
        self._clusters = [Cluster() for _ in
                          range(len(self._solutions))]
        for b in self._bricks:
            self._clusters[b.peak_group.sol_id].update(b)

        # make this once
        self.assigned_spectrum_data = self.get_assigned_spectrum_data() # plot data

    def _peakGroups_bricks_clusters(self):
        """Generate a flow of peak groups, bricks, and clusters."""
        for sol_id, sol in enumerate(self._solutions):
            for G in sol:
                if G[0] is 'G':
                    d = self._min_interval_len / 2.0
                    try:
                        min_mz = sol.node[G]['min_mz'] - d
                        max_mz = sol.node[G]['max_mz'] + d
                    except KeyError:
                        mz = sol.node[G]['mz']
                        min_mz = mz - d
                        max_mz = mz + d
                    peak_group = PeakGroup(mz_L=min_mz,
                                           mz_R=max_mz,
                                           mz_left=min_mz,     # updated later
                                           mz_right=max_mz,    # updated later
                                           intensity=sol.node[G]['intensity'],
                                           intensity_d=sol.node[G]['intensity']/
                                                       (max_mz - min_mz),
                                           estimate=sol.node[G]['estimate'],
                                           estimate_d=sol.node[G]['estimate']/
                                                      (max_mz - min_mz),
                                           sol_id=sol_id)
                    yield ('_peak_groups', peak_group)
                    for I in sol[G]:
                        estimate = sol[G][I]['estimate']
                        for M in sol[I]:
                            if M[0] is 'M':
                                # CVXOPT sometimes returns results only
                                # slightly smaller than zero: don't want them
                                #                 |||||
                                if round(estimate) > 0:
                                    yield ('_bricks',
                                            Brick(peak_group=peak_group,
                                                  top=estimate,
                                                  bottom=0.0,
                                                  color='#440154',# updated later
                                                  intensity=estimate,
                                                  molecule=sol.node[M]['molecule']))


    def aggregate(self):
        """Aggregate the results of MassTodon."""
        pass

    def write(self, path, include_zero_intensities=True):
        """Write results to a file.

        Parameters
        ==========
        file_path : str
            Path to the output you want MassTodon to save its results.
            Supported path extensions are .txt and .csv.
        """
        file_path, file_name, file_ext = parse_path(path)

        assert file_ext in ('.csv', '.tsv'), "Writing only to csv or tsv."
        delimiter = ',' if file_ext == '.csv' else '\t'

        # detailed intormation on assignment
        path_details = "{0}{1}_precise{2}".format(file_path, file_name, file_ext)
        with open(path_details, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(('molecule',
                             'source',
                             'formula',
                             'charge',
                             'quenched charge',
                             'intensity',
                             'molecule intensity',
                             'isotopologues probability',
                             'monoisotopic m/z',
                             'average m/z',
                             'left peak group m/z',
                             'right peak group m/z'))
            for b in self._bricks:
                writer.writerow((b.molecule.name,
                                 b.molecule.source.name,
                                 # str(b.molecule.formula),
                                 b.molecule.formula.str_with_charges(b.molecule.q,
                                                                     b.molecule.g),
                                 b.molecule.q,
                                 b.molecule.g,
                                 round(b.intensity),
                                 round(b.molecule.intensity),
                                 b.intensity/b.molecule.intensity,
                                 b.molecule.monoisotopic_mz,
                                 (b.peak_group.mz_L + b.peak_group.mz_R)/2.0,
                                 b.peak_group.mz_L,
                                 b.peak_group.mz_R))

        # precise quality of fit information: on G level basis
        path_local_fits = "{0}{1}_local_fits{2}".format(file_path, file_name, file_ext)
        with open(path_local_fits, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(('mz',
                             'left peak group m/z',
                             'right peak group m/z',
                             'intensity',
                             'estimate',
                             'absolute error'))
            for c in self._peak_groups:
                writer.writerow(((c.mz_L + c.mz_R)/2.0,
                                 c.mz_L,
                                 c.mz_R,
                                 int(c.intensity),
                                 int(c.estimate),
                                 abs(int(c.intensity) - int(c.estimate))))

        # aggregated intormation on assignment
        path_molecules = "{0}{1}_molecules{2}".format(file_path,
                                                      file_name,
                                                      file_ext)

        self._molecules.sort(key=lambda m: m.intensity, reverse=True)

        with open(path_molecules, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(('name',
                             'formula',
                             'charge',
                             'quenched charge',
                             'estimated intensity',
                             'source',
                             'source fasta',
                             'source formula'))

            for m in self._molecules:
                if round(m.intensity) > 0 or include_zero_intensities:
                    writer.writerow((m.name,
                                     m.formula.str_with_charges(m.q, m.g),
                                     m.q,
                                     m.g,
                                     round(m.intensity),
                                     m.source.name,
                                     m.source.fasta,
                                     m.source.formula.str_with_charges(m.source.q)))

    def get_assigned_spectrum_data(self):
        """Make data for the plot with assigned spectrum."""

        mz_repr = make_string_represenation('mz', self._mz_digits)
        x_repr = make_string_represenation('x', self._mz_digits)

        out = {'y_range_start':  0.0,
               'x_label':       'mass/charge',
               'y_label':       'intensity',
               'mz_digits':      self._mz_digits}

        out['tools'] = "crosshair pan wheel_zoom box_zoom undo redo reset box_select save".split(" ")

        # vertical experimental bars
        out['exp_vbar'] = {'x':     list(self._spectrum.mz),
                           'top':   list(self._spectrum.intensity),
                           'width': self._min_interval_len,
                           'color': 'black',
                           'alpha': 0.1}

        # Horizontal threshold line
        out['threshold_line'] = {'intensity': self._spectrum.min_intensity,
                                 'args': [(max(min(self._spectrum.mz) - 50, 0),
                                           max(self._spectrum.mz) + 50 ),
                                          (self._spectrum.min_intensity,
                                           self._spectrum.min_intensity)],
                                 'kwds': {'line_width': 2, 'color': 'red'}}

        # bricks: divisions of estimated peaks into constituents
        lists = zip(*((b.peak_group.mz_L,
                       b.peak_group.mz_R,
                       b.peak_group.mz_left,
                       b.peak_group.mz_right,
                       (b.peak_group.mz_L + b.peak_group.mz_R)/2.0,
                       b.bottom,
                       b.top,
                       b.color,
                       b.intensity,
                       b.molecule.name,
                       b.molecule.source.name,
                       # str(b.molecule.formula),
                       b.molecule.formula.str_with_charges(b.molecule.q,
                                                           b.molecule.g),
                       b.molecule.q,
                       b.molecule.g,
                       b.molecule.intensity)
                      for b in self._bricks))

        out['rectangle_data'] = dict(zip(('mz_L',
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

        out['fat_rectangles'] = {'top':     'top',
                                 'bottom':  'bottom',
                                 'left':    'mz_left',
                                 'right':   'mz_right',
                                 'color':   'color',
                                 'alpha':    0.5}

        out['fat_rectangles_tooltips'] = [('peak intensity', "@intensity{0,0}"),
                                          ('m/z interval', mz_repr),
                                          ('name', '@molname'),
                                          ('formula','@molformula'),
                                          ('q', '@molq'),
                                          ('g', '@molg'),
                                          ('total intensity', "@molintensity{0,0}")]

        out['slim_rectangles'] = {'top': 'top',
                                  'bottom': 'bottom',
                                  'left': 'mz_L',
                                  'right': 'mz_R',
                                  'color': 'color'}
        # peak groups data: for local quality of fit
        peak_group_lists = zip(*((pg.mz_left,
                                  pg.mz_right,
                                  (pg.mz_left + pg.mz_right)/2.0,
                                  pg.intensity,
                                  pg.estimate,
                               abs(pg.intensity-pg.estimate))
                              for pg in self._peak_groups))

        out['peak_groups_data'] = dict(zip(('mz_left',
                                            'mz_right',
                                            'mz',
                                            'intensity',
                                            'estimate',
                                            'abserror'),
                                       peak_group_lists))

        out['peak_groups'] = {'x0': 'mz_left',
                              'x1': 'mz_right',
                              'y0': 'intensity',
                              'y1': 'intensity',
                              'color': 'red',
                              'line_width': 3}

        out['peak_groups_tooltips'] = [('total intensity', "@intensity{0,0}"),
                                       ('total estimate', "@estimate{0,0}"),
                                       ('absolute error', "@abserror{0,0}"),
                                       ('m/z', mz_repr)]
        # Experimental Squares
        out['experimental_squares'] = {'x': list(self._spectrum.mz),
                                       'y': list(self._spectrum.intensity),
                                       'size': 5,
                                       'color':
                                       'black',
                                       'alpha': 0.5}
        out['experimental_squares_tooltips'] = [('intensity', "@y{0,0}"),
                                                ('m/z', x_repr)]
        # Text Labels
        cluster_lists = zip(*((c.mz,
                               c.intensity,
                               c.mol_intensity,
                               c.mol_name)
                              for c in self._clusters))

        out['cluster_data'] = dict(zip(('mz',
                                        'intensity',
                                        'mol_intensity',
                                        'mol_name'),
                                        cluster_lists))

        out['labels'] = {'x': 'mz',
                         'y': 'intensity',
                         'text': 'mol_name',
                         'level': 'glyph',
                         'x_offset': 0,
                         'y_offset': 1,
                         'render_mode': 'css'}
        return out

    def to_json(self, path):
        """Export plot data to json.

        Parameters
        ==========
        path : str
            Where to save the json.

        """
        out = {'assigned_spectrum': self.assigned_spectrum_data}
        output_path, file_name, _ = parse_path(path)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        with open(path, 'w') as f:
            json.dump(out, f)

    def plot(self,
             path="assigned_spectrum.html",
             mode="inline",
             show_plot=True,
             width=None,
             height=None,
             _mult=1,
             **kwds):
        """Make a plot.

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
        PD = self.assigned_spectrum_data
        output_file(path, mode=mode)
        if not width:
            width = 800 * _mult
        if not height:
            height = 600 * _mult
        plot = figure(plot_width=width,
                      plot_height=height,
                      tools=PD['tools'])
        plot.y_range.start = PD['y_range_start']
        plot.xaxis.axis_label = PD['x_label']
        plot.yaxis.axis_label = PD['y_label']

        # The experimental data bars
        experimental_bars = plot.vbar(**PD['exp_vbar'])

        # Horizontal threshold line
        if PD['threshold_line']['intensity'] > 1.0:
            plot.line(*PD['threshold_line']['args'],
                      **PD['threshold_line']['kwds'])

        # Plotting rectangles
        source_rectangles = ColumnDataSource(PD['rectangle_data'])
        fat_rectangles = plot.quad(source=source_rectangles,
                                   **PD['fat_rectangles'])

        slim_rectangles = plot.quad(source=source_rectangles,
                                    **PD['slim_rectangles'])

        hover_fat = HoverTool(renderers=[fat_rectangles],
                              tooltips=PD['fat_rectangles_tooltips'])
        plot.add_tools(hover_fat)

        # plotting peak_group / local quality of peak fitting
        source_peak_groups = ColumnDataSource(PD['peak_groups_data'])
        peak_group_intensities = plot.segment(source=source_peak_groups,
                                              **PD['peak_groups'])

        hover_peak_groups = HoverTool(renderers=[peak_group_intensities],
                                      tooltips=PD['peak_groups_tooltips'])
        plot.add_tools(hover_peak_groups)

        # Experimental Squares
        raw_spectrum = plot.square(**PD['experimental_squares'])
        hover_squares = HoverTool(renderers=[raw_spectrum],
                                  tooltips=PD['experimental_squares_tooltips'])
        plot.add_tools(hover_squares)

        # Text labels
        source_clusters = ColumnDataSource(PD['cluster_data'])
        labels = LabelSet(source=source_clusters,
                          **PD['labels'])
        plot.add_layout(labels)
        if show_plot:
            show(plot)
        return plot
