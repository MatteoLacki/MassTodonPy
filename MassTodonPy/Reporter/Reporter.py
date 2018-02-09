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
from collections import Counter
import csv
import json
import os

from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Reporter.buffers import buffers
from MassTodonPy.Reporter.misc import Brick, Cluster, PeakGroup, ColorGenerator, make_string_represenation
from MassTodonPy.Write.csv_tsv import write_rows


def float2str(x):
    if isinstance(x, float):
        return "{:10.3f}".format(x)
    else:
        return x

def float2strPerc(x):
    if isinstance(x, float):
        return "{:10.3f}%".format(x * 100)
    else:
        return x


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
        self._solutions = solutions
        self._molecules = molecules
        self._bricks = []
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

    def iter_precise_estimates(self):
        """Iterate over precise estimates.

        Supply a flow of rows to write to a csv file with output."""

        yield ('molecule',
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
               'right peak group m/z')
        for b in self._bricks:
            yield (b.molecule.name,
                   b.molecule.source.name,
                   b.molecule.formula.str_with_charges(b.molecule.q, b.molecule.g),
                   b.molecule.q,
                   b.molecule.g,
                   round(b.intensity),
                   round(b.molecule.intensity),
                   b.intensity/b.molecule.intensity,
                   b.molecule.monoisotopic_mz,
                   (b.peak_group.mz_L + b.peak_group.mz_R)/2.0,
                   b.peak_group.mz_L,
                   b.peak_group.mz_R)

    def iter_local_fit_quality(self):
        """Iterate over quality of fit to the experimental groups G."""
        yield ('mz',
               'left peak group m/z',
               'right peak group m/z',
               'intensity',
               'estimate',
               'absolute error')
        for c in self._peak_groups:
            if int(c.intensity) > 0 and int(c.estimate) > 0:
                yield ((c.mz_L + c.mz_R)/2.0,
                       c.mz_L,
                       c.mz_R,
                       int(c.intensity),
                       int(c.estimate),
                       abs(int(c.intensity) - int(c.estimate)))

    def iter_molecule_estimates(self, include_zero_intensities=True):
        """Iterate over the estimates of the intensity of molecules."""
        self._molecules.sort(key=lambda m: m.intensity, reverse=True)
        yield ('name',
               'formula',
               'charge',
               'quenched charge',
               'estimated intensity',
               'source',
               'source fasta',
               'source formula')
        for m in self._molecules:
            if round(m.intensity) > 0 or include_zero_intensities:
                yield (m.name,
                       m.formula.str_with_charges(m.q, m.g),
                       m.q,
                       m.g,
                       round(m.intensity),
                       m.source.name,
                       m.source.fasta,
                       m.source.formula.str_with_charges(m.source.q))

    def iter_global_quality_fits(self):
        """Get global quality fits."""
        for sol in self._solutions:
            yield sol.global_fit_quality()

    def global_quality_fits_stats(self):
        """Iterate over statistics of global quality of fits."""
        T = Counter()
        for fq in self.iter_global_quality_fits():
            status = fq.pop('status')
            T += fq

        within_tolerance = {}
        within_tolerance['l1'] = T['l1'] / T['total intensity']
        within_tolerance['overestimates'] = T['overestimates'] / T['total intensity']
        within_tolerance['underestimates'] = T['underestimates'] / T['total intensity']

        thresholding = {}
        high = self._spectrum.l1()
        thresholding['l1'] = (high + T['l1'] - T['total intensity']) / high
        thresholding['overestimates'] = T['overestimates'] / high
        thresholding['underestimates'] = (high - T['total intensity'] + T['underestimates']) / high

        whole_spectrum = {}
        total_intensity = self._spectrum.l1() + self._spectrum.low_spectrum.l1()
        unassigned_intensity = total_intensity - T['total intensity']
        whole_spectrum['l1'] = (T['l1'] + unassigned_intensity) / total_intensity
        whole_spectrum['overestimates'] = T['overestimates'] / total_intensity
        whole_spectrum['underestimates'] = (T['underestimates'] + unassigned_intensity) / total_intensity

        return {'within tolerance': within_tolerance,
                'after thresholding': thresholding,
                'of the whole spectrum': whole_spectrum}


    def iter_global_quality_fits_rows(self):
        """Iterate over global fits quality information."""
        col_names = ('l1', 'l2', 'overestimates','underestimates', 'total intensity', 'status')
        yield col_names
        T = Counter()
        for fq in self.iter_global_quality_fits():
            status = fq.pop('status')
            T += fq
            fq['status'] = status
            yield tuple(float2str(fq[x]) for x in col_names)
        yield ('Total','Total','Total','Total','Total')
        yield tuple(float2str(T[x]) for x in col_names[0:-1])
        yield (' ',)
        yield ('Errors',)
        yield ('l1', 'overestimates','underestimates')

        stats = self.global_quality_fits_stats()
        stats = [(t, stats[t]) for t in ('within tolerance',
                                         'after thresholding',
                                         'of the whole spectrum')]
        for text, d in stats:
            yield (float2strPerc(d['l1']),
                   float2strPerc(d['overestimates']),
                   float2strPerc(d['underestimates']),
                   'relative to intensity',
                   text)

    def write(self, path, include_zero_intensities=True):
        """Write results to a file.

        Parameters
        ==========
        path : str
            Path to the output you want MassTodon to save its results.

        """
        write_rows(self.iter_precise_estimates(), path + 'assigned_spectrum_precise.csv')
        write_rows(self.iter_local_fit_quality(), path + 'assigned_spectrum_local_fits.csv')
        write_rows(self.iter_molecule_estimates(include_zero_intensities),
                                                path + 'estimates_of_molecule_intensities.csv')
        write_rows(self.iter_global_quality_fits_rows(), path + 'global_fit_quality.csv')

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
