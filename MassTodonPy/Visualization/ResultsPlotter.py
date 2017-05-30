#
# -*- coding: utf-8 -*-
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
#

from intervaltree import Interval as interval, IntervalTree


class ResultsPlotter(object):
    def __init__( self, mz_prec ):
        self.mz_prec = mz_prec
        self.BG = None

    def iBG(self, node_type):
        '''Iterate over all nodes of a given type in the whole big graph **BG** of solved subproblems.

        node_type - either
        '''
        assert node_type in ('G','I','M'), "specified wrong type of node. Was %s. Should be G, I, M" % node_type

        for r in self.BG:
            SG = r['SG']
            for N in SG:
                N_D = SG.node[N]
                if N_D['type'] == node_type:
                    yield N, N_D

    def add_mz_ranges_to_results(self, masstodon_res):
        '''Add information about the m/z ranges for the G nodes.'''

        self.BG = masstodon_res
        prec = self.mz_prec
        I_tree = IntervalTree(interval(I_D['mz']-prec,I_D['mz']+prec) for _,I_D in self.iBG('I'))
        I_tree.split_overlaps()
        for G, G_D in self.iBG('G'):
            min_mz, max_mz = G_D['min_mz'], G_D['max_mz']
            if min_mz == max_mz:
                ## This copes with stupid border issues.
                ## iso_intervals = I_tree[II(min_mz-prec/10., max_mz+prec/10.)]
                # set digits to mz_prec/10
                iso_intervals = I_tree[min_mz]
            else:
                iso_intervals = I_tree[interval(min_mz, max_mz)]
            if len(iso_intervals) == 1:
                mz = iso_intervals.pop()
                G_D['mz_L'] = mz.begin
                G_D['mz_R'] = mz.end
            else:
                G_D['mz_L'] = G_D['mz_R'] = None


    def G_info_iter(self, full_info=False):
        '''Iterate over all information on experimental groupings G.'''
        for r in self.BG:
            SG = r['SG']
            for G in SG:
                if SG.node[G]['type'] == 'G':
                    G_D = SG.node[G]
                    if not full_info and G_D['mz_L'] and G_D['mz_R']:
                        yield { 'mz_L': G_D['mz_L'],
                                'mz_R': G_D['mz_R'],
                                'tot_estimate': G_D['estimate'],
                                'tot_intensity':G_D['intensity'] }
                    else:
                        for I in SG[G]:
                            for M in SG[I]:
                                if SG.node[M]['type'] == 'M':
                                    M_D = SG.node[M]
                                    yield { 'formula':  M_D['formula'],
                                            'molType':  M_D['molType'],
                                            'q':        M_D['q'],
                                            'g':        M_D['g'],
                                            'mz_L':     G_D['mz_L'],
                                            'mz_R':     G_D['mz_R'],
                                            'estimate': SG.edge[G][I]['estimate'],
                                            'tot_estimate': G_D['estimate'],
                                            'tot_intensity':G_D['intensity']    }
