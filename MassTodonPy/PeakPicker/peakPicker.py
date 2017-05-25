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

from    intervaltree    import Interval as II, IntervalTree as Itree
import  networkx        as     nx
from    math            import sqrt
import  numpy           as     np
from    collections     import Counter, defaultdict


inf = float('inf')

class MultiCounter(Counter):
    '''Callable Counter for special field operations.

    Offers stroke to any functional purist.
    '''
    def __call__(self, key):
        val = self[key]
        self[key] += 1
        return str(key) + str(val)


def contains_experimental_peaks(cc):
    '''Check if a given problem contains any experimental peaks.'''
    return any(isinstance(N, float) for N in cc)

def trim_unlikely_molecules(cc, minimal_prob=0.7):
    '''Trim molecules whose isotopic envelopes cover potentially less than the minimal_prob threshold.'''
    nodes_to_remove = []
    for M in cc:
        if cc.node[M]['type'] == 'M':  # we are looking at a molecule node
            total_prob = 0.0
            for I in cc[M]:
                if len(cc[I]) > 1:
                    total_prob += cc.node[I]['intensity']
            if total_prob < minimal_prob:  # g\onna remove some stupid nodes.
                for I in cc[M]:
                    nodes_to_remove.append(I)
                nodes_to_remove.append(M)
    cc.remove_nodes_from(nodes_to_remove)

class PeakPicker(object):
    '''Class for peak picking.'''
    def __init__(   self,
                    Forms,
                    IsoCalc,
                    mzPrec = 0.05 ):

        self.Forms  = Forms
        self.IsoCalc= IsoCalc
        self.mzPrec = mzPrec
        self.cnts   = MultiCounter() # TODO finish it.

    def represent_as_Graph(self, massSpectrum):
        '''Prepare the Graph based on mass spectrum and the formulas.'''
        
        Exps = Itree( II( mz-self.mzPrec, mz+self.mzPrec, (mz, intensity) ) for mz, intensity in zip(*massSpectrum) )

        Graph = nx.Graph()
        for M_type, M_formula, _, M_q, M_g in self.Forms.makeMolecules():
            I_mzs, I_intensities = self.IsoCalc.isoEnvelope(atomCnt_str=M_formula, q=M_q, g=M_g)
            M = self.cnts('M')
            Graph.add_node(M,
                           formula = M_formula,
                           type = 'M',
                           q = M_q,
                           g = M_g,
                           molType = M_type )
            for I_mz, I_intensity in zip(I_mzs, I_intensities):
                I = self.cnts('I')
                Graph.add_node(I,
                               mz = I_mz,
                               intensity = I_intensity,
                               type = 'I' )
                Graph.add_edge(M, I)
                for E_interval in Exps[I_mz]:
                    # experimental peaks are stored by their m/z
                    # they must have been aggregated and so each m/z corresponds to one peak
                    E_mz, E_intensity = E_interval.data
                    if not E_mz in Graph:
                        Graph.add_node(E_mz,
                                       intensity = E_intensity, type = 'E' )
                    Graph.add_edge(I, E_mz)
        return Graph


    def get_problems(self, massSpectrum, minimal_prob_per_molecule=.7):
        '''Enumerate deconvolution problems.'''
        Graph = self.represent_as_Graph(massSpectrum)
        for cc in nx.connected_component_subgraphs(Graph):
            if contains_experimental_peaks(cc):
                trim_unlikely_molecules(cc, minimal_prob_per_molecule)
                for small_graph in nx.connected_component_subgraphs(cc):
                    if len(small_graph) > 1:
                        self.add_G_nodes(small_graph)
                        yield small_graph


    def add_G_nodes(self, small_graph):
        '''Collect experimental peaks into groups of experimental data G.

        This is done without any loss in resolution, but surely in a stupid way.'''
        iso_vals = Itree()    # m/z intervals around isotopologues
        E_to_remove = []
        Gs_intensity = Counter()
        Gs_min_mz = defaultdict(lambda:inf)
        Gs_max_mz = Counter()

        for E in small_graph:
            if small_graph.node[E]['type'] == 'E':
                I_of_G = frozenset(small_graph[E])
                Gs_intensity[I_of_G] += small_graph.node[E]['intensity']
                Gs_min_mz[I_of_G] = min(Gs_min_mz[I_of_G], E)
                Gs_max_mz[I_of_G] = max(Gs_max_mz[I_of_G], E)
                E_to_remove.append(E)
                for I in I_of_G:
                    I_mz = small_graph.node[I]['mz']
                    L, R = I_mz-self.mzPrec, I_mz+self.mzPrec
                    iso_vals.addi(L, R)
        small_graph.remove_nodes_from(E_to_remove)

        iso_vals.split_overlaps()
        for Gcnt, I_of_G in enumerate(Gs_intensity):
            min_mz, max_mz = Gs_min_mz[I_of_G], Gs_max_mz[I_of_G]
            if min_mz == max_mz:
                mz = iso_vals[min_mz]
            else:
                mz = iso_vals[II(min_mz,max_mz)]
            assert len(mz) == 1, "There should be only one interval containing the interval [min_mz,max_mz] among the 'tesselation' of the m/z axis."
            mz = mz.pop()
            G = self.cnts('G')
            small_graph.add_node(G, intensity=Gs_intensity[I_of_G], type='G', mz=mz)
            for I in I_of_G:
                small_graph.add_edge(I, G)

        # zero intensity G nodes
        newGnodes = []
        newGIedges= []
        for I in small_graph:
            if small_graph.node[I]['type'] == 'I':
                if len(small_graph[I]) == 1:
                    I_mz = small_graph.node[I]['mz']
                    mz = II(I_mz, I_mz)
                    G = self.cnts('G')
                    newGnodes.append( (G,{'intensity':0.0, 'type':'G', 'mz':mz }) )
                    newGIedges.append( (G,I) )
        small_graph.add_nodes_from(newGnodes)
        small_graph.add_edges_from(newGIedges)
