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
import  networkx        as      nx
from    math            import  sqrt
import  numpy           as      np
from    collections     import Counter

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
            if total_prob < minimal_prob:  # gonna remove some stupid nodes.
                for I in cc[M]:
                    nodes_to_remove.append(I)
                nodes_to_remove.append(M)
    cc.remove_nodes_from(nodes_to_remove)


def add_zero_intensity_G_nodes(P):
    '''Pair unpaired isotopologue peaks I with zero-intensity experimental groups G.
    '''
    cnts = Counter(P.node[N]['type'] for N in P)
    Gcnt = cnts['G']+1
    newGnodes = []
    newGIedges= []
    for I in P:
        if P.node[I]['type'] == 'I':
            if len(P[I]) == 1:
                G = 'G' + str(Gcnt)
                newGnodes.append( (G,{'intensity':0.0, 'type':'G'}) )
                newGIedges.append( (G,I) )
                Gcnt += 1
    P.add_nodes_from(newGnodes)
    P.add_edges_from(newGIedges)


def create_G_nodes(SFG):
    '''Collect experimental peaks into groups of experimental data G.

    This is done without any loss in resolution.'''
    E2remove = []
    Gs = Counter()
    for E in SFG:
        if SFG.node[E]['type'] == 'E':
            G = frozenset(SFG[E])
            Gs[G] += SFG.node[E]['intensity']
            E2remove.append(E)
    SFG.remove_nodes_from(E2remove)
    for Gcnt, (Is, G_intensity) in enumerate(Gs.items()):
        G = 'G' + str(Gcnt)
        SFG.add_node(G, intensity=G_intensity, type='G')
        for I in Is:
            SFG.add_edge(I, G)
    add_zero_intensity_G_nodes(SFG)


class PeakPicker(object):
    '''Class for peak picking.'''
    def __init__(   self,
                    Forms,
                    IsoCalc,
                    mzPrec = 0.05 ):
        self.Forms  = Forms
        self.IsoCalc= IsoCalc
        self.mzPrec = mzPrec


    def represent_as_BFG(self, massSpectrum):
        '''Prepare the Big Graph based on mass spectrum and the formulas.'''
        ePeaks  = Itree( II( mz-self.mzPrec, mz+self.mzPrec, (mz, intensity) ) for mz, intensity in zip(*massSpectrum))
        iso_cnt = 0
        BFG     = nx.Graph()
        for cnt, (molType, formula, aaNo, q, g) in enumerate(self.Forms.makeMolecules()):
            iso_mzs, iso_intensities = self.IsoCalc.isoEnvelope( atomCnt_str=formula, q=q, g=g )
            M = 'M'+str(cnt)
            BFG.add_node(M, formula=formula, type='M', q=q, g=g, molType=molType )
            for isoMZ, isoI in zip(iso_mzs, iso_intensities):
                I = 'I' + str(iso_cnt)
                BFG.add_node(I, mz=isoMZ, intensity=isoI, type='I')
                BFG.add_edge(M, I)
                for exp in ePeaks[isoMZ]:
                    expMZ, expI = exp.data
                    if not expMZ in BFG:
                        BFG.add_node(expMZ, intensity=expI, type='E')
                    BFG.add_edge(I, expMZ)
                iso_cnt += 1
        return BFG


    def get_problems(self, massSpectrum, minimal_prob_per_molecule=.7):
        '''Enumerate deconvolution problems.'''
        BFG = self.represent_as_BFG(massSpectrum)
        for cc in nx.connected_component_subgraphs(BFG):
            # considers only good connected components
            if contains_experimental_peaks(cc):
                trim_unlikely_molecules(cc, minimal_prob_per_molecule)
                for SFG in nx.connected_component_subgraphs(cc):
                    if len(SFG) > 1:
                        create_G_nodes(SFG)
                        yield SFG
