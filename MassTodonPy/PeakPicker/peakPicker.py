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

from    intervaltree    import Interval as II, IntervalTree as Itree
import  networkx        as      nx
from    math            import  sqrt
import  numpy           as      np
from    collections     import Counter

def contains_experimental_peaks(cc):
    return any(isinstance(N, float) for N in cc)

def trim_unlikely_molecules(cc, minimal_prob=0.7):
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

def getGraphs(ccs, minimal_prob=.7):
    for cc in ccs:
        # consider only good connected components
        if contains_experimental_peaks(cc):
            trim_unlikely_molecules(cc)
            for G in nx.connected_component_subgraphs(cc):
                if len(G) > 1:
                    yield G

def group_experimental_peaks(P):
    E2remove = []
    Gs = Counter()
    for E in P:
        if P.node[E]['type'] == 'E':
            G = frozenset(P[E])
            Gs[G] += P.node[E]['intensity']
            E2remove.append(E)
    P.remove_nodes_from(E2remove)
    for Gcnt, (Is, G_intensity) in enumerate(Gs.items()):
        G = 'G' + str(Gcnt)
        P.add_node(G, intensity=G_intensity, type='G')
        for I in Is:
            P.add_edge(I, G)


class PeakPicker():
    '''Class for peak picking.'''

    def __init__(self,
            formulator,
            isotopeCalculator,
            chebyshevCoverage       = 0.99,
            jointProbabilityIsoSpec = 0.999,
            precisionDigits         = 2,
            precisionMass           = 0.05   ):

        self.cheb   = 1.0 - chebyshevCoverage
        self.forms  = formulator
        self.isoCalc= isotopeCalculator
        self.jP     = jointProbabilityIsoSpec
        self.prec   = precisionDigits
        self.mzPrec = precisionMass

    def BFG_representation(self, massSpectrum):
        ePeaks  = Itree( II( mz-self.mzPrec, mz+self.mzPrec, (mz, intensity) )
                            for mz, intensity in np.rollaxis(massSpectrum,0))
        iso_cnt = 0
        BFG     = nx.Graph()
        for cnt, (mType, formula, aaNo, q, g) in enumerate(self.forms.makeMolecules()):
            iso_mzs, iso_intensities = self.isoCalc.isoEnvelope(formula, self.jP, q, g)
            M = mType + '_' + str(q) + '_' + str(g)
            BFG.add_node(M, formula=formula, type='M')
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
