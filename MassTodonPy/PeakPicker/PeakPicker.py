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

from collections import Counter
from collections import defaultdict
from intervaltree import Interval as II
from intervaltree import IntervalTree
import networkx as nx
from six.moves import zip
from time import time

inf = float('inf')


def trim_unlikely_molecules(cc, minimal_prob=0.7):
    """Trim theoretic molecules that are unlikely to be in real spectrum.

    If the joint probability of their theoretical peaks close to some
    experimental peaks is less than minimal_prob,
    then do not assume that they are in the spectrum.
    This is highly conservative for such an evil leftist like me.
    """
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
    return cc


class PeakPicker(object):
    """Class for peak picking."""

    def __init__(self,
                 molecules,
                 isotopic_calculator,
                 mz_precision_digits=0.05,
                 _verbose=False):
        """Initialize peak picker."""
        self.molecules = molecules
        self.isotopic_calculator = isotopic_calculator
        self.mz_precision_digits = mz_precision_digits
        self.M = 0  # the number of molecule nodes
        self.I = 0  # the number of isotopologue nodes
        self._verbose = _verbose
        self.stats = Counter()
        self.Used_Exps = set()
        self.G_stats = []  # logger

    def represent_as_Graph(self, spectrum):
        """Prepare the Graph based on mass spectrum and the formulas.

        Parameters
        ----------
        spectrum : tuple
            The experimental spectrum, a tuple consisting of
            a m/z numpy array and an intensities numpy array.
        Returns
        -------
        Graph : graph
            The deconvolution graph.

        """
        prec = self.mz_precision_digits
        Exps = IntervalTree(II(mz - prec, mz + prec, (mz, intensity))
                            for mz, intensity in spectrum)
        # TODO eliminate the above using binary search solo

        Graph = nx.Graph()
        no_exp_per_iso = Counter()

        if self._verbose:  # logger...
            print('Representing problem as a graph.')

        T0 = time()

        for mol in self.molecules:
            I_mzs, I_intensities =\
                self.isotopic_calculator.get_envelope(formula=mol.formula,
                                                      q=mol.q,
                                                      g=mol.g,
                                                      memoize=True)

            # Why do we make a node before establishing if it should be there?
            # Because it's nice to know, that the formula was somewhere there?
            # But it's among the formulas list anyway...
            # But it's more elegant not to store the list, but make a generator
            # that will populate the networkx graph.
            # Follow the logic of a memory efficient problem representation.

            # Maybe subclass networkx graph to get additional functions?
            # Like DeconvolutionGraph.add_formula(Formula).
            Graph.add_node(self.M,
                           formula=M_formula,
                           type='M',
                           q=M_q,
                           g=M_g,
                           molType=M_type)

            self.I += len(I_mzs)
            self.M += 1

            for I_mz, I_intensity in zip(I_mzs, I_intensities):
                I = self.cnts('I')

                # Should the DeconvolutionGraph contain all isotopologues?
                # Like why? If they are unnecessary then better simply to keep
                # the molecule.
                # Also, add some dictionary to the graph's molecule nodes.

                # We can trim them later. It's easier to add them initially.
                Graph.add_node(I,
                               mz=I_mz,
                               intensity=I_intensity,
                               type='I')

                Graph.add_edge(M, I)
                for E_interval in Exps[I_mz]:
                    # logger
                    no_exp_per_iso[len(Exps[I_mz])] += 1

                    # Better to collect peaks
                    E_mz, E_intensity = E_interval.data

                    # This is another shit thing:
                    # unused peaks are all experimental with zero degree
                    self.Used_Exps.add(E_interval.data)

                    if E_mz not in Graph:
                        Graph.add_node(E_mz,
                                       intensity=E_intensity,
                                       type='E')
                    Graph.add_edge(I, E_mz)
                    self.stats['E-I No'] += 1
            self.stats['M No'] += 1

        # Logger
        T1 = time()

        # Make one object for collecting these informations:
        # a global singleton stats. It might be possible that logger can do it.
        self.stats['graph construction T'] = T1 - T0

        # Logger
        if self._verbose:
            print('\tFinished graph representation in',
                  self.stats['graph construction T'])
            print('\tIsotopologues number was', self.stats['I No'])
            print('\tEnvelopes number was', self.stats['M No'])
            print()

        return Graph

    def get_problems(self,
                     spectrum,
                     min_prob_per_molecule=.7):
        """Enumerate deconvolution problems.

        Parameters
        ----------
        spectrum : tuple
            The experimental spectrum.
        min_prob_per_molecule :
            The minimal probability an envelope has to scoop
            to be included in the graph.

        """
        Graph = self.represent_as_Graph(spectrum)
        problems = []

        for cc in nx.connected_component_subgraphs(Graph):
            # less M and I, not E
            reduced_cc = trim_unlikely_molecules(cc, min_prob_per_molecule)

            for SG in nx.connected_component_subgraphs(reduced_cc):
                if len(SG) > 1:
                    problems.append(self.__add_G_nodes(SG))
                else:
                    mz, data = SG.nodes(data=True)[0]
                    if data['type'] == 'E':
                        self.Used_Exps.remove((mz, data['intensity']))
                        # these peaks do not get used

        self.stats['total intensity of experimental peaks paired with isotopologues'] = sum(SG.node[G]['intensity'] for SG in problems for G in SG if SG.node[G]['type'] == 'G')

        return problems


    def __add_G_nodes(self, small_graph):
        '''Collect experimental peaks into groups of experimental data G.'''
        T0 = time()

        E_to_remove = []
        Gs_intensity = Counter()
        Gs_min_mz = defaultdict(lambda: inf)
        Gs_max_mz = defaultdict(lambda: 0.0)
        prec = self.mz_precision_digits

        for E in small_graph:
            if small_graph.node[E]['type'] == 'E':
                I_of_G = frozenset(small_graph[E])
                Gs_intensity[I_of_G] += small_graph.node[E]['intensity']
                Gs_min_mz[I_of_G] = min(Gs_min_mz[I_of_G], E)
                Gs_max_mz[I_of_G] = max(Gs_max_mz[I_of_G], E)
                E_to_remove.append(E)

        small_graph.remove_nodes_from(E_to_remove)

        for I_of_G in Gs_intensity:
            G = self.cnts('G')
            small_graph.add_node(G, intensity=Gs_intensity[I_of_G],
                                    type        = 'G',
                                    min_mz      = Gs_min_mz[I_of_G],
                                    max_mz      = Gs_max_mz[I_of_G]     )
            for I in I_of_G:
                small_graph.add_edge(I, G)

        # zero intensity G nodes
        newGnodes = []
        newGIedges= []
        for I in small_graph:
            if small_graph.node[I]['type'] == 'I':
                if len(small_graph[I]) == 1: # only the molecule node in the neighbourhood
                    G = self.cnts('G')
                    mz = small_graph.node[I]['mz']
                    # newGnodes.append( (G,{'intensity':0.0, 'type':'G', 'min_mz':mz-self.mz_precision_digits, 'max_mz':mz+self.mz_precision_digits }) )
                    ## Two E peaks can get together ....
                    newGnodes.append( (G,{'intensity':0.0, 'type':'G', 'min_mz':mz, 'max_mz':mz}) )
                    newGIedges.append( (G,I) )
        small_graph.add_nodes_from(newGnodes)
        small_graph.add_edges_from(newGIedges)
        T1 = time()

        self.G_stats.append(T1-T0)

        return small_graph
