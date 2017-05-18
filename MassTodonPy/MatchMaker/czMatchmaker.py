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

from    collections     import  Counter
import  networkx        as      nx

class czMatchMaker(object):
    '''Virtual class of all matchmakers.'''

    def __init__(self, MassTodonResults, Q, fasta, accept_nonOptimalDeconv = False, min_acceptEstimIntensity = 100., verbose=False):
        self.MassTodonResults = MassTodonResults
        self.Q = Q
        self.fasta = fasta
        self.verbose = verbose
        self.min_acceptEstimIntensity = min_acceptEstimIntensity
        self.accept_nonOptimalDeconv = accept_nonOptimalDeconv

    def define_fragment(self, molecule):
        '''Defines what should be considered a node in the c-z matching graphs.'''
        raise NotImplementedError

    def add_edges(self, graph):
        '''Defines what should be considered a node in the c-z matching graphs.'''
        raise NotImplementedError

    def get_graph_analyze_precursor(self):
        '''Generate the graph of pairings, find its connected components, find the number of PTR and ETnoD reactions on precursors.'''
        unreacted_precursors = ETnoDs_on_precursors = PTRs_on_precursors = 0.0
        graph = nx.Graph()
        Q = self.Q
        for res in self.MassTodonResults:
            if self.accept_nonOptimalDeconv or res['status']=='optimal': #TODO what to do otherwise? Nothing for now.
                for mol in res['alphas']:
                    estimate = mol['estimate']
                    if estimate > self.min_acceptEstimIntensity:
                        if mol['molType']=='precursor':
                            q, g = mol['q'], mol['g']
                            if q==Q and g==0:
                                unreacted_precursors = estimate
                            else:
                                ETnoDs_on_precursors += g * estimate
                                PTRs_on_precursors   += (Q-q-g) * estimate
                        else:
                            frag = self.define_fragment(mol)
                            if not frag in graph:
                                graph.add_node( frag, intensity=0 )
                            graph.node[frag]['intensity'] += int(estimate)
        ETnoDs_on_precursors = int(ETnoDs_on_precursors)
        PTRs_on_precursors   = int(PTRs_on_precursors)
        graph = self.add_edges(graph)
        return graph, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors


    def optimize(self):
        raise NotImplementedError


    def analyze_counts(self, Counts):
        raise NotImplementedError


    def pair(self, verbose=False):
        '''Pair molecules minimizing the number of reactions and calculate the resulting probabilities.'''
        Counts = Counter()
        graph, ETnoDs_on_precursors, PTRs_on_precursors, unreacted_precursors = self.get_graph_analyze_precursor()

        Counts['ETnoD_precursor'] = ETnoDs_on_precursors
        Counts['PTR_precursor']   = PTRs_on_precursors

        OptimInfo = []
        for G in nx.connected_component_subgraphs(graph):
            if verbose:
                cnt, info = self.optimize(G)
                Counts += cnt
                OptimInfo.append(info)
            else:
                Counts += self.optimize(G)

        # self.analyze_counts(Counts)

        return Counts
