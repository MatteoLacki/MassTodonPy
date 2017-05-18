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
                                unreacted_precursors = int(estimate)
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

    def get_probs(self, Counts, Probs, tag1, tag2, name1=None, name2=None):
        if not name1:
            name1 = tag1
        if not name2:
            name2 = tag2
        total = Counts[tag1]+Counts[tag2]
        if total > 0.0:
            Probs[name1] = float(Counts[tag1])/total
            Probs[name2] = 1.0 - Probs[name1]
        return Probs

    def get_etnod_ptr_probs(self, Counts, Probs):
        Probs = self.get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = self.get_probs(Counts,Probs,'ETnoD','PTR')
        return Probs

    def analyze_counts(self, Counts):
        Counts['total_frags'] = float(sum(Counts[k] for k in Counts if isinstance(k,(int,long))))
        Probs = Counter()
        if Counts['total_frags'] > 0.0:
            for k in Counts:
                if isinstance(k, (int,long)):
                    Probs[k] = float(Counts[k])/Counts['total_frags']
        Counts['total_reactions'] = sum(Counts[k] for k in Counts if k != 'unreacted_precursors')
        Probs = self.get_probs( Counts, Probs, 'unreacted_precursors', 'total_reactions', 'anion_did_not_approach_cation', 'anion_approached_cation' )
        if Counts['total_reactions'] > 0.0:
            Probs['fragmentation'] = float(Counts['total_frags'])/Counts['total_reactions']
            Probs['no fragmentation'] = 1.0 - Probs['fragmentation']
        Counts['ETnoD'] = Counts['ETnoD_frag'] + Counts['ETnoD_precursor']
        Counts['PTR']   = Counts['PTR_frag'] + Counts['PTR_precursor']
        Probs = self.get_etnod_ptr_probs(Counts,Probs)
        return Probs, Counts


    def pair(self, verbose=False):
        '''Pair molecules minimizing the number of reactions and calculate the resulting probabilities.'''
        Counts = Counter()

        graph, Counts['ETnoD_precursor'], Counts['PTR_precursor'], Counts['unreacted_precursors'] = self.get_graph_analyze_precursor()

        OptimInfo = []
        for G in nx.connected_component_subgraphs(graph):
            if self.verbose:
                cnt, info = self.optimize(G)
                Counts += cnt
                OptimInfo.append(info)
            else:
                Counts += self.optimize(G)

        Probs, Counts = self.analyze_counts(Counts)

        if self.verbose:
            return Probs, Counts, OptimInfo
        else:
            return Probs, Counts
