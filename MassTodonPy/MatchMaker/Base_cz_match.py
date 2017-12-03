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

from __future__ import absolute_import, division, print_function
from collections import Counter, namedtuple
from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers
from future.builtins import super
import networkx as nx
from networkx import connected_component_subgraphs as connected_components


BaseNode = namedtuple('N', ['type', 'no', 'bp', 'q'])


class Base_cz_match(object):
    """Base class for matching c and z ions' intensities.

    Parameters
    ==========
    results : list
        A list of raw results of **MassTodon.run()**.
    precursor : Precursor
        A precursor for the matching problem.
    minimal_intensity : float
        The minimal intenstity of experimental peaks in the deconvolution graph.

    """
    def __init__(self,
                 results,
                 precursor,
                 minimal_intensity=100.0):
        self.results = results
        self.precursor = precursor
        self.minimal_intensity = minimal_intensity
        self._get_graph()
        self._match()
        self._analyze_intensities()

    def _iter_results(self):
        for res in self.results:
            if res['status'] is not 'ValueError':
                for mol in res['alphas']:
                    estimate = int(mol['estimate'])
                    if estimate > self.minimal_intensity:
                        mol = mol['molecule']
                        yield mol, estimate

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mol_type = molecule.name[0]
        no = int(molecule.name[1:])
        if mol_type is 'c':
            bp = int(molecule.name[1:])
        else:
            bp = len(self.precursor.fasta) - int(molecule.name[1:])
        return BaseNode(mol_type, no, bp, molecule.q)

    def _add_edges_to_graph(self):
        """Add edges between c and z fragments."""
        for C in self.graph:
            if C.type is 'c':
                for Z in self.graph:
                    if Z.type is 'z':
                        if C.bp is Z.bp and C.q + Z.q < self.precursor.q:
                            self.graph.add_edge(C, Z)

    def _get_graph(self):
        """Prepare the matching graph."""
        Q = self.precursor.q
        self.graph = nx.Graph()
        self.not_reacted = 0
        self.ETnoD = 0
        self.PTR = 0
        for mol, estimate in self._iter_results():
            if mol.name is 'precursor':
                q = mol.q
                g = mol.g
                if q is Q and g is 0:
                    self.not_reacted = estimate
                else:
                    self.ETnoD += g * estimate
                    self.PTR += (Q - q - g) * estimate
            else:
                frag = self._get_node(mol)
                if not frag in self.graph:
                    self.graph.add_node(frag, intensity=0)
                    self.graph.node[frag]['intensity'] += estimate
        self._add_edges_to_graph()

    def _optimize(self, G):
        """Match the intensities in a cluster."""
        Q = self.precursor.q
        
        # The intensity of the lavish case:
        #   all fragments lost their cofragments
        lavish = sum((Q - 1 - N.q) * I for N, I in
                     G.nodes.data('intensity'))
        self.lavish += lavish
        
        if len(G) > 1:  # not trivial
            FG = nx.DiGraph()
            FG.add_node('S')  # start
            FG.add_node('T')  # terminus/sink
            for C, C_intensity in G.nodes.data("intensity"):
                if C.type is 'c':
                    FG.add_node(C)
                    FG.add_edge('S', C, capacity=C_intensity)
                    for Z in G[C]:
                        Z_intensity = G.node[Z]['intensity']
                        FG.add_node(Z)
                        FG.add_edge(C, Z)
                        FG.add_edge(Z, 'T', capacity=Z_intensity)
            total_flow, flows = nx.maximum_flow(FG, 'S', 'T')
            self.ETnoD_PTR_fragments += lavish - (Q - 1) * total_flow
            self.ETD += total_flow
            for N in G:
                G.add_edge(N, N)
            for N in flows:
                for M in flows[N]:
                    if N is 'S':  # M is a C fragment
                        G[M][M]['flow'] = G.node[M]['intensity'] - flows[N][M]
                    elif M is 'T':  # N -> Z
                        G[N][N]['flow'] = G.node[N]['intensity'] - flows[N][M]
                    else:  # N -> C and M -> Z
                        G[N][M]['flow'] = flows[N][M]
        else:  # trivial
            self.ETnoD_PTR_fragments += lavish
            N, N_intensity = list(G.nodes.data('intensity'))[0]
            G.add_edge(N, N, flow=N_intensity)

        for N, M, intensity in G.edges.data('flow'):
            self.broken_bond[M.bp] += intensity

    def _analyze_intensities(self):
        total_broken_bonds = sum(v for k, v in self.broken_bond.items())
        self.f_probs = {k: v/total_broken_bonds
                        for k, v in self.broken_bond.items()}


    def _match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities."""
        self.lavish = 0  # total intensity of all lavish pairings
        self.ETD = 0     # total intensity of ETD
        self.ETnoD_PTR_fragments = 0 # total intensity of reactions on fragments
        self.broken_bond = Counter()
        for G in connected_components(self.graph):
            self._optimize(G)
