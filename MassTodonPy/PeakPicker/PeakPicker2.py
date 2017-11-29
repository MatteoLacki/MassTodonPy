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
import networkx as nx
from networkx import connected_component_subgraphs as components

from MassTodonPy.Data.Constants import infinity

def get_deconvolution_problems(molecules,
                               spectrum,
                               mz_tol=.05,
                               min_prob_per_molecule=.7,
                               joint_probability=.999,
                               mz_precision=3):
    """Get the sequence of deconvolution problems."""
    deconvolution_graph = nx.Graph()

    M_cnt = 0
    I_cnt = 0
    E_cnt = 0

    for mol in molecules:
        M = 'M' + str(M_cnt)
        M_cnt += 1

        # this graph might be included in the deconvolution graph
        graph = nx.Graph()
        graph.add_node(M, molecule=mol)
        mol.graph_tag = M

        for mz_I, prob in mol.isotopologues(joint_probability,
                                            mz_precision):
            I = 'I' + str(I_cnt)
            I_cnt += 1
            graph.add_node(I, mz=mz_I, probability=prob)
            graph.add_edge(M, I)

            for mz_E, intensity in spectrum[mz_I - mz_tol, mz_I + mz_tol]:
                E = 'E' + str(E_cnt)
                E_cnt += 1
                graph.add_node(E, mz=mz_E, intensity=intensity)
                graph.add_edge(I, E)

        # total probability of isotopologues next to experimental peaks
        total_prob = sum(d['probability']
                         for n, d in graph.nodes(data=True)
                         if n[0] is 'I' and graph.degree[n] > 1)

        if total_prob >= min_prob_per_molecule: # update Deconvolution Graph
            deconvolution_graph.add_nodes_from(graph.nodes(data=True))
            deconvolution_graph.add_edges_from(graph.edges(data=True))

    for problem in components(deconvolution_graph):
        yield add_G_nodes(problem)


def add_G_nodes(graph):
    """Aggregate experimental nodes and remove experimental nodes."""
    # key - a frozenset of neighbouring isotopologues
    G_intensity = Counter()
    G_min_mz = defaultdict(lambda: infinity)
    G_max_mz = defaultdict(lambda: 0.0)

    for E, E_data in graph.nodes(data=True):
        if E[0] is 'E':
            isotopologues = frozenset(graph[E]) # unmutable!
            G_intensity[isotopologues] += E_data['intensity']
            G_min_mz[isotopologues] = min(G_min_mz[isotopologues], E_data['mz'])
            G_max_mz[isotopologues] = max(G_max_mz[isotopologues], E_data['mz'])

    # removing experimental peaks entirely
    graph.remove_nodes_from([E for E in graph.nodes() if E[0] is 'E'])

    # add G nodes with positive intensity
    G_cnt = 0
    for isotopologues in G_intensity:
        G = 'G' + str(G_cnt)
        G_cnt += 1
        graph.add_node(G,
                       intensity=G_intensity[isotopologues],
                       min_mz=G_min_mz[isotopologues],
                       max_mz=G_max_mz[isotopologues])
        for I in isotopologues:
            graph.add_edge(I, G)

    # add G nodes with zero intensity: without the experimental support
    new_nodes = []  # to avoid changing dict size during iteration
    new_edges = []  # explicitly construct lists of nodes and edges
    for I in graph.nodes():
        if I[0] is 'I':
            if graph.degree[I] is 1: # no experimental support
                G = 'G' + str(G_cnt)
                G_cnt += 1
                new_nodes.append((G, {'mz': graph.node[I]['mz'],
                                      'intensity': 0.0}))
                new_edges.append((I, G))
    graph.add_nodes_from(new_nodes)
    graph.add_edges_from(new_edges)
    return graph
