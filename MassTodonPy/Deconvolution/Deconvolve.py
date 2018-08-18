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

from __future__  import absolute_import, division, print_function
from collections import Counter
from collections import defaultdict
import cvxopt
import imp
import networkx as      nx
from networkx   import  connected_component_subgraphs

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Deconvolution.DeconvolutionProblem import DeconvolutionProblem
# from MassTodonPy.Deconvolution.Wanda_Ciacho_DeconvolutionProblem import GaussianDeconvolutionProblem
# from MassTodonPy.Misc.wrapper import cvxopt_wrapper


#TODO: build it in parts! Change it to supplier of 
# connected components?
def build_deconvolution_graph(molecules,
                              spectrum,
                              mz_tol=.05,
                              min_prob_per_molecule=.7,
                              method="Matteo",
                              _merge_sister_Is=True,
                              _verbose=False,
                              **isospec_args):
    """Build a deconvolution graph.

    It will later on get trimmed to reasonable sizes.

    Parameters
    ==========
    molecules : iterable
        An iterable of the Molecule objects.
    spectrum : Spectrum
        An instance of the Spectrum class.
    mz_tol : float or function
        The tolerance in the m/z axis.
        Ultimately turned into a function mz_tol(mz),
        reporting the tolerance interval - a tuple '(mz_L, mz_R)',
        where mz_L = mz - mz_tol and mz_R = mz + mz_tol.
        Accepts user-provided functions too: in that case, the header should be
        'def mz_tol(mz):' and it should return a tuple '(mz_L, mz_R)'.
    min_prob_per_molecule : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.
    method: string
        Input 'Matteo' for MassTodon paper deconvolution.
        Input 'Ciacho_Wanda' for gaussian kernel deconvolution (not yet implemented)
    Returns
    =======
    deconvolution_graph: networkx.Graph
        A graph build of molecule: nodes M, isotopologues I, experimental peaks E.
    """
    I_cnt = 0
    deconvolution_graph = nx.Graph()
    if isinstance(mz_tol, float):
        _mz_tol = mz_tol
        def mz_tol(mz):
            """Calculates the tolerance interval around m/z."""
            return (mz - _mz_tol, mz + _mz_tol)

    # build up the deconvolution graph
    for M_cnt, mol in enumerate(molecules):
        M = 'M' + str(M_cnt)
        # mol_graph represents how a molecule links to the spectrum
        # and might be included in the deconvolution graph
        mol_graph = nx.Graph()
        mol_graph.add_node(M, molecule=mol)
        # this is the reverse mapping between molecules and the deconvolution_graph
        mol.graph_tag = M

        # if mol.name == 'precursor' and mol.q in (15, 16, 17, 18, 19):
        #     print(mol.name)
        # add isotopologues
        for mz_I, prob in mol.isotopologues(**isospec_args):
            I = 'I' + str(I_cnt)
            I_cnt += 1
            mol_graph.add_node(I, mz=mz_I, probability=prob)
            mol_graph.add_edge(M, I)
            mz_L, mz_R = mz_tol(mz_I)  # tolerance interval
            # if mol.name == 'precursor' and mol.q in (15, 16, 17, 18, 19):
            #     print(mz_L, mz_R)
            # link with real peaks ---------------> show indices ↓ ?
            for E_cnt, mz_E, intensity in spectrum[mz_L, mz_R, True]:
                # get E_cnt additionally to m/z and intensity <--|
                # if mol.name == 'precursor' and mol.q in (15, 16, 17, 18, 19):
                #     print(mz_E, '-->', intensity)
                E = 'E' + str(E_cnt)
                mol_graph.add_node(E, mz=mz_E, intensity=intensity)
                mol_graph.add_edge(I, E)

        # total probability of isotopologues around experimental peaks
        total_prob = sum(d['probability']
                         for n, d in mol_graph.nodes(data=True)
                         if n[0] == 'I' and mol_graph.degree[n] > 1)

        # if mol.name == 'precursor' and mol.q in (15, 16, 17, 18, 19):
        #     print(len(mol_graph))
        #     print(total_prob)
        #     print()


        # plant the mol_graph into the deconvolution_graph?
        if total_prob >= min_prob_per_molecule:
            #TODO: add another intensity-based criterion here.
            #TODO: basically, check how much of a substance could there be,
            #      if the was the only possible source of these ions.
            #      and accept if it is more than some number.
            #      This way silly solutions should be eliminated.
            if method == 'Matteo' and _merge_sister_Is:
                # Glue Is sharing common Es.
                _glue_sister_isotopologues(mol_graph)

            # if mol.name == 'precursor':
            #     print(mol.name, mol.q, mol.g , total_prob)

            deconvolution_graph.add_nodes_from(mol_graph.nodes(data=True))
            deconvolution_graph.add_edges_from(mol_graph.edges(data=True))

    return deconvolution_graph


def deconvolve(deconvolution_graph,
               method="Matteo",
              _verbose=False,
             **deconvolution_problem_kwds):
    """Get the sequence of deconvolution problems.

    Parameters
    ==========
    deconvolution_graph: networkx.Graph
        A graph build of molecule: nodes M, isotopologues I, experimental peaks E.
    method: string
        Input 'Matteo' for MassTodon paper deconvolution.
        Input 'Ciacho_Wanda' for gaussian kernel deconvolution (not yet implemented)
    deconvolution_problem_kwds : dict
        Arguments for the deconvolution solver.
    """
    if method == 'Matteo':
        if _verbose:
            print(Counter(n[0] for n in deconvolution_graph.nodes))
        _add_G_remove_E(deconvolution_graph)
        graphs = connected_component_subgraphs(deconvolution_graph)
        for g in graphs:
            iters = 0
            solved = False
            max_iters = 10
            problem = DeconvolutionProblem(g, 
                                         **deconvolution_problem_kwds)
            yield problem
    elif method == 'Ciacho_Wanda':
        graphs = connected_component_subgraphs(deconvolution_graph)
        for g in graphs:
            problem = GaussianDeconvolutionProblem(g,
                                                 **deconvolution_problem_kwds)
            problem.solve()
            yield problem
    else:
        raise NotImplemented("Choose 'Matteo' or 'Ciacho_Wanda'.")


def _glue_sister_isotopologues(deconvolution_graph):
    """Merge I nodes that share parents.

    The final deconvolution graph has I with the id of the smallest I in the group
    of Is that share common E nodes. I know it could be made more clear, but trust me,
    I know what I'm doing (https://cdn-images-1.medium.com/max/455/1*snTXFElFuQLSFDnvZKJ6IA.png).

    Parameters
    ==========
    deconvolution_graph: networkx.Graph
        A graph build of molecule: nodes M, isotopologues I, experimental peaks E.
    """
    visited = {}
    I_2_delete = []
    for I in deconvolution_graph:
        if I[0] == 'I':
            Gs = frozenset(G for G in deconvolution_graph[I] if G[0] != 'M')
            if Gs:
                if not Gs in visited:
                    visited[Gs] = I
                    deconvolution_graph.node[I]['min_mz'] = deconvolution_graph.node[I]['mz']
                    deconvolution_graph.node[I]['max_mz'] = deconvolution_graph.node[I]['mz']
                    del deconvolution_graph.node[I]['mz']
                else:
                    _I = visited[Gs]
                    deconvolution_graph.node[_I]['min_mz'] = min(deconvolution_graph.node[_I]['min_mz'],
                                                       deconvolution_graph.node[I]['mz'])
                    deconvolution_graph.node[_I]['max_mz'] = max(deconvolution_graph.node[_I]['max_mz'],
                                                       deconvolution_graph.node[I]['mz'])
                    deconvolution_graph.node[_I]['probability'] += deconvolution_graph.node[I]['probability']
                    I_2_delete.append(I)
    deconvolution_graph.remove_nodes_from(I_2_delete)


def _add_G_remove_E(deconvolution_graph):
    """Merge experimental nodes E into experimental group nodes G.

    Parameters
    ==========
    deconvolution_graph: networkx.Graph
        The graph consisting of a molecule node M,
        its isotopologues I, and their experimental peaks E.
    """
    G_intensity = Counter()
    G_min_mz = defaultdict(lambda: infinity)
    G_max_mz = defaultdict(lambda: 0.0)
    for E, E_data in deconvolution_graph.nodes(data=True):
        if E[0] == 'E':
            isotopologues = frozenset(deconvolution_graph[E]) # unmutable!
            G_intensity[isotopologues] += E_data['intensity']
            G_min_mz[isotopologues] = min(G_min_mz[isotopologues], E_data['mz'])
            G_max_mz[isotopologues] = max(G_max_mz[isotopologues], E_data['mz'])

    # removing experimental peaks 'E'
    deconvolution_graph.remove_nodes_from([E for E in deconvolution_graph.nodes() if E[0] == 'E'])

    # add G nodes with positive intensity
    G_cnt = 0
    for isotopologues in G_intensity:
        G = 'G' + str(G_cnt)
        G_cnt += 1
        deconvolution_graph.add_node(G,
                                     intensity=G_intensity[isotopologues],
                                     min_mz=G_min_mz[isotopologues],
                                     max_mz=G_max_mz[isotopologues])
        for I in isotopologues:
            deconvolution_graph.add_edge(I, G)

    # add G nodes with zero intensity: without the experimental support
    new_nodes = []  # to avoid "changing dict size" during iteration
    new_edges = []  # explicitly construct lists of nodes and edges
    for I in deconvolution_graph.nodes():
        if I[0] == 'I' and deconvolution_graph.degree[I] == 1: # no experimental support
            G = 'G' + str(G_cnt)
            G_cnt += 1
            new_nodes.append((G, {'mz': deconvolution_graph.node[I]['mz'],
                                  'intensity': 0.0}))
            new_edges.append((I, G))
    deconvolution_graph.add_nodes_from(new_nodes)
    deconvolution_graph.add_edges_from(new_edges)
