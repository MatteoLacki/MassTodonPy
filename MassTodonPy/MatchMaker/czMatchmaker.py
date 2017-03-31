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

def update_nators(mol, Q, nominator=0.0, denominator=0.0):
    ETnoDs  = mol['g']
    PTRs    = Q-mol['q']-mol['g']
    nominator   += PTRs*mol['estimate']
    denominator += (ETnoDs+PTRs)*mol['estimate']
    return nominator, denominator

def min_cost_flow(G, Q, verbose=False):
    '''Finds the minimal number of reactions necessary to explain the MassTodon results.

    Uses the min_cost_flow algorithm.'''
    FG = nx.DiGraph() # flow graph
    FG.add_node('T', demand=0) # we gonna play with integers to eliminate real numbers' imprecisions.
    for N in G:
        intensity = int(G.node[N]['intensity'])
        FG.add_node( N, demand = -intensity )
        FG.node['T']['demand'] += intensity
    for N, M in G.edges_iter():
        if N == M:
            FG.add_edge( N, 'T', weight=Q-1-N[1] )
        else:
            FG.add_node((N,M))
            FG.add_edge( N, (N,M))
            FG.add_edge( M, (N,M))
            FG.add_edge( (N,M), 'T', weight=Q-1-N[1]-M[1])
    res = nx.min_cost_flow_cost(FG)
    if verbose:
        res = (res, nx.min_cost_flow(FG))
    return res

def reaction_analist_basic(MassTodonResults, fasta, Q):
    '''Estimate probabilities of reactions out of the MassTodon results.

    Divide the molecules using a c-z molecules that can be obtained by minimizing the total number of reactions needed to express the MassTodon output.'''
    no_reactions = denominator = nominator = 0.0
    L = len(fasta)
    IDG = nx.Graph() # the intensity division graph
    minimal_estimated_intensity = 100.

    for mols, error, status in MassTodonResults:
        if status=='optimal': #TODO what to do otherwise?
            for mol in mols:
                if mol['estimate'] > minimal_estimated_intensity:
                    if mol['molType']=='precursor':
                        if mol['q']==Q and mol['g']==0:
                            no_reactions = mol['estimate']
                        else:
                            nominator, denominator = update_nators(mol, Q, nominator, denominator)
                    else:
                        frag = (mol['molType'],mol['q'])
                        IDG.add_node(frag, intensity=mol['estimate'] )
                        IDG.add_edge(frag,frag)

    prob_PTR = nominator/denominator
    prob_ETnoD = 1.0 - prob_PTR
    reactions_on_precursors = denominator

    for C, qC in IDG: # adding edges between c and z fragments
        if C[0]=='c':
            for Z, qZ in IDG:
                if Z[0]=='z':
                    bpC = int(C[1:])
                    bpZ = L - int(Z[1:])
                    if bpC==bpZ and qC + qZ < Q-1:
                        IDG.add_edge((C,qC),(Z,qZ))

    res = [ min_cost_flow(cc, Q, verbose=True) for cc in nx.connected_component_subgraphs(IDG) ]

    fragmentations_no_aas = Counter()
    reactions_on_frags_other_than_fragmentation = 0
    for reactionNo, flow in res:
        for N in flow:
            for M in flow[N]:
                if M=='T':
                    if N[0][0] == 'c':
                        breakPoint = int(N[0][1:])
                    if N[0][0] == 'z':
                        breakPoint = L-int(N[0][1:])
                    fragmentations_no_aas[ breakPoint ] += flow[N][M]
        reactions_on_frags_other_than_fragmentation += reactionNo

    fragmentations_no_total = sum(fragmentations_no_aas.values())

    prob_no_reaction = float(no_reactions)/ (no_reactions+reactions_on_precursors+reactions_on_frags_other_than_fragmentation+fragmentations_no_total)
    prob_reaction = 1.0 - prob_no_reaction

    prob_fragmentation = float(fragmentations_no_total)/( fragmentations_no_total+reactions_on_frags_other_than_fragmentation+reactions_on_precursors )
    prob_no_fragmentation = 1.0 - prob_fragmentation

    probs_fragmentation_on_aas = [ float(fragmentations_no_aas[i])/fragmentations_no_total for i in xrange(len(fasta))]

    results = { 'prob_PTR'          :   prob_PTR,
                'prob_ETnoD'        :   prob_ETnoD,
                'prob_no_reaction'  :   prob_no_reaction,
                'prob_reaction'     :   prob_reaction,
                'prob_fragmentation':   prob_fragmentation,
                'prob_no_fragmentation': prob_no_fragmentation,
                'probs_fragmentation_on_aas':   probs_fragmentation_on_aas }
    return results
