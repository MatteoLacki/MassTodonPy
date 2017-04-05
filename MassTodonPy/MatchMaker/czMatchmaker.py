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

def minimal_cost(G, Q):
    '''Finds the minimal number of reactions necessary to explain the MassTodon results.

    Uses the max flow algorithm in all but trivial cases.
    '''
        # The number of reactions other than fragmentation
        # that would result from not having any edges between C and Z ions.
    no_edges_reactions_cnt = sum( (Q-1-N[1])*G.node[N]['intensity'] for N in G)

    if len(G)>1:
        FG = nx.DiGraph()
        FG.add_node('S') # start
        FG.add_node('T') # terminus/sink
        for C in G:
            if C[0][0]=='c':
                Cintensity = G.node[C]['intensity']
                FG.add_node(C)
                FG.add_edge( 'S', C, capacity=Cintensity )
                for Z in G[C]:
                    Zintensity = G.node[Z]['intensity']
                    FG.add_node(Z)
                    FG.add_edge(C,Z)
                    FG.add_edge( Z, 'T', capacity=Zintensity )
        flow_val, flows = nx.maximum_flow(FG,'S','T')
        min_cost = no_edges_reactions_cnt - (Q-1)*flow_val
        for N in G:
            G.add_edge(N,N)
        for N in flows:
            for M in flows[N]:
                if N=='S': # M is a C fragment
                    G.edge[M][M]['flow'] = G.node[M]['intensity']-flows[N][M]
                elif M=='T': # N is a Z fragment
                    G.edge[N][N]['flow'] = G.node[N]['intensity']-flows[N][M]
                else: # N is a C and M a Z fragment
                    G.edge[N][M]['flow'] = flows[N][M]
    else:
        min_cost = no_edges_reactions_cnt
        G.nodes(data=True)
        N = G.nodes()[0]
        G.add_edge( N, N, flow=G.node[N]['intensity'] )
    return min_cost, G


def reaction_analist_basic(MassTodonResults, fasta, Q):
    '''Estimate probabilities of reactions out of the MassTodon results.

    Divide the molecules using a c-z molecules that can be obtained by minimizing the total number of reactions needed to express the MassTodon output.'''

    no_reactions = ETnoD_cnt = PTR_cnt = 0.0
    L   = len(fasta)
    BFG = nx.Graph()
    minimal_estimated_intensity = 100.

    for mols, error, status in MassTodonResults:
        if status=='optimal': #TODO what to do otherwise?
            for mol in mols:
                if mol['estimate'] > minimal_estimated_intensity: # a work-around the stupidity of the optimization methods
                    if mol['molType']=='precursor':
                        if mol['q']==Q and mol['g']==0:
                            no_reactions = mol['estimate']
                        else:
                            ETnoD_cnt  += mol['g'] * mol['estimate']
                            PTR_cnt    += (Q-mol['q']-mol['g']) * mol['estimate']
                    else:
                        frag = ( mol['molType'], mol['q'] )
                        if not frag in BFG:
                            BFG.add_node( frag, intensity=0 ) # intensities will be integers
                        BFG.node[frag]['intensity'] += int(mol['estimate']) # convert the intensities to ints

    reactions_on_precursors = ETnoD_cnt + PTR_cnt
    prob_PTR   = float(PTR_cnt)/reactions_on_precursors
    prob_ETnoD = 1.0 - prob_PTR

    for C, qC in BFG: # adding edges between c and z fragments
        if C[0]=='c':
            for Z, qZ in BFG:
                if Z[0]=='z':
                    bpC = int(C[1:])
                    bpZ = L - int(Z[1:])
                    if bpC==bpZ and qC + qZ < Q-1:
                        BFG.add_edge((C,qC),(Z,qZ))


    fragmentations_no_aas = Counter()
    reactions_on_frags_other_than_fragmentation = 0

    for cc in nx.connected_component_subgraphs(BFG):
        reactionNo, G = minimal_cost(cc, Q)
        reactions_on_frags_other_than_fragmentation += reactionNo
        for N in G:
            for M in G[N]:
                if M[0][0]=='z':
                    fragmented_AA = L-int(M[0][1:])
                else:
                    fragmented_AA = int(M[0][1:])
                fragmentations_no_aas[ fragmented_AA ] += G[N][M]['flow']

    fragmentations_no_total = sum(fragmentations_no_aas.values())

    prob_no_reaction = float(no_reactions)/ (no_reactions+reactions_on_precursors+reactions_on_frags_other_than_fragmentation+fragmentations_no_total)
    prob_reaction = 1.0 - prob_no_reaction

    prob_fragmentation = float(fragmentations_no_total)/( fragmentations_no_total+reactions_on_frags_other_than_fragmentation+reactions_on_precursors )
    prob_no_fragmentation = 1.0 - prob_fragmentation

    probs_fragmentation_on_aas = [ float(fragmentations_no_aas[i])/fragmentations_no_total for i in xrange(len(fasta)+1)]

    results = { 'prob_PTR'          :   prob_PTR,
                'prob_ETnoD'        :   prob_ETnoD,
                'prob_no_reaction'  :   prob_no_reaction,
                'prob_reaction'     :   prob_reaction,
                'prob_fragmentation':   prob_fragmentation,
                'prob_no_fragmentation': prob_no_fragmentation,
                'probs_fragmentation_on_aas':   probs_fragmentation_on_aas }

    return results
