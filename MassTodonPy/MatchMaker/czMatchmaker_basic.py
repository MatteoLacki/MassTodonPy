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

from    collections  import  Counter
import  networkx     as      nx
from    czMatchmaker import czMatchMaker

class czMatchMakerBasic(czMatchMaker):
    def define_fragment(self, molecule):
        return molecule['molType'], molecule['q']


    def add_edges(self, graph):
        for C, qC in graph: # adding edges between c and z fragments
            if C[0]=='c':
                for Z, qZ in graph:
                    if Z[0]=='z':
                        bpC = int(C[1:])
                        bpZ = len(self.fasta) - int(Z[1:])
                        if bpC==bpZ and qC + qZ < self.Q-1:
                            graph.add_edge((C,qC),(Z,qZ))
        return graph

    def optimize(self, G):
        '''Finds the minimal number of reactions necessary to explain the MassTodon results.

        Uses the max flow algorithm in all but trivial cases.
        '''
        no_edges_reactions_cnt = sum( (self.Q-1-N[1])*G.node[N]['intensity'] for N in G)
        Counts = Counter()
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
            Counts['reactions'] = no_edges_reactions_cnt - (self.Q-1)*flow_val
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
            Counts['reactions'] = no_edges_reactions_cnt
            G.nodes(data=True)
            N = G.nodes()[0]
            G.add_edge( N, N, flow=G.node[N]['intensity'] )
            flow_val= flows = None

        for N in G:
            for M in G[N]:
                if M[0][0]=='z':
                    fragmented_AA = len(self.fasta) - int(M[0][1:])
                else:
                    fragmented_AA = int(M[0][1:])
                Counts[ fragmented_AA ] += G[N][M]['flow']

        if self.verbose:
            return Counts, (flow_val, flows)
        else:
            return Counts
