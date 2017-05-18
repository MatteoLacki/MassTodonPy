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
import  networkx        as     nx
from    collections     import Counter
from    czMatchmaker    import czMatchMaker


class czMatchMakerIntermediate(czMatchMaker):
    def define_fragment(self, molecule):
        molG = molecule['g']
        molQ = molecule['q']
        if molG == - 1:     # HTR product
            molG += 1
        if molG + molQ == self.Q:# HTR product
            molG -= 1
        return molecule['molType'], molQ, molG

    def etnod_ptr_on_c_z_pairing(self, q0, g0, q1, g1):
        '''Get the number of ETnoD and PTR reactions on a regular edge.'''
        Netnod  = g0 + g1
        Nptr    = self.Q - 1 - g0 - g1 - q0 - q1
        return Netnod, Nptr

    def get_break_point(self, nType ):
        '''Get the amino acid number that was cleft.'''
        if nType[0] == 'c':
            bP = int(nType[1:])
        else:
            bP = len(self.fasta) - int(nType[1:])
        return bP

    def add_edges(self, graph):
        for Ctype, qC, gC in graph: # adding edges between c and z fragments
                if Ctype[0]=='c':
                    for Ztype, qZ, gZ in graph:
                        if Ztype[0]=='z':
                            bpC = self.get_break_point(Ctype)
                            bpZ = self.get_break_point(Ztype)
                            if bpC==bpZ and qC + qZ + gC + gZ <= self.Q-1:
                                ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing( qC, gC, qZ, gZ )
                                graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD=ETnoD_cnt, PTR=PTR_cnt )
        return graph

    def optimize(self, G):
        if len(G) > 1:
            Jsum = sum( G.node[N]['intensity'] for N in G)
            bP = self.get_break_point( next(G.nodes_iter())[0])
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
            TotalFrags  = Jsum - flow_val
            TotalETnoD  = 0
            TotalPTR    = 0
            for C, Z in G.edges_iter():
                if C[0][0]=='z':
                    C,Z = Z,C
                TotalETnoD  += flows[C][Z]*G.edge[C][Z]['ETnoD']
                TotalPTR    += flows[C][Z]*G.edge[C][Z]['PTR']
        else:
            (nType, nQ, nG), Data =  G.nodes(data=True)[0]
            I   = Data['intensity']
            bP  = self.get_break_point(nType)
            TotalPTR   = 0
            TotalETnoD = 0
            TotalFrags = I
            flow_val = flows = None

        Counts = Counter({'ETnoD_frag':TotalETnoD, 'PTR_frag':TotalPTR, bP: TotalFrags})
        if self.verbose:
            return Counts, (flow_val, flows)
        else:
            return Counts

    def get_etnod_ptr_probs(self, Counts, Probs):
        Probs = self.get_probs(Counts,Probs,'ETnoD_precursor','PTR_precursor')
        Probs = self.get_probs(Counts,Probs,'ETnoD','PTR')
        Probs = self.get_probs(Counts,Probs,'ETnoD_frag','PTR_frag')
        return Probs
