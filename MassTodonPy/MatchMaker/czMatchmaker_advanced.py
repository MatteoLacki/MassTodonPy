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
from    cvxopt      import matrix, spmatrix, sparse, spdiag, solvers
from    czMatchmaker_intermediate import czMatchMakerIntermediate

solvers.options['show_progress'] = False

class czMatchMakerAdvanced(czMatchMakerIntermediate):
    def __init__(self, MassTodonResults, Q, fasta,
                 accept_nonOptimalDeconv  = False,
                 min_acceptEstimIntensity = 100.,
                 verbose=False,
                 L1=0.0, L2=0.01
        ):
        self.MassTodonResults = MassTodonResults
        self.Q = Q
        self.fasta = fasta
        self.min_acceptEstimIntensity   = min_acceptEstimIntensity
        self.accept_nonOptimalDeconv    = accept_nonOptimalDeconv
        self.L1 = L1
        self.L2 = L2
        self.verbose = verbose


    def add_edges(self, graph):
        for Ctype, qC, gC in graph: # adding edges between c and z fragments
                # add self loops
            N = (Ctype, qC, gC)
            graph.add_edge( N, N, ETnoD_PTR_cnt=self.Q-1-qC)
            if Ctype[0]=='c':
                for Ztype, qZ, gZ in graph:
                    if Ztype[0]=='z':
                        bpC = self.get_break_point(Ctype)
                        bpZ = self.get_break_point(Ztype)
                        if bpC==bpZ and qC + qZ + gC + gZ <= self.Q-1:
                            ETnoD_cnt, PTR_cnt = self.etnod_ptr_on_c_z_pairing(qC, gC, qZ, gZ)
                            graph.add_edge( (Ctype,qC,gC), (Ztype,qZ,gZ), ETnoD_PTR_cnt=self.Q-1-qC-qZ, ETnoD=ETnoD_cnt, PTR=PTR_cnt )
        return graph


    def diag(self, val, dim):
        '''Make a sparse identity matrix multiplied by a scalar val.'''
        return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])


    def incidence_matrix(self, G, Jdim, Idim):
        '''Make a sparse incidence matrix of the graph G.'''
        L = spmatrix([], [], [], size=(Jdim,Idim) )
        NodesNo = dict([ (N,i) for i,N in enumerate(G)])
        for j, (N0, N1) in enumerate(G.edges()):
            L[NodesNo[N0],j] = 1
            L[NodesNo[N1],j] = 1
        return L


    def optimize(self, G):
        '''Find regularized max flow.'''
        if len(G) > 1:
            bP= self.get_break_point(next(G.nodes_iter())[0])
            J = matrix([ float(G.node[N]['intensity']) for N in G ])
            Jdim = J.size[0]
            R = matrix([ D['ETnoD_PTR_cnt'] for N0, N1, D in G.edges(data=True) ])
            R = R + self.L1*1.0
            Idim = R.size[0]
            P = self.diag(self.L2*2, Idim)
            L = self.incidence_matrix(G, Jdim, Idim)
            Gmat = self.diag(-1.0, Idim)
            h = matrix(0.0, size=(Idim, 1))
            S = solvers.qp(P=P, q=R, G=Gmat, h=h, A=L, b=J)
            I = S['x']
            TotalFrags = sum(I)
            TotalPTR   = 0.0
            TotalETnoD = 0.0
            for i, (N0, N1) in enumerate(G.edges_iter()):
                if N0 != N1:
                    TotalETnoD  += I[i]*G.edge[N0][N1]['ETnoD']
                    TotalPTR    += I[i]*G.edge[N0][N1]['PTR']
        else:
            (nType, nQ, nG), Data = G.nodes(data=True)[0]
            I   = Data['intensity']
            bP  = self.get_break_point(nType)
            TotalPTR   = 0
            TotalETnoD = 0
            TotalFrags = I
            S = {'status':'trivial', 'x':I}

        res = Counter({'ETnoD_frag':TotalETnoD, 'PTR_frag':TotalPTR, bP: TotalFrags})
        if self.verbose:
            return res, S
        else:
            return res
