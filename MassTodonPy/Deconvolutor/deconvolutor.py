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

from    math        import sqrt
from    collections import Counter
from    cvxopt      import matrix, spmatrix, sparse, spdiag, solvers
solvers.options['show_progress'] = False


def normalize_rows(M):
    '''Divide rows of a matrix by their sums.'''
    for i in xrange(M.size[0]):
        row_hopefully = M[i,:]
        M[i,:] = row_hopefully/sum(row_hopefully)


def number_graph(SFG):
    '''Add numbers to graph nodes and GI edges and return the counts thereof.'''
    cnts = Counter()
    for N in SFG: # ordering some of the SFG graph nodes and edges
        Ntype = SFG.node[N]['type']
        SFG.node[N]['cnt'] = cnts[Ntype]
        cnts[Ntype] += 1
        if Ntype == 'G':
            for I in SFG.edge[N]:
                SFG.edge[N][I]['cnt'] = cnts['GI']
                cnts['GI'] += 1
    return cnts


def get_P_q(SFG, M_No, varNo, L2=0.0, spectral_norm=False):
    '''Prepare cost function 0.5 <x|P|x> + <q|x>.'''
    q_list = []
    P_list = []
    for G_name in SFG:
        if SFG.node[G_name]['type']=='G':
            G_intensity = SFG.node[G_name]['intensity']
            G_degree    = len(SFG[G_name])
            q_list.append(matrix( -G_intensity, size=(G_degree,1) ))
            ones = matrix(1.0, (G_degree,1))
            P_list.append(ones * ones.T)

    q_list.append( matrix(0.0, (M_No,1)) )
    q_vec = matrix(q_list)
    P_spectral_norm = 1.0   # spectral normalization
    if spectral_norm:
        P_spectral_norm = max(p_mat.size[0] for  p_mat in P_list) # spec(11')=dim 1
    P_list.append( matrix(0.0, (M_No, M_No)) )
    P_mat = spdiag(P_list)/P_spectral_norm
    q_vec = q_vec/P_spectral_norm
    Id_mat = spmatrix(1.0, xrange(varNo), xrange(varNo))
    if L2:          # L2 regularization
        P_mat = P_mat+L2*max(P_mat)*Id_mat
    return P_mat, q_vec


def get_initvals(varNo):
    '''Initial values of flows.'''
    initvals= {}
    initvals['x'] = matrix( 0.0, ( varNo, 1)  )
    return initvals


def get_G_h(varNo):
    '''Prepare for conditions Gx <= h'''
    G_mat = spmatrix(-1.0, xrange(varNo), xrange(varNo))
    h_vec = matrix(0.0, size=(varNo,1) )
    return G_mat, h_vec


def get_A_b(SFG, varNo, I_No, GI_No):
    '''Prepare for conditions Ax = b'''
    A_x=[]; A_i=[]; A_j=[]
    for M in SFG:
        if SFG.node[M]['type']=='M':
            M_cnt = SFG.node[M]['cnt']
            for I in SFG[M]:
                i_cnt = SFG.node[I]['cnt']
                A_x.append(-SFG.node[I]['intensity'])# probability
                A_i.append( i_cnt )
                A_j.append( M_cnt + GI_No )
                for G in SFG[I]:
                    if not G == M:
                        A_x.append( 1.0 )
                        A_i.append( i_cnt )
                        A_j.append( SFG.edge[G][I]['cnt'] )
    A_mat = spmatrix( A_x, A_i, A_j, size=(I_No, varNo) )
    normalize_rows(A_mat)
    b_vec = matrix( 0.0, (I_No, 1)  )
    return A_mat, b_vec


class Deconvolutor(object):
    '''Class for deconvolving individual Small Graphs.'''
    def __init__(self, SFG):
        self.SFG    = SFG
        self.cnts   = number_graph(self.SFG)
        self.varNo  = self.cnts['GI'] + self.cnts['M']
        self.M_No   = self.cnts['M']
        self.GI_No  = self.cnts['GI']
        self.I_No   = self.cnts['I']

    def run(self):
        '''Perform deconvolution.'''
        raise NotImplementedError


class Deconvolutor_Min_Sum_Squares(Deconvolutor):
    def run(self, L2=0.000001, spectral_norm=False):
        '''Perform deconvolution that minimizes the mean square error.'''

        P, q = get_P_q(self.SFG, self.M_No, self.varNo, L2, spectral_norm)
        x0   = get_initvals(self.varNo)
        G, h = get_G_h(self.varNo)
        A, b = get_A_b(self.SFG, self.varNo, self.I_No, self.GI_No)
        sol  = solvers.qp(P, q, G, h, A, b, initvals=x0)
        Xopt = sol['x']
        # reporting results
        alphas = []
        for N_name in self.SFG:
            N = self.SFG.node[N_name]
            if N['type'] == 'M':
                N['estimate'] = Xopt[self.GI_No + N['cnt']]
                alphas.append(N.copy())
            if N['type'] == 'G':
                for I_name in self.SFG[N_name]:
                    NI = self.SFG.edge[N_name][I_name]
                    NI['estimate'] = Xopt[NI['cnt']]
        # fit error: evaluation of the cost function at the minimizer
        error = 0.0
        for G_name in self.SFG:
            if self.SFG.node[G_name]['type'] == 'G':
                I_intensity = self.SFG.node[G_name]['intensity']
                outflow = 0.0
                for I_name in self.SFG[G_name]:
                    GI = self.SFG.edge[G_name][I_name]
                    outflow += sol['x'][GI['cnt']]
                error += (I_intensity - outflow)**2
        error = sqrt(error)

        return alphas, error, sol['status']

class Deconvolutor_Max_Flow(Deconvolutor):
    def run(self):
        pass

def deconvolve(SFG, method='MSE', **args):
    deconvolutor = {
        'MSE':      Deconvolutor_Min_Sum_Squares,
        'MaxFlow':  Deconvolutor_Max_Flow
    }[method](SFG)
    alphas, error, status = deconvolutor.run(**args)
    return alphas, error, status
