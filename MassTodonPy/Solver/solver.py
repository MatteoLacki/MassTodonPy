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
from    collections import Counter
from    cvxopt import matrix, spmatrix, sparse, spdiag, solvers

def normalize_rows(M):
    '''Divide rows of a matrix by their sums.'''
    for i in xrange(M.size[0]):
        row_hopefully = M[i,:]
        M[i,:] = row_hopefully/sum(row_hopefully)


def prepare_deconvolution(SFG, L2_percent=0.0):
    '''Prepare input for CVXOPT routines (as sparse as possible).'''
    cnts = Counter()
    for N in SFG: # ordering some of the SFG graph nodes and edges
        Ntype = SFG.node[N]['type']
        SFG.node[N]['cnt'] = cnts[Ntype]
        cnts[Ntype] += 1
        if Ntype == 'G':
            for I in SFG.edge[N]:
                SFG.edge[N][I]['cnt'] = cnts['GI']
                cnts['GI'] += 1
    varNo  = cnts['GI']+cnts['M']
    squared_G_intensity = 0.0
    total_G_intensity   = 0.0
    q_list = []
    P_list = []
    for G in SFG:
        if SFG.node[G]['type']=='G':
            G_intensity = SFG.node[G]['intensity']
            squared_G_intensity += G_intensity
            total_G_intensity   += G_intensity
            G_degree = len(SFG[G])
            q_list.append(  matrix( -G_intensity, size=(G_degree,1) )  )
            ones = matrix( 1.0, (G_degree,1))
            P_list.append( 2.0 * ones * ones.T ) # 2 because of .5 x'Px parametrization.
    q_list.append(  matrix( 0.0, ( cnts['M'], 1 ) )  )
    q_vec = matrix(q_list)
    P_spectral_norm = 2.0 * max(p_mat.size[0] for  p_mat in P_list) # spec(11')=dim 1
    P_list.append(  matrix( 0.0, ( cnts['M'], cnts['M'] ) )  )
    P_mat = spdiag(P_list)/P_spectral_norm      # normalization with spectral norm
    q_vec = q_vec/P_spectral_norm               # normalization with spectral norm
    G_mat = spmatrix(-1.0, xrange(varNo), xrange(varNo))
    if L2_percent:
        P_mat = P_mat + L2_percent*max(P_mat)*G_mat # L2 regularization. Sort of.
    h_vec = matrix(0.0, size=(varNo,1) )
    A_x = [];A_i = [];A_j = []
    for M in SFG:
        if SFG.node[M]['type']=='M':
            M_cnt = SFG.node[M]['cnt']
            for I in SFG[M]:
                i_cnt = SFG.node[I]['cnt']
                A_x.append(-SFG.node[I]['intensity'] ) # probability (not intensity)
                A_i.append( i_cnt )
                A_j.append( M_cnt + cnts['GI'] )
                for G in SFG[I]:
                    if not G == M:
                        A_x.append( 1.0 )
                        A_i.append( i_cnt )
                        A_j.append( SFG.edge[G][I]['cnt'] )
    A_mat   = spmatrix( A_x, A_i, A_j, size=( cnts['I'], varNo ) )
    normalize_rows(A_mat)
    b_vec   = matrix( 0.0, ( cnts['I'], 1)  )
    initvals= {}
    initvals['x'] = matrix( 0.0, ( varNo, 1)  )
    return total_G_intensity, cnts, varNo, P_mat, q_vec, G_mat, h_vec, A_mat, b_vec, initvals
