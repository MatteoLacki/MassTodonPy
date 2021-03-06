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
from    random      import randint
from    cvxopt      import matrix, spmatrix, sparse, spdiag, solvers, setseed
import  traceback

class Error_in_update_scaling(Exception):
    pass


def diag(val, dim):
    return spdiag([spmatrix(val,[0],[0]) for i in xrange(dim)])


def normalize_rows(M):
    '''Divide rows of a matrix by their sums.'''
    for i in xrange(M.size[0]):
        row_hopefully = M[i,:]
        M[i,:] = row_hopefully/sum(abs(row_hopefully))


# def normalize_rows2(M):
#     '''Divide rows of a matrix by their sums.'''
#     for i in xrange(M.size[0]):
#         row_hopefully = M[i,:]
#         M[i,:] = row_hopefully/sum(abs(row_hopefully))


def number_graph(SG):
    '''Add numbers to graph nodes and GI edges and return the counts thereof.'''
    cnts = Counter()
    for N in SG: # ordering some of the SG graph nodes and edges
        Ntype = SG.node[N]['type']
        SG.node[N]['cnt'] = cnts[Ntype]
        cnts[Ntype] += 1
        if Ntype == 'G':
            for I in SG.edge[N]:
                SG.edge[N][I]['cnt'] = cnts['GI']
                cnts['GI'] += 1
    return cnts


def get_P_q(SG, M_No, var_No, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001):
    '''
    Prepare cost function
    0.5 <x|P|x> + <q|x> + L1_x * sum x + L2_x * sum x^2 + L1_alpha * sum alpha + L2_alpha * sum alpha^2
    '''
    q_list = []
    P_list = []
    for G_name in SG:
        if SG.node[G_name]['type']=='G':
            G_intensity = SG.node[G_name]['intensity']
            G_degree    = len(SG[G_name])
            q_list.append( matrix( -G_intensity + L1_x, size=(G_degree,1)) ) # L1 penalty for x
            ones = matrix(1.0, (G_degree,1))
            P_g  = ones * ones.T + diag(L2_x,G_degree) # L2 penalty for x
            P_list.append(P_g)
    q_list.append(matrix(L1_alpha, (M_No,1))) # L1 penalty for alphas
    q_vec = matrix(q_list)
    P_list.append(diag(L2_alpha, M_No))       # L2 penalty for alphas
    P_mat = spdiag(P_list)
    return P_mat, q_vec


def get_initvals(var_No):
    '''Initial values of flows.'''
    initvals= {}
    initvals['x'] = matrix( 0.0, (var_No, 1)  )
    return initvals


def get_G_h(var_No):
    '''Prepare for conditions Gx <= h'''
    G_mat = diag(-1.0, var_No)
    h_vec = matrix(0.0, size=(var_No, 1))
    return G_mat, h_vec


def get_A_b(SG, M_No, I_No, GI_No):
    '''Prepare for conditions Ax = b'''
    A_x=[]; A_i=[]; A_j=[]
    for M in SG:
        if SG.node[M]['type']=='M':
            M_cnt = SG.node[M]['cnt']
            for I in SG[M]:
                i_cnt = SG.node[I]['cnt']
                A_x.append(-SG.node[I]['intensity'])# probability
                A_i.append( i_cnt )
                A_j.append( M_cnt + GI_No )
                for G in SG[I]:
                    if not G == M:
                        A_x.append( 1.0 )
                        A_i.append( i_cnt )
                        A_j.append( SG.edge[G][I]['cnt'] )
    A_mat = spmatrix( A_x, A_i, A_j, size=(I_No, M_No+GI_No ) )
    normalize_rows(A_mat)
    b_vec = matrix( 0.0, (I_No, 1)  )
    return A_mat, b_vec


class Deconvolutor(object):
    '''Class for deconvolving individual Small Graphs.'''
    def __init__(self, SG):
        self.SG = SG
        cnts = number_graph(self.SG)
        self.set_names(cnts)
        solvers.options['show_progress'] = False
        solvers.options['maxiters'] = 1000


    def iSG(self, node_type):
        '''Iterate over all nodes of a given type in the small graph **SG**.

        node_type - either
        '''
        assert node_type in ('G','I','M'), "specified wrong type of node. Was %s. Should be G, I, M" % node_type

        for N in self.SG:
            N = self.SG.node[N]
            if N['type'] == node_type:
                yield N

    def set_names(self, cnts):
        self.var_No = cnts['GI'] + cnts['M']
        self.M_No   = cnts['M']
        self.GI_No  = cnts['GI']
        self.I_No   = cnts['I']
        self.G_No   = cnts['G']

    def run(self, **args):
        '''Perform deconvolution.'''
        raise NotImplementedError

    def get_L1_error(self):
        L1_error = sum(abs(G['estimate']-G['intensity']) for G in self.iSG('G'))
        return float(L1_error)

    def get_sum_of_node_intensities(self):
        node_intensity = sum(G['intensity'] for G in self.iSG('G'))
        return float(node_intensity)

    def get_L2_error(self):
        L2_error = sqrt(sum((G['estimate']-G['intensity'])**2 for G in self.iSG('G')))
        return float(L2_error)

    def get_L1_signed_error(self, sign=1.0):
        L1_sign = sum( max(sign*(G['intensity']-G['estimate']), 0)
            for G in self.iSG('G'))
        return float(L1_sign)

class Deconvolutor_Min_Sum_Squares(Deconvolutor):
    def run(self, L1_x=.001, L2_x=.001, L1_alpha=.001, L2_alpha=.001, verbose=False):
        '''Perform deconvolution that minimizes the mean square error.'''


        if verbose:
            print('Preparing matrices')
            print(L1_x, L2_x, L1_alpha, L2_alpha)

        P, q = get_P_q(self.SG, self.M_No, self.var_No, L1_x, L2_x, L1_alpha, L2_alpha)
        x0   = get_initvals(self.var_No)
        G, h = get_G_h(self.var_No)
        A, b = get_A_b(self.SG, self.M_No, self.I_No, self.GI_No)

        setseed(randint(0,1000000))
        # this is to test from different points
        # apparently this is used by the asynchroneous BLAS library
        # I hate the asynchroneous BLAS library
        try:
            if verbose:
                print('optimizing')
            self.sol = solvers.qp(P, q, G, h, A, b, initvals=x0)
            if verbose:
                print('finished')
            Xopt = self.sol['x']
            #################### reporting results
            alphas = []
            for N_name in self.SG:
                N = self.SG.node[N_name]
                if N['type'] == 'M':
                    N['estimate'] = Xopt[self.GI_No + N['cnt']]
                    alphas.append(N.copy())
                if N['type'] == 'G':
                    N['estimate'] = 0.0
                    for I_name in self.SG[N_name]:
                        NI = self.SG.edge[N_name][I_name]
                        NI['estimate'] = Xopt[NI['cnt']]
                        N['estimate'] += Xopt[NI['cnt']]

            res = { 'alphas':   alphas,
                    'L1_error': self.get_L1_error(),
                    'L2_error': self.get_L2_error(),
                    'underestimates': self.get_L1_signed_error(sign=1.0),
                    'overestimates':  self.get_L1_signed_error(sign=-1.0),
                    'status':   self.sol['status'],
                    'SG':       self.SG
            }
            if verbose:
                res['param']= {'P':P,'q':q,'G':G,'h':h,'A':A,'b':b,'x0':x0}
                res['sol']  = self.sol
        except ValueError as ve:
            print ve
            res = { 'SG':               self.SG,
                    'status':           'ValueError',
                    'L1_error':         self.get_sum_of_node_intensities(),
                    'overestimates':    0.0 }
            res['underestimates'] = res['L1_error']
            if verbose:
                res['param']= {'P':P,'q':q,'G':G,'h':h,'A':A,'b':b,'x0':x0}
                res['exception'] = ve

            print res['L1_error']
            # traceback.print_exc()

        return res

#TODO: make this a valid option to use.
class Deconvolutor_Max_Flow(Deconvolutor):
    def set_names(self, cnts):
        # self.var_No = cnts['GI'] + cnts['M'] + cnts['G'] #TODO: why cnts['G'] appeared?
        self.var_No = cnts['GI'] + cnts['M']
        self.M_No   = cnts['M']
        self.GI_No  = cnts['GI']
        self.I_No   = cnts['I']
        self.G_No   = cnts['G']

    def run(self, lam=10.0, mu=0.0, eps=0.0, s0_val=0.001, verbose=False):
        G_tmp, h_tmp = get_G_h(self.var_No)
        h_eps = matrix(0.0, (self.G_No,1))
        G_x=[]; G_i=[]; G_j=[]
        for G_name in self.SG:
            G = self.SG.node[G_name]
            if G['type'] == 'G':
                g_cnt = G['cnt']
                h_eps[g_cnt] = G['intensity'] + eps
                for I_name in self.SG[G_name]:
                    i_cnt = self.SG.node[I_name]['cnt']
                    G_x.append(1.0)
                    G_i.append(g_cnt)
                    G_j.append(i_cnt)
                G_x.append(-1.0)
                G_i.append(g_cnt)
                G_j.append(self.GI_No + self.M_No + g_cnt)
        G_tmp2 = spmatrix( G_x, G_i, G_j, size=( self.G_No, self.var_No ) )
        G = sparse([G_tmp, G_tmp2])
        h = matrix([h_tmp, h_eps])

        A_tmp, b = get_A_b(self.SG, self.M_No, self.I_No, self.GI_No)
        A_eps = spmatrix([],[],[], (self.I_No, self.G_No))
        A = sparse([[A_tmp],[A_eps]])

        x0 = get_initvals(self.var_No)
        x0['s'] = matrix(s0_val, size=h.size)

        c = matrix([matrix( -1.0,       size=(self.GI_No,1) ),
                    matrix( mu,         size=(self.M_No,1)  ),
                    matrix( (1.0+lam),  size=(self.G_No,1)  )   ])

        setseed(randint(0,1000000))
            # the BLAS library performs calculations in parallel
            # and their order depends upon the seed.
            # It's possible to run the solver with the same input and
            # have different outputs.
            # Some of them have the optimal flag.
            # So, it is a poor man's solution to numeric instability problem.
        self.sol = solvers.conelp(c=c,G=G,h=h,A=A,b=b,primalstart=x0)
        Xopt = self.sol['x']

        alphas = [] # reporting results
        for N_name in self.SG:
            N = self.SG.node[N_name]
            if N['type'] == 'M':
                N['estimate'] = Xopt[self.GI_No + N['cnt']]
                alphas.append(N.copy())
            if N['type'] == 'G':
                N['estimate'] = 0.0
                for I_name in self.SG[N_name]:
                    NI = self.SG.edge[N_name][I_name]
                    assert Xopt[NI['cnt']] > -0.01, "Optimization need to be checked: seriously negative result obtained."
                    estimate = max(Xopt[NI['cnt']], 0.0)
                    NI['estimate'] = Xopt[NI['cnt']]
                    g_cnt = self.GI_No + self.M_No + N['cnt']
                    N['relaxation'] = Xopt[g_cnt]
                    N['estimate'] += Xopt[NI['cnt']]
        res = { 'alphas':   alphas,
                'L1_error': self.get_L1_error(),
                'L2_error': self.get_L2_error(),
                'underestimates': self.get_L1_signed_error(sign=1.0),
                'overestimates':  self.get_L1_signed_error(sign=-1.0),
                'status':sol['status'],
                'SG': self.SG     }
        if verbose:
            res['param']= {'c':c,'G':G,'h':h,'A':A,'b':b,'x0':x0}
            res['sol']  = self.sol
        return res


def deconvolve(SG, args, method):
    deconvolutor = {
        'MSE':      Deconvolutor_Min_Sum_Squares,
        'MaxFlow':  Deconvolutor_Max_Flow
    }[method](SG)
    return deconvolutor.run(**args)
