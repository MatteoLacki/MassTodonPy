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

from __future__ import absolute_import, division, print_function
from collections import Counter
from cvxopt import matrix, spmatrix, sparse, spdiag, solvers, setseed
from future.builtins import super
from math import sqrt
import networkx as nx
from random import randint

from MassTodonPy.Deconvolutor.Misc import diag, normalize_rows

class DeconvolutionProblem(nx.Graph):
    """Prepare and solve one deconvolution problem."""
    def __init__(self, data=None, max_times=10,
                 L1_x=0.001, L2_x=0.001,
                 L1_alpha=0.001, L2_alpha=0.001,
                 show_progress=False, maxiters=1000,
                 **attr):
        global solvers
        super().__init__(data, **attr)
        self.count_nodes_and_edges()
        self.max_times = max_times
        self.L1_x = L1_x
        self.L2_x = L2_x
        self.L1_alpha = L1_alpha
        self.L2_alpha = L2_alpha
        self.get_P_q()
        self.get_A_b()
        self.get_G_h()
        self.get_initvals()
        solvers.options['show_progress'] = show_progress
        solvers.options['maxiters'] = maxiters

    def count_nodes_and_edges(self):
        """Add numbers to graph nodes and edges."""
        self.cnts = Counter()
        for N in self:
            self.node[N]['cnt'] = self.cnts[N[0]]
            self.cnts[N[0]] += 1
            if N[0] is not 'I':  # number edges
                for I in self[N]:
                    self[N][I]['cnt'] = self.cnts[N[0]+I[0]]
                    self.cnts[N[0]+I[0]] += 1
        self.cnts['vars'] = self.cnts['GI'] + self.cnts['M']

    def node_iter(self, node_type):
        """Iterate over the nodes of a given type."""
        for N in self:
            if N[0] is node_type:
                yield N

    def get_P_q(self):
        """Get the cost function parameters: matrix P and vector q.

        Notes
        =====
        0.5 <x|P|x> + <q|x> + L1_x * sum x + L2_x * sum x^2 + L1_alpha * sum alpha + L2_alpha * sum alpha^2
        """
        q_list = []
        P_list = []
        for G in self.node_iter('G'):
            intensity = self.node[G]['intensity']
            degree = self.degree(G)
                                            # L1 penalty for x
            q_list.append(matrix(-intensity + self.L1_x,
                                 size=(degree, 1)))
            ones = matrix(1.0, (degree, 1))
                                      # L2 penalty for x
            P_g  = ones * ones.T + diag(self.L2_x, degree)
            P_list.append(P_g)
                           # L1 penalty for alphas
        q_list.append(matrix(self.L1_alpha, (self.cnts['M'], 1)))
        self.q = matrix(q_list)
                         # L2 penalty for alphas
        P_list.append(diag(self.L2_alpha, self.cnts['M']))
        self.P = spdiag(P_list)

    def get_initvals(self):
        '''Initial values of flows.'''
        self.initvals= {}
        self.initvals['x'] = matrix( 0.0, (self.cnts['vars'], 1)  )

    def get_G_h(self):
        '''Prepare for conditions Gx <= h'''
        self.G = diag(-1.0, self.cnts['vars'])
        self.h = matrix(0.0, size=(self.cnts['vars'], 1))

    def get_A_b(self):
        """Prepare for conditions Ax = b."""
        A_x=[]; A_i=[]; A_j=[]
        for M in self.node_iter('M'):
            M_cnt = self.node[M]['cnt']
            for I in self[M]:
                I_cnt = self.node[I]['cnt']
                A_x.append(-self.node[I]['probability'])
                A_i.append(I_cnt)
                A_j.append(M_cnt + self.cnts['GI'])
                for G in self[I]:
                    if not G == M:
                        A_x.append(1.0)
                        A_i.append(I_cnt)
                        A_j.append(self[G][I]['cnt'])
        self.A = spmatrix(A_x, A_i, A_j, size=(self.cnts['I'],
                                               self.cnts['vars']))
        self.A = normalize_rows(self.A)
        self.b = matrix(0.0, (self.cnts['I'], 1))

    def get_L1_error(self):
        return sum(abs(self.node[G]['estimate'] - self.node[G]['intensity'])
                   for G in self.node_iter('G'))

    def get_L2_error(self):
        return sqrt(sum((self.node[G]['estimate']-self.node[G]['intensity'])**2
                        for G in self.node_iter('G')))

    def get_L1_signed_error(self, sign=1.0):
        return sum(max(sign*(self.node[G]['intensity'] -
                             self.node[G]['estimate']),
                       0) for G in self.node_iter('G'))

    def sum_intensities(self):
        return sum(self.node[G]['intensity']
                   for G in self.node_iter('G'))

    def solve(self):
        """Solve the given deconvolution problem."""
        stop = False
        solved = False
        iteration = 1
        while not stop or iteration is self.max_times:
            setseed(randint(0,1000000))
            self.solution = solvers.qp(self.P, self.q, self.G,
                                       self.h, self.A, self.b,
                                       initvals=self.initvals)
            if self.solution['status'] is 'optimal':
                stop = True
                solved = True
            iteration += 1
        if not solved:
            print("Tried {} times and no optimum reached.".format(self.max_times))

        # filling up the graph
        X = self.solution['x']
        self.alphas = []
        for M in self.node_iter('M'):
            self.node[M]['estimate'] = X[self.cnts['GI'] +
                                         self.node[M]['cnt']]
            self.alphas.append(self.node[M])
        for G in self.node_iter('G'):
            self.node[G]['estimate'] = 0.0
            for I in self[G]:
                idx = self[G][I]['cnt']
                self.node[G]['estimate'] += X[idx]
                self[G][I]['estimate'] = X[idx]

    def report(self):
        if self.solution['status'] is 'optimal':
            res = {'alphas': self.alphas,
                   'L1_error': self.get_L1_error(),
                   'L2_error': self.get_L2_error(),
                   'underestimates': self.get_L1_signed_error(sign=1.0),
                   'overestimates': self.get_L1_signed_error(sign=-1.0),
                   'status': self.solution['status']}
        else:
            res = {'status': 'ValueError',
                   'L1_error': self.sum_intensities(),
                   'overestimates': 0.0}
        return res
